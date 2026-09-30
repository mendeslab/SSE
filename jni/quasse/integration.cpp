// Native T/X calculations. Java chooses time steps and handles grids and normalization.
#include "SSE_QuaSSENativeIntegrator.h"
#include "native_support.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <sstream>

using namespace quasse_native;

namespace {

// Reusable storage, not a cache: every call replaces all original arrays and the kernel.
class Integration {
public:
    ComplexFFT fft;
    const int dimensions;
    std::vector<double> original, result, kernel;
    std::vector<double> birth, death, dFactors, lowerRoot, rootOffsetRate;
    std::vector<double> decayFullStep, decayIntegralFullStep, decayHalfStep, decayIntegralHalfStep;

    // Caller holds preparation_mutex, including destruction if an allocation fails.
    Integration(int size, int dims) : fft(size), dimensions(dims),
        original(std::size_t(2) * size * dims), result(original.size()), kernel(2 * size),
        birth(size), death(size), dFactors(size), lowerRoot(size), rootOffsetRate(size),
        decayFullStep(size), decayIntegralFullStep(size), decayHalfStep(size), decayIntegralHalfStep(size) {}

    // At fixed trait x, λ is speciation, μ is extinction, and ψ is non-removing fossil sampling.
    // E is the probability of no sampled descendants; D is the descendant-subtree likelihood.
    // Exactly solves E′ = μ − (λ + μ + ψ)E + λE² and D′ = [2λE − (λ + μ + ψ)]D backward
    // in time, without trait transport. The combined Strang calculation remains approximate.
    // This extends Fitzjohn's propagate_t (src/quasse-eqs-fftC.c) and fftR.propagate.t
    // (R/model-quasse-fftR.R), whose reaction equations have ψ = 0.
    //
    // integrate() prepares v = λ − μ − ψ, κ = sqrt(v² + 4λψ), β = (λ + μ + ψ + κ)/2.
    // For λ > 0, r₋ and r₊ are the lower and upper roots of the reaction polynomial E′ = 0.
    // β = λr₊ and s = μ/β = r₋; the latter is stored in lowerRoot.
    // rootOffsetRate stores η = β − λ = λ(r₊ − 1). The same coefficients cover λ = 0;
    // all-zero rates use s = 0. For reaction duration h, the supplied arrays contain
    // w = exp(−κh) and q = ∫₀ʰ exp(−κt) dt, named decayFactor and decayIntegral.
    //
    // For OLD E = e, A = η + λ(1 − e), B = w + Aq, and Enew = (we + sAq)/B.
    // Aq is denominatorTerm and B is denominator. Differentiating Enew in e gives J = w/B².
    // J satisfies J′ = [2λE − (λ + μ + ψ)]J, J(0) = 1, so Dnew = Dold * J is exact too.
    // Both updates use old E. Real useful bins are updated first; saved dFactors then update D.
    // E-only execution needs no D factor. Negative D is cleared without multiplying;
    // NaN follows the multiplication branch, as before. Padding/imaginary entries are unchanged.
    void propagateT(int useful, const std::vector<double>& decayFactor,
                    const std::vector<double>& decayIntegral,
                    double psi, double duration) {
        for (int i = 0; i < useful; ++i) {
            const double e = original[2 * i];
            const double denominatorTerm = (rootOffsetRate[i] + birth[i] * (1 - e)) * decayIntegral[i];
            const double denominator = decayFactor[i] + denominatorTerm;
            const double nextE = (decayFactor[i] * e + lowerRoot[i] * denominatorTerm) / denominator;
            const double dFactor = dimensions > 1 ? (decayFactor[i] / denominator) / denominator : 1;
            if (!(denominator > 0) || !std::isfinite(denominator) || !std::isfinite(nextE)
                    || !(dFactor > 0) || !std::isfinite(dFactor)) {
                std::ostringstream message;
                message << "Invalid QuaSSE reaction at bin " << i << ": lambda=" << birth[i]
                        << ", mu=" << death[i] << ", psi=" << psi << ", dt=" << duration << ", old E=" << e;
                throw std::runtime_error(message.str());
            }
            original[2 * i] = nextE;
            if (dimensions > 1) dFactors[i] = dFactor;
        }
        for (int d = 1; d < dimensions; ++d) {
            auto* values = original.data() + std::size_t(d) * 2 * fft.size;
            for (int i = 0; i < useful; ++i) {
                if (values[2 * i] < 0) values[2 * i] = 0;
                else values[2 * i] *= dFactors[i];
            }
        }
    }

    // Adjacent Strang half-Ts compose to T(dt) for fixed rates: T/2, X, T, ..., X, T/2.
    // X still restores its input boundaries. Combine only within this fixed-resolution segment;
    // zero steps remain the identity and final values reside in original after each swap and T.
    void integrate(double psi, double dt, int steps, int left, int right) {
        if (steps == 0 || dt == 0) return;
        const int useful = fft.size - left - right - 1;
        // Prepare coefficients once per segment: rates and dt stay fixed in the loop below,
        // but each reaction step still calculates its E-dependent denominator and D multiplier.
        // For positive v, η = (κ − v)/2 would cancel; use η = 2λψ/(κ + v) instead.
        // q = (1 − w)/κ uses expm1, with q = h at κ = 0 (including λ = μ, ψ = 0).
        // Full-step coefficients follow w(2h) = w(h)² and q(2h) = q(h)(1 + w(h)).
        // Refresh every segment, even when only dt or ψ changes; no exponentials or square
        // roots occur inside the timestep loop. These arrays cache storage, not model state.
        for (int i = 0; i < useful; ++i) {
            const double netRateMinusSampling = birth[i] - death[i] - psi;
            const double kappa = std::hypot(netRateMinusSampling, 2 * std::sqrt(birth[i]) * std::sqrt(psi));
            const double upperRootRate = (birth[i] + death[i] + psi + kappa) / 2;
            lowerRoot[i] = upperRootRate > 0 ? death[i] / upperRootRate : 0;
            rootOffsetRate[i] = netRateMinusSampling > 0
                ? 2 * birth[i] * psi / (kappa + netRateMinusSampling)
                : (kappa - netRateMinusSampling) / 2;
            decayHalfStep[i] = std::exp(-kappa * (dt / 2));
            decayIntegralHalfStep[i] = kappa > 0 ? -std::expm1(-kappa * (dt / 2)) / kappa : dt / 2;
            decayFullStep[i] = decayHalfStep[i] * decayHalfStep[i];
            decayIntegralFullStep[i] = decayIntegralHalfStep[i] * (1 + decayHalfStep[i]);
        }
        propagateT(useful, decayHalfStep, decayIntegralHalfStep, psi, dt / 2);
        for (int step = 0; step < steps; ++step) {
            propagate(left, right);
            original.swap(result);
            const bool last = step == steps - 1;
            propagateT(useful, last ? decayHalfStep : decayFullStep,
                       last ? decayIntegralHalfStep : decayIntegralFullStep, psi, last ? dt / 2 : dt);
        }
    }

    // Match SSEUtils: forward, complex product, normalized inverse, real clipping, real restoration.
    // Convolution reads input[i-offset], so left/right kernel extents restore opposite boundaries.
    // Imaginary entries are transformed but never clipped or restored; even negative flanks survive.
    void propagate(int left, int right) {
        const int length = 2 * fft.size;
        const int useful = fft.size - left - right - 1;
        auto* input = reinterpret_cast<double*>(fft.input.get());
        auto* output = reinterpret_cast<double*>(fft.output.get());
        for (int d = 0; d < dimensions; ++d) {
            const auto* source = original.data() + std::size_t(d) * length;
            std::copy_n(source, length, input);
            fft.transform(false);
            for (int i = 0; i < length; i += 2) {
                const double a = output[i], b = output[i + 1];
                input[i] = a * kernel[i] - b * kernel[i + 1];
                input[i + 1] = a * kernel[i + 1] + b * kernel[i];
            }
            fft.transform(true);
            // Validated flanks do not overlap: restore edges, clip only the interior, then zero padding.
            // Kernel extents restore opposite edges; imaginary entries remain untouched in every range.
            for (int i = 0; i < right; ++i) output[2 * i] = source[2 * i];
            // Keep x first to preserve NaNs and negative zero; unconditional stores allow vectorization.
            for (int i = right; i < useful - left; ++i)
                output[2 * i] = std::max(output[2 * i], 0.0);
            for (int i = useful - left; i < useful; ++i) output[2 * i] = source[2 * i];
            for (int i = useful; i < fft.size; ++i) output[2 * i] = 0;
            std::copy_n(output, length, result.data() + std::size_t(d) * length);
        }
    }

    // Validate and copy current inputs for either X-only or whole-segment execution.
    bool readInputs(JNIEnv* env, jobjectArray rows, jdoubleArray kernelArray, int left, int right) {
        const int length = 2 * fft.size;
        const std::int64_t flanks = std::int64_t(left) + right;
        if (left < 0 || right < 0 || fft.size - flanks - 1 < flanks
                || !rows || env->GetArrayLength(rows) != dimensions
                || !kernelArray || env->GetArrayLength(kernelArray) != length) {
            throw_java(env, "java/lang/IllegalArgumentException", "Invalid QuaSSE arrays or padding.");
            return false;
        }
        for (int d = 0; d < dimensions; ++d) {
            auto row = static_cast<jdoubleArray>(env->GetObjectArrayElement(rows, d));
            if (env->ExceptionCheck()) return false;
            bool valid = row && env->GetArrayLength(row) == length && !env->IsSameObject(row, kernelArray);
            for (int previous = 0; valid && previous < d; ++previous) {
                jobject other = env->GetObjectArrayElement(rows, previous);
                if (env->ExceptionCheck()) return false;
                valid = !env->IsSameObject(row, other);
                env->DeleteLocalRef(other);
            }
            if (!valid) {
                env->DeleteLocalRef(row);
                throw_java(env, "java/lang/IllegalArgumentException", "QuaSSE requires distinct rows of 2*N doubles.");
                return false;
            }
            env->GetDoubleArrayRegion(row, 0, length, original.data() + std::size_t(d) * length);
            env->DeleteLocalRef(row);
            if (env->ExceptionCheck()) return false;
        }
        env->GetDoubleArrayRegion(kernelArray, 0, length, kernel.data());
        if (env->ExceptionCheck()) return false;
        return true;
    }

    // Write only completed results into the caller's existing rows.
    void writeOutput(JNIEnv* env, jobjectArray rows, const std::vector<double>& values) {
        const int length = 2 * fft.size;
        for (int d = 0; d < dimensions; ++d) {
            auto row = static_cast<jdoubleArray>(env->GetObjectArrayElement(rows, d));
            if (env->ExceptionCheck()) return;
            env->SetDoubleArrayRegion(row, 0, length, values.data() + std::size_t(d) * length);
            env->DeleteLocalRef(row);
            if (env->ExceptionCheck()) return;
        }
    }
};

} // namespace

// Allocate independently owned resources; reject sizes whose array lengths/products overflow.
JNIEXPORT jlong JNICALL Java_SSE_QuaSSENativeIntegrator_create(JNIEnv* env, jclass, jint size, jint dims) {
    if (size <= 0 || size > std::numeric_limits<jint>::max() / 2 || dims <= 0
            || std::size_t(size) > std::numeric_limits<std::size_t>::max() / sizeof(double) / 2 / dims) {
        throw_java(env, "java/lang/IllegalArgumentException", "Invalid QuaSSE array size or dimension count.");
        return 0;
    }
    try {
        std::lock_guard lock(preparation_mutex);
        return static_cast<jlong>(reinterpret_cast<std::intptr_t>(new Integration(size, dims)));
    } catch (...) {
        translate_exception(env);
        return 0;
    }
}

// Validate/copy all inputs and finish calculation before writing any Java result row.
JNIEXPORT void JNICALL Java_SSE_QuaSSENativeIntegrator_propagateXNative(JNIEnv* env, jclass, jlong handle,
        jobjectArray rows, jdoubleArray kernel, jint left, jint right) {
    try {
        auto* integration = reinterpret_cast<Integration*>(static_cast<std::intptr_t>(handle));
        if (!integration) {
            throw_java(env, "java/lang/IllegalStateException", "QuaSSE integrator has been closed.");
            return;
        }
        if (!integration->readInputs(env, rows, kernel, left, right)) return;
        integration->propagate(left, right);
        integration->writeOutput(env, rows, integration->result);
    } catch (...) {
        translate_exception(env);
    }
}

// Copy rates once per segment; no JNI calls or allocation occur in the repeated T/X loop.
JNIEXPORT void JNICALL Java_SSE_QuaSSENativeIntegrator_integrateSegmentNative(JNIEnv* env, jclass, jlong handle,
        jobjectArray rows, jdoubleArray birth, jdoubleArray death, jdoubleArray kernel,
        jdouble psi, jdouble dt, jint steps, jint left, jint right) {
    try {
        auto* integration = reinterpret_cast<Integration*>(static_cast<std::intptr_t>(handle));
        if (!integration) {
            throw_java(env, "java/lang/IllegalStateException", "QuaSSE integrator has been closed.");
            return;
        }
        if (!integration->readInputs(env, rows, kernel, left, right)) return;
        const int useful = integration->fft.size - left - right - 1;
        if (steps < 0 || integration->dimensions < 1 || !birth || !death
                || env->GetArrayLength(birth) != useful || env->GetArrayLength(death) != useful) {
            throw_java(env, "java/lang/IllegalArgumentException", "Segment requires E with optional D, useful-bin rates and nonnegative steps.");
            return;
        }
        if (steps == 0 || dt == 0) return;
        env->GetDoubleArrayRegion(birth, 0, useful, integration->birth.data());
        if (env->ExceptionCheck()) return;
        env->GetDoubleArrayRegion(death, 0, useful, integration->death.data());
        if (env->ExceptionCheck()) return;
        integration->integrate(psi, dt, steps, left, right);
        integration->writeOutput(env, rows, integration->original);
    } catch (...) {
        translate_exception(env);
    }
}

// Java excludes execution on this owner; the shared lock excludes other FFT preparation/destruction.
JNIEXPORT void JNICALL Java_SSE_QuaSSENativeIntegrator_destroy(JNIEnv* env, jclass, jlong handle) {
    try {
        std::lock_guard lock(preparation_mutex);
        delete reinterpret_cast<Integration*>(static_cast<std::intptr_t>(handle));
    } catch (...) {
        translate_exception(env);
    }
}
