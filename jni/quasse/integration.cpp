// Native X propagation only. Java retains T, time stepping, grids, and normalization.
#include "SSE_QuaSSENativeIntegrator.h"
#include "native_support.h"
#include <algorithm>
#include <vector>

using namespace quasse_native;

namespace {

// Reusable storage, not a cache: every call replaces all original arrays and the kernel.
class Integration {
public:
    ComplexFFT fft;
    const int dimensions;
    std::vector<double> original, result, kernel;

    // Caller holds preparation_mutex, including destruction if an allocation fails.
    Integration(int size, int dims) : fft(size), dimensions(dims),
        original(std::size_t(2) * size * dims), result(original.size()), kernel(2 * size) {}

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
            for (int i = 0; i < fft.size; ++i) {
                if (output[2 * i] < 0 || i >= useful) output[2 * i] = 0;
                if (i < right || (i >= useful - left && i < useful)) output[2 * i] = source[2 * i];
            }
            std::copy_n(output, length, result.data() + std::size_t(d) * length);
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
        const int length = 2 * integration->fft.size;
        const std::int64_t flanks = std::int64_t(left) + right;
        if (left < 0 || right < 0 || integration->fft.size - flanks - 1 < flanks
                || !rows || env->GetArrayLength(rows) != integration->dimensions
                || !kernel || env->GetArrayLength(kernel) != length) {
            throw_java(env, "java/lang/IllegalArgumentException", "Invalid QuaSSE arrays or padding.");
            return;
        }
        for (int d = 0; d < integration->dimensions; ++d) {
            auto row = static_cast<jdoubleArray>(env->GetObjectArrayElement(rows, d));
            if (env->ExceptionCheck()) return;
            bool valid = row && env->GetArrayLength(row) == length && !env->IsSameObject(row, kernel);
            for (int previous = 0; valid && previous < d; ++previous) {
                jobject other = env->GetObjectArrayElement(rows, previous);
                if (env->ExceptionCheck()) return;
                valid = !env->IsSameObject(row, other);
                env->DeleteLocalRef(other);
            }
            if (!valid) {
                env->DeleteLocalRef(row);
                throw_java(env, "java/lang/IllegalArgumentException", "QuaSSE requires distinct rows of 2*N doubles.");
                return;
            }
            env->GetDoubleArrayRegion(row, 0, length, integration->original.data() + std::size_t(d) * length);
            env->DeleteLocalRef(row);
            if (env->ExceptionCheck()) return;
        }
        env->GetDoubleArrayRegion(kernel, 0, length, integration->kernel.data());
        if (env->ExceptionCheck()) return;
        integration->propagate(left, right);
        for (int d = 0; d < integration->dimensions; ++d) {
            auto row = static_cast<jdoubleArray>(env->GetObjectArrayElement(rows, d));
            if (env->ExceptionCheck()) return;
            env->SetDoubleArrayRegion(row, 0, length, integration->result.data() + std::size_t(d) * length);
            env->DeleteLocalRef(row);
            if (env->ExceptionCheck()) return;
        }
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
