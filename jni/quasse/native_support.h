#pragma once

// Internal FFT resources and JNI error handling shared by transform and integration entry points.
#include <jni.h>
#include <fftw3.h>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>

namespace quasse_native {

// FFTW only guarantees concurrent execution, not preparation/destruction. This lock covers
// all preparation/destruction in this module, including Java Cleaner calls. No global cleanup.
// This does not coordinate unrelated FFTW libraries; fftw_cleanup() could invalidate their plans.
inline std::mutex preparation_mutex;

// Owns aligned arrays and FFTW's prepared transform instructions (called "plans").
// Java's per-object lock excludes concurrent execution/cleanup; the lock above covers preparation.
class ComplexFFT {
public:
    const int size;
    std::unique_ptr<fftw_complex, decltype(&fftw_free)> input{nullptr, fftw_free};
    std::unique_ptr<fftw_complex, decltype(&fftw_free)> output{nullptr, fftw_free};
    fftw_plan forward = nullptr;
    fftw_plan inverse = nullptr;

    // Caller holds preparation_mutex. Arrays and partially prepared transforms are released on failure.
    explicit ComplexFFT(int n) : size(n) {
        input.reset(fftw_alloc_complex(size));
        output.reset(fftw_alloc_complex(size));
        if (!input || !output) throw std::bad_alloc();
        // Measure candidate algorithms once for repeated execution. Preparation overwrites arrays;
        // callers supply actual input only after construction.
        forward = fftw_plan_dft_1d(size, input.get(), output.get(), FFTW_FORWARD, FFTW_MEASURE);
        inverse = fftw_plan_dft_1d(size, input.get(), output.get(), FFTW_BACKWARD, FFTW_MEASURE);
        if (!forward || !inverse) {
            // A failed constructor does not run this destructor; release any completed plan here.
            if (forward) fftw_destroy_plan(forward);
            if (inverse) fftw_destroy_plan(inverse);
            throw std::runtime_error("FFTW could not prepare the requested transforms.");
        }
    }

    // Caller holds preparation_mutex; member-owned arrays are freed after the transforms.
    ~ComplexFFT() {
        fftw_destroy_plan(forward);
        fftw_destroy_plan(inverse);
    }

    ComplexFFT(const ComplexFFT&) = delete;
    ComplexFFT& operator=(const ComplexFFT&) = delete;

    // FFTW forward is unscaled; normalize its inverse by 1/N to match SST's convention.
    // Multiplying by the reciprocal avoids division for every component.
    void transform(bool backwards) {
        fftw_execute(backwards ? inverse : forward);
        if (backwards) {
            auto* values = reinterpret_cast<double*>(output.get());
            const double inverseSize = 1.0 / size;
            for (int i = 0; i < 2 * size; ++i) values[i] *= inverseSize;
        }
    }
};

// Preserve pending JNI exceptions, including allocation failures while finding the exception class.
inline void throw_java(JNIEnv* env, const char* type, const char* message) {
    if (env->ExceptionCheck()) return;
    jclass exception = env->FindClass(type);
    if (exception) env->ThrowNew(exception, message);
}

// Called only inside catch: translate C++ failures to Java exceptions, never through JVM frames.
inline void translate_exception(JNIEnv* env) {
    try {
        throw;
    } catch (const std::bad_alloc&) {
        throw_java(env, "java/lang/OutOfMemoryError", "Cannot allocate QuaSSE FFTW resources.");
    } catch (const std::exception& error) {
        throw_java(env, "java/lang/IllegalStateException", error.what());
    } catch (...) {
        throw_java(env, "java/lang/IllegalStateException", "Unexpected QuaSSE FFTW failure.");
    }
}

} // namespace
