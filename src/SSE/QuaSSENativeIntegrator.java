package SSE;

import java.lang.ref.Cleaner;
import java.lang.ref.Reference;

/**
 * Native QuaSSE X propagation, with independently owned arrays and FFTW transforms.
 * Each call copies current E/D and kernel inputs; no model state is cached between calls.
 * Per-step copying is temporary; whole-segment integration can later transfer inputs once.
 * T propagation remains in Java. Close explicitly; Cleaner releases abandoned instances.
 */
public final class QuaSSENativeIntegrator implements AutoCloseable {
    private static final Cleaner CLEANER = Cleaner.create();
    static {
        try {
            System.loadLibrary("sse_quasse");
        } catch (UnsatisfiedLinkError cause) {
            UnsatisfiedLinkError error = new UnsatisfiedLinkError(
                    "Cannot load sse_quasse. Run ant native and put its output directory on "
                    + "LD_LIBRARY_PATH (or java.library.path in the IDE). " + cause.getMessage());
            error.initCause(cause);
            throw error;
        }
    }

    private final NativeState state;
    private final Cleaner.Cleanable cleanable;

    /** Allocate storage for dimensions distinct arrays of 2*size alternating real/imaginary values. */
    public QuaSSENativeIntegrator(int size, int dimensions) {
        state = new NativeState();
        state.handle = create(size, dimensions);
        try {
            cleanable = CLEANER.register(this, state);
        } catch (RuntimeException | Error failure) {
            state.run();
            throw failure;
        }
    }

    /** Convolve, clip real values, and restore real boundaries; leave kernelSpectrum unchanged.
     * Rows must be distinct from each other and the kernel. Padding must leave nonoverlapping flanks.
     */
    public void propagateX(double[][] esDs, double[] kernelSpectrum, int leftPadding, int rightPadding) {
        try {
            state.propagateX(esDs, kernelSpectrum, leftPadding, rightPadding);
        } finally {
            Reference.reachabilityFence(this);
        }
    }

    @Override
    public void close() {
        cleanable.clean();
    }

    // No reference to the owner; the same lock protects execution and explicit/GC cleanup.
    private static final class NativeState implements Runnable {
        private long handle;

        // Check and use the native pointer while destruction is excluded.
        synchronized void propagateX(double[][] esDs, double[] kernel, int left, int right) {
            if (handle == 0) throw new IllegalStateException("QuaSSE integrator has been closed.");
            propagateXNative(handle, esDs, kernel, left, right);
        }

        // Repeated cleanup is harmless, including close followed by garbage collection.
        @Override
        public synchronized void run() {
            if (handle == 0) return;
            destroy(handle);
            handle = 0;
        }
    }

    private static native long create(int size, int dimensions);
    private static native void propagateXNative(long handle, double[][] esDs, double[] kernel,
                                                int leftPadding, int rightPadding);
    private static native void destroy(long handle);
}
