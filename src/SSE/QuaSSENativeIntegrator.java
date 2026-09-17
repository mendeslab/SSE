package SSE;

import java.lang.ref.Cleaner;
import java.lang.ref.Reference;

/**
 * Native QuaSSE calculations, with independently owned arrays and FFTW transforms.
 * Segment calls copy current E/D, rates and kernel once, then repeat Strang steps in native storage.
 * No model state is cached between calls. Close explicitly; Cleaner releases abandoned instances.
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

    /** Repeat T(dt/2), X(dt), T(dt/2); row zero is E and remaining rows are D.
     * Rates have N-left-right-1 entries; the supplied kernel must represent the full dt.
     * Requires at least one D row and nonnegative steps; E/D storage follows propagateX's contract.
     * Inputs other than E/D are unchanged. Zero steps are the identity, with no normalization.
     */
    public void integrateSegment(double[][] esDs, double[] birthRates, double[] deathRates,
            double[] kernelSpectrum, double dt, int steps, int leftPadding, int rightPadding) {
        try {
            state.integrateSegment(esDs, birthRates, deathRates, kernelSpectrum, dt, steps, leftPadding, rightPadding);
        } finally {
            Reference.reachabilityFence(this);
        }
    }

    // No reference to the owner; the same lock protects execution and explicit/GC cleanup.
    private static final class NativeState implements Runnable {
        private long handle;

        // Check and use the native pointer while destruction is excluded.
        synchronized void propagateX(double[][] esDs, double[] kernel, int left, int right) {
            if (handle == 0) throw new IllegalStateException("QuaSSE integrator has been closed.");
            propagateXNative(handle, esDs, kernel, left, right);
        }

        // Hold the ownership lock for the whole segment, excluding explicit and automatic cleanup.
        synchronized void integrateSegment(double[][] esDs, double[] birth, double[] death,
                double[] kernel, double dt, int steps, int left, int right) {
            if (handle == 0) throw new IllegalStateException("QuaSSE integrator has been closed.");
            integrateSegmentNative(handle, esDs, birth, death, kernel, dt, steps, left, right);
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
    private static native void integrateSegmentNative(long handle, double[][] esDs, double[] birth,
            double[] death, double[] kernel, double dt, int steps, int leftPadding, int rightPadding);
}
