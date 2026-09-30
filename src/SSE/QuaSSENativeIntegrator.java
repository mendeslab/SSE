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

    /**
     * Repeat T(dt/2), X(dt), T(dt/2); row zero is E and remaining rows are D.
     * E is the probability of no sampled descendants; D is a descendant-subtree likelihood.
     * Allows E alone or E plus D rows; storage follows propagateX's contract. Rates and the
     * kernel are copied once per call and reused, while the E-dependent D factors change each step.
     * Inputs other than E/D are unchanged. Zero steps or zero dt are the identity.
     *
     * @param esDs distinct rows of 2*N alternating real/imaginary entries, updated in place
     * @param birthRates speciation rates λ, with N-leftPadding-rightPadding-1 entries
     * @param deathRates extinction rates μ, with the same useful-bin count
     * @param kernelSpectrum Fourier transform of the real transport kernel for the full dt
     * @param psi nonnegative fossil sampling rate, constant across bins
     * @param dt full Strang-step duration, increasing backward from the present
     * @param steps nonnegative number of steps; this method performs no likelihood normalization
     * @param leftPadding kernel support in bins to the left; flanks must not overlap
     * @param rightPadding kernel support in bins to the right; flanks must not overlap
     */
    public void integrateSegment(double[][] esDs, double[] birthRates, double[] deathRates,
            double[] kernelSpectrum, double psi, double dt, int steps, int leftPadding, int rightPadding) {
        try {
            state.integrateSegment(esDs, birthRates, deathRates, kernelSpectrum, psi, dt, steps, leftPadding, rightPadding);
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
                double[] kernel, double psi, double dt, int steps, int left, int right) {
            if (handle == 0) throw new IllegalStateException("QuaSSE integrator has been closed.");
            integrateSegmentNative(handle, esDs, birth, death, kernel, psi, dt, steps, left, right);
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
            double[] death, double[] kernel, double psi, double dt, int steps, int leftPadding, int rightPadding);
}
