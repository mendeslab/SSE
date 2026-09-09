package SSE.fft;

import java.lang.ref.Cleaner;
import java.lang.ref.Reference;

/**
 * Calls jni/quasse/fft.cpp using reusable native arrays and prepared FFTW transforms.
 * Each call copies Java input into native storage and the result back into Java output.
 * No model state is cached. Explicit close releases resources; Cleaner is a safety net
 * for abandoned objects. Calls on one instance are serialized with cleanup.
 */
public final class FftwFFT implements ComplexFFT {
    private static final Cleaner CLEANER = Cleaner.create();
    // Loading is deferred until this class is initialized; failure never selects SST silently.
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

    /** Create transforms for positive size N, with 2*N fitting in an int; reject invalid sizes. */
    public FftwFFT(int size) {
        state = new NativeState();
        state.handle = create(size);
        // If cleanup registration fails, release the resources already created.
        try {
            cleanable = CLEANER.register(this, state);
        } catch (RuntimeException | Error failure) {
            state.run();
            throw failure;
        }
    }

    // Keep this owner reachable until execution ends, even though only its cleanup state is used below.
    @Override
    public void forward(double[] input, double[] output) {
        try {
            state.transform(input, output, false);
        } finally {
            Reference.reachabilityFence(this);
        }
    }

    // The native inverse applies 1/N so callers use the same convention as SST.
    @Override
    public void inverse(double[] input, double[] output) {
        try {
            state.transform(input, output, true);
        } finally {
            Reference.reachabilityFence(this);
        }
    }

    /** Release resources once, coordinating with transforms on this instance. */
    @Override
    public void close() {
        cleanable.clean();
    }

    // This object never references the Java owner. Its lock serializes use and explicit/GC cleanup.
    private static final class NativeState implements Runnable {
        // Private C++ pointer encoded as a long; zero means no live resources.
        private long handle;

        // Check and use the pointer under the same lock that protects destruction.
        synchronized void transform(double[] input, double[] output, boolean inverse) {
            if (handle == 0) throw new IllegalStateException("FFT has been closed.");
            execute(handle, input, output, inverse);
        }

        // Clear the pointer after successful destruction; subsequent cleanup is a no-op.
        @Override
        public synchronized void run() {
            if (handle == 0) return;
            destroy(handle);
            handle = 0;
        }
    }

    // Allocate arrays and prepare both directions; return their owning C++ object's pointer.
    private static native long create(int size);
    // Validate/copy arrays, execute, and copy back; inverse includes division by N.
    private static native void execute(long handle, double[] input, double[] output, boolean inverse);
    // Destroy that object; NativeState prevents concurrent use or repeated destruction.
    private static native void destroy(long handle);
}
