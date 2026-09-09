package SSE;

import java.lang.ref.Cleaner;
import java.lang.ref.Reference;

/** Reuses FFTW arrays and prepared transforms; no model or likelihood state lives in C++. */
public final class FftwFFT implements ComplexFFT {
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

    // Register cleanup only after successful creation, releasing native resources if registration fails.
    public FftwFFT(int size) {
        state = new NativeState();
        state.handle = create(size);
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

    @Override
    public void close() {
        cleanable.clean();
    }

    // This object never references the Java owner. Its lock serializes use and explicit/GC cleanup.
    private static final class NativeState implements Runnable {
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

    private static native long create(int size);
    private static native void execute(long handle, double[] input, double[] output, boolean inverse);
    private static native void destroy(long handle);
}
