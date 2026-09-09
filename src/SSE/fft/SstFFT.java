package SSE.fft;

import org.shared.fft.JavaFftService;

/**
 * Uses the Shared Scientific Toolbox (SST), bundled in lib/sst.jar.
 * Default Java implementation, also retained for comparison with FFTW.
 * Delegates to JavaFftService; keeping this object does not make SST reuse its
 * internal scratch arrays or transform preparation between calls.
 */
public final class SstFFT implements ComplexFFT {
    private final int[] dimensions;
    private final int size;
    private final JavaFftService service;
    private boolean closed;

    /** Create a one-dimensional transform of positive size N, with 2*N fitting in an int. */
    public SstFFT(int size) {
        this(size, new JavaFftService());
    }

    // Retain service-based SSEUtils callers while the active solver uses ComplexFFT directly.
    public SstFFT(int size, JavaFftService service) {
        this(new int[] {size}, service);
    }

    // Preserve the older utility's multidimensional SST argument without changing active QuaSSE layout.
    public SstFFT(int[] dimensions, JavaFftService service) {
        int count = 1;
        for (int length : dimensions) {
            if (length <= 0 || count > Integer.MAX_VALUE / 2 / length)
                throw new IllegalArgumentException("FFT size must be positive and 2*N must fit in an int.");
            count *= length;
        }
        this.dimensions = dimensions.clone();
        size = count;
        this.service = service;
    }

    // Preserve SST's forward convention and caller-owned arrays.
    @Override
    public void forward(double[] input, double[] output) {
        check(input, output);
        service.fft(dimensions, input, output);
    }

    // SST already divides its inverse result by N.
    @Override
    public void inverse(double[] input, double[] output) {
        check(input, output);
        service.ifft(dimensions, input, output);
    }

    // Check closure as well as the shared array contract before entering SST.
    private void check(double[] input, double[] output) {
        if (closed) throw new IllegalStateException("FFT has been closed.");
        ComplexFFT.checkArrays(input, output, size);
    }

    /** Mark closed for the shared contract; this implementation owns no native resources. */
    @Override
    public void close() {
        closed = true;
    }
}
