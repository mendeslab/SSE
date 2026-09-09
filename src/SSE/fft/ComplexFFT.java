package SSE.fft;

/**
 * Lets QuaSSE switch FFT implementations without duplicating integration code.
 * N complex values occupy 2*N doubles: [real0, imaginary0, real1, imaginary1, ...].
 * Both operations preserve input and write into a distinct, caller-supplied output array.
 * Null, aliased or incorrectly sized arrays cause IllegalArgumentException;
 * transforms after close() cause IllegalStateException.
 */
public interface ComplexFFT extends AutoCloseable {
    /** Forward transform with negative exponential sign and no scaling. */
    void forward(double[] input, double[] output);
    /** Inverse transform with positive exponential sign, divided by N. */
    void inverse(double[] input, double[] output);
    /** Release owned resources. Repeated calls are harmless. */
    @Override
    void close();

    // Both implementations enforce the same layout and input-preservation contract.
    static void checkArrays(double[] input, double[] output, int size) {
        if (input == null || output == null || input.length != 2 * size || output.length != 2 * size)
            throw new IllegalArgumentException("FFT arrays must each contain 2*N doubles.");
        if (input == output)
            throw new IllegalArgumentException("FFT input and output must be distinct arrays.");
    }
}
