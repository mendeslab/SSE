package SSE;

/** Fixed-size, interleaved complex transforms; forward is unscaled, inverse divides by N. */
public interface ComplexFFT extends AutoCloseable {
    // Input and output must be distinct arrays of 2*N doubles; input is preserved.
    void forward(double[] input, double[] output);
    void inverse(double[] input, double[] output);
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
