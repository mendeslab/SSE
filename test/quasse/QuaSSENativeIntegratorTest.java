package quasse;

import SSE.QuaSSENativeIntegrator;
import SSE.SSEUtils;
import SSE.fft.FftwFFT;
import org.junit.Test;
import java.util.Random;
import static org.junit.Assert.*;

public class QuaSSENativeIntegratorTest {
    // Full likelihoods do not expose imaginary entries, clipping, or repeated-input ownership.
    // These comparisons remain needed while Java and native X are independently implemented.
    @Test
    public void matchesJavaX() {
        Random random = new Random(127);
        for (int size : new int[]{32, 128}) {
            try (FftwFFT fft = new FftwFFT(size);
                 QuaSSENativeIntegrator nativeX = new QuaSSENativeIntegrator(size, 3)) {
                for (int repeat = 0; repeat < 6; ++repeat) {
                    double[][] javaValues = new double[3][2 * size];
                    double[][] nativeValues = new double[3][];
                    for (int d = 0; d < 3; ++d) {
                        for (int i = 0; i < 2 * size; ++i) javaValues[d][i] = random.nextDouble() - 0.4;
                        nativeValues[d] = javaValues[d].clone();
                    }
                    double[] kernel = new double[2 * size];
                    for (int i = 0; i < kernel.length; ++i) kernel[i] = random.nextDouble();
                    double[] savedKernel = kernel.clone();
                    int left = repeat % 3, right = 3 - left;
                    SSEUtils.propagateEandDinXQuaSSE(javaValues, new double[3][2 * size], kernel,
                            new double[3][2 * size], size, left, right, 1, 2, fft);
                    nativeX.propagateX(nativeValues, kernel, left, right);
                    for (int d = 0; d < 3; ++d) assertArrayEquals(javaValues[d], nativeValues[d], 1e-12);
                    assertArrayEquals(savedKernel, kernel, 0);
                }
            }
        }
    }

    // JNI must reject malformed storage before writing results and serialize close with use.
    // Model tests cannot exercise this public boundary; obsolete if native ownership is removed.
    @Test
    public void validatesBeforeWritingAndCloses() {
        assertThrows(IllegalArgumentException.class, () -> new QuaSSENativeIntegrator(0, 2));
        assertThrows(IllegalArgumentException.class, () -> new QuaSSENativeIntegrator(16, 0));
        double[] first = new double[32];
        first[0] = 7;
        double[] kernel = new double[32];
        QuaSSENativeIntegrator nativeX = new QuaSSENativeIntegrator(16, 2);
        try (nativeX) {
            for (double[][] bad : new double[][][]{null, {first}, {first, null},
                    {first, new double[31]}, {first, first}, {first, kernel}}) {
                assertThrows(IllegalArgumentException.class, () -> nativeX.propagateX(bad, kernel, 1, 2));
                assertEquals(7, first[0], 0);
            }
            double[][] valid = {first, new double[32]};
            assertThrows(IllegalArgumentException.class, () -> nativeX.propagateX(valid, null, 1, 2));
            assertThrows(IllegalArgumentException.class, () -> nativeX.propagateX(valid, kernel, -1, 2));
            assertThrows(IllegalArgumentException.class,
                    () -> nativeX.propagateX(valid, kernel, Integer.MAX_VALUE, 2));
            nativeX.propagateX(valid, kernel, 1, 2);
        }
        nativeX.close();
        assertThrows(IllegalStateException.class, () -> nativeX.propagateX(null, null, 0, 0));
    }
}
