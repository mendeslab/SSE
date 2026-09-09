package quasse;

import SSE.QuaSSENativeIntegrator;
import SSE.SSEUtils;
import SSE.fft.FftwFFT;
import org.junit.Test;
import java.util.Random;
import static org.junit.Assert.*;

public class QuaSSENativeIntegratorTest {
    // Model likelihoods can hide incorrect intermediate values and segment-to-segment reuse.
    // Compare full arrays and one-step calls while both Java and native integration are supported.
    @Test
    public void matchesJavaSegments() {
        Random random = new Random(128);
        for (int size : new int[]{32, 128}) {
            try (FftwFFT fft = new FftwFFT(size);
                 QuaSSENativeIntegrator nativeT = new QuaSSENativeIntegrator(size, 3)) {
                for (int steps : new int[]{0, 1, 2, 7}) {
                    int left = steps % 3, right = 3 - left, useful = size - 4;
                    double[] birth = new double[useful], death = new double[useful];
                    double[] kernel = new double[2 * size], spectrum = new double[2 * size];
                    kernel[0] = .7;
                    kernel[2] = .2;
                    kernel[2 * size - 2] = .1;
                    fft.forward(kernel, spectrum);
                    double[][] expected = new double[3][2 * size];
                    double[][] actual = new double[3][], singleSteps = new double[3][];
                    for (int i = 0; i < useful; ++i) {
                        birth[i] = .15 + .1 * random.nextDouble();
                        death[i] = .01 + .05 * random.nextDouble();
                    }
                    for (int d = 0; d < 3; ++d) {
                        for (int i = 0; i < 2 * size; ++i)
                            expected[d][i] = random.nextDouble() - (d == 0 ? 0 : .2);
                        actual[d] = expected[d].clone();
                        singleSteps[d] = expected[d].clone();
                    }
                    double[] savedBirth = birth.clone(), savedDeath = death.clone(), savedKernel = spectrum.clone();
                    double[][] scratch = new double[3][2 * size], transformed = new double[3][2 * size];
                    nativeT.integrateSegment(actual, birth, death, spectrum, .005, steps, left, right);
                    for (int step = 0; step < steps; ++step) {
                        SSEUtils.propagateEandDinTQuaSSEInPlaceSSTJavaFftService(
                                expected, scratch, birth, death, .005, useful, 2);
                        SSEUtils.propagateEandDinXQuaSSE(expected, transformed, spectrum, scratch,
                                size, left, right, 1, 2, fft);
                        nativeT.integrateSegment(singleSteps, birth, death, spectrum, .005, 1, left, right);
                    }
                    for (int d = 0; d < 3; ++d) {
                        assertArrayEquals(expected[d], actual[d], steps == 0 ? 0 : 1e-12);
                        assertArrayEquals(singleSteps[d], actual[d], 0);
                    }
                    assertArrayEquals(savedBirth, birth, 0);
                    assertArrayEquals(savedDeath, death, 0);
                    assertArrayEquals(savedKernel, spectrum, 0);
                }
            }
        }
    }

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
            double[] rates = new double[12];
            for (double[] bad : new double[][]{null, new double[11], new double[13]}) {
                assertThrows(IllegalArgumentException.class,
                        () -> nativeX.integrateSegment(valid, bad, rates, kernel, .01, 1, 1, 2));
                assertThrows(IllegalArgumentException.class,
                        () -> nativeX.integrateSegment(valid, rates, bad, kernel, .01, 1, 1, 2));
                assertEquals(7, first[0], 0);
            }
            assertThrows(IllegalArgumentException.class,
                    () -> nativeX.integrateSegment(valid, rates, rates, kernel, .01, -1, 1, 2));
            nativeX.integrateSegment(valid, rates, rates, kernel, Double.NaN, 0, 1, 2);
            assertEquals(7, first[0], 0);
            try (QuaSSENativeIntegrator onlyE = new QuaSSENativeIntegrator(16, 1)) {
                assertThrows(IllegalArgumentException.class,
                        () -> onlyE.integrateSegment(new double[][]{first}, rates, rates, kernel, .01, 0, 1, 2));
            }
        }
        nativeX.close();
        assertThrows(IllegalStateException.class, () -> nativeX.propagateX(null, null, 0, 0));
        assertThrows(IllegalStateException.class,
                () -> nativeX.integrateSegment(null, null, null, null, .01, 1, 0, 0));
    }
}
