package test;

import SSE.fft.ComplexFFT;
import SSE.fft.FftwFFT;
import SSE.fft.SstFFT;
import org.junit.Assert;
import org.junit.Test;

import java.util.Random;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Checks transform conventions and resource lifetime independently of likelihood calculations.
 * ant test uses SST; ant test-native sets test.quasse.fft=fftw and requires the native library.
 */
public class ComplexFFTTest {
    // Select native execution only in the explicit target; absence must fail, never skip or fall back.
    private ComplexFFT create(int size) {
        return "fftw".equals(System.getProperty("test.quasse.fft")) ? new FftwFFT(size) : new SstFFT(size);
    }

    // Individual transforms catch sign/scale errors that round trips and likelihoods can hide.
    // This group is needed while implementations differ; remove if only one transform remains.
    @Test
    public void conventionsAndReuse() {
        Random random = new Random(7291);
        for (int n : new int[] {1, 7, 32, 128, 1024, 4096}) {
            try (ComplexFFT fft = create(n); ComplexFFT sst = new SstFFT(n)) {
                double[] input = new double[2*n], actual = new double[2*n], expected = new double[2*n];
                for (int sample = 0; sample < 4; sample++) {
                    Arrays.fill(input, 0);
                    if (sample == 1) {
                        for (int i = 0; i < input.length; i += 2) input[i] = 1;
                    } else if (sample == 2) {
                        input[2 * (1 % n)] = 1;
                    } else if (sample == 3) {
                        for (int i = 0; i < input.length; i++) input[i] = random.nextDouble() - .5;
                    }
                    double[] saved = input.clone();
                    fft.forward(input, actual);
                    sst.forward(input, expected);
                    Assert.assertArrayEquals(expected, actual, 2e-12*n);
                    Assert.assertArrayEquals(saved, input, 0);
                    if (sample == 2) {
                        // A unit impulse at index 1 has spectrum exp(-2πik/N), or 1 when N=1.
                        for (int k = 0; k < n; k++) {
                            Assert.assertEquals(Math.cos(2*Math.PI*k/n), actual[2*k], 1e-12);
                            Assert.assertEquals(-Math.sin(2*Math.PI*k/n), actual[2*k+1], 1e-12);
                        }
                    }
                    double[] roundTrip = new double[2*n];
                    fft.inverse(actual, roundTrip);
                    Assert.assertArrayEquals(input, roundTrip, 1e-12);
                    fft.inverse(input, actual);
                    sst.inverse(input, expected);
                    Assert.assertArrayEquals(expected, actual, 1e-12);
                    Assert.assertArrayEquals(saved, input, 0);
                }
            }
        }
    }

    // Likelihood tests never deliberately misuse FFT arrays or resources. Keep these checks while
    // Java owns native memory; they protect against invalid access and double destruction.
    @Test
    public void invalidArgumentsAndClosure() {
        Assert.assertThrows(IllegalArgumentException.class, () -> create(0));
        Assert.assertThrows(IllegalArgumentException.class, () -> create(-1));
        Assert.assertThrows(IllegalArgumentException.class, () -> create(Integer.MAX_VALUE));
        ComplexFFT fft = create(8);
        try {
            double[] values = new double[16];
            Assert.assertThrows(IllegalArgumentException.class, () -> fft.forward(values, values));
            Assert.assertThrows(IllegalArgumentException.class, () -> fft.forward(null, values));
            Assert.assertThrows(IllegalArgumentException.class, () -> fft.inverse(values, new double[14]));
            Assert.assertThrows(IllegalArgumentException.class, () -> fft.forward(new double[18], values));
        } finally {
            fft.close();
            fft.close();
        }
        Assert.assertThrows(IllegalStateException.class, () -> fft.forward(new double[16], new double[16]));
    }

    // Independent BEAST analyses may prepare/release FFTW concurrently, unlike the serial reference
    // tests. Keep this small stress check while the JNI library manages shared FFTW preparation.
    @Test
    public void concurrentCreationAndRelease() throws Exception {
        ExecutorService threads = Executors.newFixedThreadPool(2);
        try {
            // Each worker repeatedly prepares, executes and closes independent transforms.
            Runnable work = () -> {
                for (int i = 0; i < 20; i++) {
                    try (ComplexFFT fft = create(32)) {
                        double[] input = new double[64], output = new double[64];
                        input[0] = 1;
                        fft.forward(input, output);
                        for (int k = 0; k < 32; k++) Assert.assertEquals(1, output[2*k], 1e-12);
                    }
                }
            };
            Future<?> first = threads.submit(work), second = threads.submit(work);
            first.get();
            second.get();
        } finally {
            threads.shutdownNow();
        }
    }
}
