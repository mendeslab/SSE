package SSE;

import SSE.fft.ComplexFFT;
import SSE.fft.SstFFT;

/*
 * diversitree's (src/quasse-eqs-fftC.c) variables explained:
 *
 * nx = total number of bins for quantitative ch (default = 1024)
 *
 * nkl and nkr:
 * in the normal kernel (giving the pdf of trait value change),
 * specify how many bins on the right (and left, respectively -- and yes, it's reversed!)
 * are non-zero,
 *
 * and
 *
 * in the D and E arrays, they determine how many bins on the left and on the right
 * are not updated by propagate.x()
 *
 * npad = nkl + nkr
 *
 * ndat (my 'nUsefulTraitBins') = nx - npad
 *
 * nd = number of dimensions for each plan (i.e., number of equations per quantitative ch)
 * x = flat, contains all E's followed by all D's (size 2 * nx; my 2D-array 'esDs')
 * wrk = scratch space, used to store and restore x values in left and right padding
 *
 * vars = deep copy of E's, used locally in functions, used locally in propagate_t
 * d = deep copy of x skipping E's (contains D's), used locally in propagate_t
 * dd = deep copy of wrk, only used for D's, used locally in propagate_t to save values and restore them later, also for storing terms used in the math
 */

import org.jtransforms.fft.DoubleFFT_1D;
import org.shared.array.ComplexArray;
import org.shared.array.RealArray;
import org.shared.fft.JavaFftService;

import java.util.Arrays;

public class SSEUtils {

    private static final double SQRT2PI = Math.sqrt(2 * Math.PI);

    /**
     * Propagate E and D through a reaction step in the contiguous real layout used by JTransforms.
     * Row zero holds E; each subsequent row holds a D component. The shared reaction calculation
     * below explains the equations and their numerical evaluation. A zero duration is a no-op.
     *
     * @param esDsAtNode E/D rows updated in place, with useful real bins at consecutive indices
     * @param scratchAtNode workspace; row 1 holds one D multiplier per useful bin when D is present
     * @param birthRate speciation rates λ, one per useful trait bin
     * @param deathRate extinction rates μ, one per useful trait bin
     * @param psi nonnegative fossil sampling rate, constant across bins
     * @param dt reaction duration, increasing backward from the present
     * @param nUsefulTraitBins number of real bins to update, excluding the zero-padded tail
     * @param nDimensionsD number of D rows after E; zero permits E-only rows and workspace
     */
    public static void propagateEandDinTQuaSSEInPlace(double[][] esDsAtNode, double[][] scratchAtNode,
            double[] birthRate, double[] deathRate, double psi, double dt, int nUsefulTraitBins, int nDimensionsD) {
        propagateReaction(esDsAtNode, scratchAtNode, birthRate, deathRate, psi, dt, nUsefulTraitBins, nDimensionsD, 1);
    }

    /**
     * Propagate E and D through a reaction step in the interleaved complex FFT layout.
     * Row zero holds E; each subsequent row holds a D component. Real values occupy even indices;
     * imaginary values and the padded tail are unchanged. A zero duration is a no-op.
     *
     * @param esDsAtNode E/D rows updated in place, alternating real and imaginary entries
     * @param scratchAtNode workspace; D factors in row 1 are contiguous, not interleaved
     * @param birthRate speciation rates λ, one per useful trait bin
     * @param deathRate extinction rates μ, one per useful trait bin
     * @param psi nonnegative fossil sampling rate, constant across bins
     * @param dt reaction duration, increasing backward from the present
     * @param nUsefulTraitBins number of real bins to update, excluding the zero-padded tail
     * @param nDimensionsD number of D rows after E; zero permits E-only rows and workspace
     */
    public static void propagateEandDinTQuaSSEInPlaceSSTJavaFftService(double[][] esDsAtNode, double[][] scratchAtNode,
            double[] birthRate, double[] deathRate, double psi, double dt, int nUsefulTraitBins, int nDimensionsD) {
        propagateReaction(esDsAtNode, scratchAtNode, birthRate, deathRate, psi, dt, nUsefulTraitBins, nDimensionsD, 2);
    }

    // At fixed trait x, λ is speciation, μ is extinction, and ψ is non-removing fossil sampling.
    // E is the probability of no sampled descendants; D is the descendant-subtree likelihood.
    // With constant nonnegative rates, this exactly solves the reaction ODEs backward in time:
    //     E′ = μ − (λ + μ + ψ)E + λE²
    //     D′ = [2λE − (λ + μ + ψ)]D.
    // Trait drift/diffusion is handled separately; the combined Strang calculation is not exact.
    // This extends Fitzjohn's fftR.propagate.t (R/model-quasse-fftR.R) and propagate_t
    // (src/quasse-eqs-fftC.c), whose reaction equations have ψ = 0.
    //
    // Rate coefficients: v = λ − μ − ψ, κ = sqrt(v² + 4λψ), β = (λ + μ + ψ + κ)/2.
    // For λ > 0, r₋ and r₊ are the lower and upper roots of the reaction polynomial E′ = 0.
    // β = λr₊ and s = μ/β = r₋. The code calls these
    // upperRootRate and lowerRoot. The rootOffsetRate η = β − λ = λ(r₊ − 1).
    // These coefficients also give the linear flow when λ = 0; all-zero rates use s = 0.
    //
    // For duration h = dt, define w = exp(−κh) and q = ∫₀ʰ exp(−κt) dt = (1 − w)/κ.
    // The code calls them decayFactor and decayIntegral. With the OLD E value E₀ = ithE:
    //     A = η + λ(1 − E₀), B = w + Aq
    //     Enew = (wE₀ + sAq)/B, Dnew = Dold * w/B².
    // Here Aq is denominatorTerm and B is denominator. Differentiating Enew in E₀ gives
    // J = w/B²; J′ = [2λE − (λ + μ + ψ)]J, J(0) = 1, so J is exactly the D multiplier.
    // Both updates must therefore use the same old E, even though E is overwritten first.
    //
    // Evaluate q with expm1; its κ = 0 limit is h (including λ = μ with ψ = 0).
    // For positive v, rationalize η as 2λψ/(κ + v) to avoid subtracting nearly equal rates.
    // See LSU/TODO/AnalyticSolution.tex and AnalyticSolution2.tex for the Riccati derivation.
    //
    // stride is 1 for contiguous real data and 2 for interleaved complex data. Only real useful
    // bins are updated. D multipliers are contiguous in scratchAtNode[1], regardless of stride;
    // E-only calls neither calculate those multipliers nor access a D workspace row.
    private static void propagateReaction(double[][] esDsAtNode, double[][] scratchAtNode, double[] birthRate,
            double[] deathRate, double psi, double dt, int nUsefulTraitBins, int nDimensionsD, int stride) {
        if (dt == 0) return;
        for (int i = 0; i < nUsefulTraitBins; ++i) {
            double ithLambda = birthRate[i];
            double ithMu = deathRate[i];
            double ithE = esDsAtNode[0][stride * i];
            double netRateMinusSampling = ithLambda - ithMu - psi;
            double kappa = Math.hypot(netRateMinusSampling, 2 * Math.sqrt(ithLambda) * Math.sqrt(psi));
            double upperRootRate = (ithLambda + ithMu + psi + kappa) / 2;
            double lowerRoot = upperRootRate > 0 ? ithMu / upperRootRate : 0;
            double rootOffsetRate = netRateMinusSampling > 0
                    ? 2 * ithLambda * psi / (kappa + netRateMinusSampling)
                    : (kappa - netRateMinusSampling) / 2;
            double decayFactor = Math.exp(-kappa * dt);
            double decayIntegral = kappa > 0 ? -Math.expm1(-kappa * dt) / kappa : dt;
            double denominatorTerm = (rootOffsetRate + ithLambda * (1 - ithE)) * decayIntegral;
            double denominator = decayFactor + denominatorTerm;
            double nextE = (decayFactor * ithE + lowerRoot * denominatorTerm) / denominator;
            double dFactor = nDimensionsD > 0 ? (decayFactor / denominator) / denominator : 1;
            if (!(denominator > 0) || !Double.isFinite(denominator) || !Double.isFinite(nextE)
                    || !(dFactor > 0) || !Double.isFinite(dFactor))
                throw new ArithmeticException("Invalid QuaSSE reaction at bin " + i + ": lambda=" + ithLambda
                        + ", mu=" + ithMu + ", psi=" + psi + ", dt=" + dt + ", old E=" + ithE);
            esDsAtNode[0][stride * i] = nextE;
            if (nDimensionsD > 0) scratchAtNode[1][i] = dFactor;
        }
        // Preserve existing negative-D roundoff treatment; padding and imaginary entries are untouched.
        for (int d = 1; d <= nDimensionsD; ++d)
            for (int i = 0; i < nUsefulTraitBins; ++i) {
                int j = stride * i;
                if (esDsAtNode[d][j] < 0) esDsAtNode[d][j] = 0;
                else esDsAtNode[d][j] *= scratchAtNode[1][i];
            }
    }

    /*
     * Propagate E's and D's along the x-axis (quantitative trait).
     *
     * The normal kernel (fftfY) comes in FFTed already.
     * The math work (FFTs/iFFts and convolution) is done in convolveInPlace.
     * This function here basically reorganizes the arrays, keeping track
     * of the flanking bins (not letting those be affected by propagate in x),
     * squashing negative numbers back to 0.0, and setting the last
     * (nLeftFlankBins + nRightFlankBins) to 0.0.
     *
     * Uses JTransforms for FFT/inverse-FFT
     *
     * @param   esDs    2D-array containing E's followed by D's (e.g., for QuaSSE, esDs[0]=Es, esDs[1]=Ds)
     * @param   scratch 2D-array for storing/restoring flanking values and other math terms
     * @param   fftfY  Normal kernel (already FFT-ed) giving the density for changes in value of the quantitative trait
     * @param   nXbins  total number of bins resulting from discretizing quantitative trait-change normal kernel (fY and each row of esDs will have these many nXbins)
     * @param   nLeftFlankBins  how many bins on the right side of kernel are non-zero
     * @param   nRightFlankBins  how many bins on the left side of kernel are non-zero
     * @param   nDimensionsE number of E equations (dimensions in plan) to solve for each quantitative ch
     * @param   nDimensionsD number of D equations (dimensions in plan) to solve for each quantitative ch
     * @param   fft instance of DoubleFFT_1D that will carry out FFT and inverse-FFT
     */
    public static void propagateEandDinXQuaLike(double[][] esDsAtNode, double[][] scratchAtNode, double[] fftfY, int nXbins, int nLeftFlankBins, int nRightFlankBins, int nDimensionsE, int nDimensionsD, DoubleFFT_1D fft) {

        // debugging
        // System.out.println("\nEntering propagateEandDinXQuaLike");
        // System.out.println("esDsAtNode = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode = " + Arrays.toString(scratchAtNode[1]));

        // Convolution reads input[i-offset]: preserve nRightFlankBins on the left and
        // nLeftFlankBins on the right. Keep this backend's existing E/D restoration policy.
        int nPad = nLeftFlankBins + nRightFlankBins + 1;
        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int j = 0; j < nRightFlankBins; ++j) {
                scratchAtNode[ithDim][j] = esDsAtNode[ithDim][j];
            }
        }

        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int j = (nXbins - nPad - nLeftFlankBins); j < (nXbins - nPad); ++j) {
                scratchAtNode[ithDim][j] = esDsAtNode[ithDim][j];
            }
        }

        // debugging
        // System.out.println("After copying left and right flanks");
        // System.out.println("esDsAtNode = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode = " + Arrays.toString(scratchAtNode[1]));

        // fY is already FFTed (in QuaSSEDistribution -> populatefY()), then FFT scratch, inverse-FFT scratch, result is left in scratch
        convolveInPlace(esDsAtNode, fftfY, nDimensionsE, nDimensionsD, fft);

        // debugging
        // System.out.println("After convolve");
        // System.out.println("esDsAtNode = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode = " + Arrays.toString(scratchAtNode[1]));

        // System.out.println("esDsAtNode[1] after convolve: " + Arrays.toString(esDsAtNode[1]));

        int nItems2Copy = nXbins - nLeftFlankBins - nRightFlankBins;
        // System.out.println("nItems2Copy=" + nItems2Copy);
        double scaleBy = 1.0 / nXbins;

        /*
         * Note that while we are grabbing every other element (so we can get just the real parts)
         * as we iterate over esDsAtNode, everyOtherToHeadInPlace moves these elements to the head
         * of esDsAtNode without any skipping, i.e.,
         * [ ... skipFirst_nLeftFlankBins ..., real1, real2, real3, ...]
         */
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            everyOtherToHeadInPlace(esDsAtNode[ithDim], nXbins, nRightFlankBins, nLeftFlankBins, 2, scaleBy); // grabbing real part and scaling by 1/nXbins

            for (int i=0; i<nXbins; ++i) {
                // if negative value, set to 0.0
                // if at last (nLeftFlankBins + nRightFlankBins) items, set to 0.0
                if (esDsAtNode[ithDim][i] < 0.0 || i >= (nItems2Copy-1)) esDsAtNode[ithDim][i] = 0.0;
            }
        }

        // Restore the boundary strips whose convolutions would read unavailable input.
        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int j = 0; j < nRightFlankBins; ++j) {
                esDsAtNode[ithDim][j] = scratchAtNode[ithDim][j];
            }
        }

        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int j = (nXbins - nPad - nLeftFlankBins); j < (nXbins - nPad); ++j) {
                esDsAtNode[ithDim][j] = scratchAtNode[ithDim][j];
            }
        }

        // debugging
        // System.out.println("After scaling and zero-ing");
        // System.out.println("esDsAtNode (D's) = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode = " + Arrays.toString(scratchAtNode[1]));
    }

    /*
     * Version for SST using ComplexArray and RealArray (see unit tests in PropagatesQuaSSETest)
     */
    public static void propagateEandDinXQuaLikeSSTModalFftService(double[][] esDsAtNode, ComplexArray fftFYCA, double[][] scratchAtNode, RealArray scratchRA, int nXbins, int nLeftFlankBins, int nRightFlankBins, int nDimensionsE, int nDimensionsD, DoubleFFT_1D fft) {
        // Convolution reads input[i-offset]: preserve nRightFlankBins on the left and
        // nLeftFlankBins on the right. Keep this backend's existing E/D restoration policy.
        int nPad = nLeftFlankBins + nRightFlankBins + 1;
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++)
            for (int j = 0; j < nRightFlankBins; ++j) scratchAtNode[ithDim][j] = esDsAtNode[ithDim][j];
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++)
            for (int j = (nXbins - nPad - nLeftFlankBins); j < (nXbins - nPad); ++j)
                scratchAtNode[ithDim][j] = esDsAtNode[ithDim][j];

        SSEUtils.convolveInPlaceSSTModalFftService(esDsAtNode, fftFYCA, scratchRA, nDimensionsE, nDimensionsD);

        int nItems2Copy = nXbins - nLeftFlankBins - nRightFlankBins;
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int i=0; i<nXbins; ++i) {
                if (esDsAtNode[ithDim][i] < 0.0 || i >= (nItems2Copy-1)) esDsAtNode[ithDim][i] = 0.0;
            }
        }

        // Restore the left nRightFlankBins and right nLeftFlankBins saved before convolution.
        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int j = 0; j < nRightFlankBins; ++j) {
                esDsAtNode[ithDim][j] = scratchAtNode[ithDim][j];
            }
        }
        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int j = (nXbins - nPad - nLeftFlankBins); j < (nXbins - nPad); ++j) {
                esDsAtNode[ithDim][j] = scratchAtNode[ithDim][j];
            }
        }
    }

    public static void propagateEandDinXQuaLikeSSTJavaFftService(double[][] esDsAtNode, double[][] fftEsDsAtNode, double[] fftFY, double[][] scratchAtNode, int nXbins, int nLeftFlankBins, int nRightFlankBins, int nDimensionsE, int nDimensionsD, JavaFftService ffts) {
        // Compatibility entry point for service-based callers; remove when they use ComplexFFT.
        propagateEandDinXQuaSSE(esDsAtNode, fftEsDsAtNode, fftFY, scratchAtNode, nXbins,
                nLeftFlankBins, nRightFlankBins, nDimensionsE, nDimensionsD, new SstFFT(nXbins, ffts));
    }

    // E/D and spectra alternate real/imaginary entries. Only transforms are delegated;
    // convolution, clipping and boundary restoration stay in Java for both implementations.
    public static void propagateEandDinXQuaSSE(double[][] esDsAtNode, double[][] fftEsDsAtNode,
            double[] fftFY, double[][] scratchAtNode, int nXbins, int nLeftFlankBins,
            int nRightFlankBins, int nDimensionsE, int nDimensionsD, ComplexFFT fft) {
        // Convolution reads input[i-offset]: preserve nRightFlankBins on the left and
        // nLeftFlankBins on the right. Keep this backend's existing E/D restoration policy.
        int nPad = nLeftFlankBins + nRightFlankBins + 1;
        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++)
            // note the additional i index here (as compared to the JTransforms version) to grab every other in esDsAtNode
            for (int j=0, i=0; j < nRightFlankBins; ++j, i+=2) scratchAtNode[ithDim][j] = esDsAtNode[ithDim][i];
        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++)
            // same as above with i index; note that we multiply by 2 because esDsAtNode is twice the length of nXbins because it's real complex real complex
            for (int j=(nXbins - nPad - nLeftFlankBins), i=2 * (nXbins - nPad - nLeftFlankBins); j < (nXbins - nPad); ++j, i+=2)
                scratchAtNode[ithDim][j] = esDsAtNode[ithDim][i];

        // debugging
        // System.out.println("After copying left and right flanks");
        // System.out.println("esDsAtNode = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode = " + Arrays.toString(scratchAtNode[1]));

        SSEUtils.convolveInPlace(esDsAtNode, fftEsDsAtNode, fftFY, nDimensionsE, nDimensionsD, fft);

        // number of real elements we will grab and move to the head of the array
        int nItems2Copy = nXbins - nLeftFlankBins - nRightFlankBins;

        // scale and move real elements in place to the head of the array, ignoring flanking elements
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            /*
             * Note that while we are grabbing every other element (so we can get just the real parts)
             * as we iterate over esDsAtNode, everyOtherToHeadInPlace moves these elements to the head
             * of esDsAtNode without any skipping, i.e.,
             * [ ... skipFirst_nLeftFlankBins ..., real1, real2, real3, ...]
             */
            // TODO: I believe we should not place the real elements consecutively in esDsAtNode
            // but rather keep them at every other index when using SST's JavaFftService; this class
            // requires alternating real elements any time FFT is done...
            //
            // this line below is not necessary (I added so I could pass unit tests, but instead
            // should update unit tests and remove this line)

            // everyOtherToHeadInPlace(esDsAtNode[ithDim], nXbins, nLeftFlankBins, nRightFlankBins, 2, 1.0); // grabbing real part and scaling by 1/nXbins

            for (int i=0; i<(nXbins*2); i += 2) {
                // if negative value, set to 0.0
                // if at last (nLeftFlankBins + nRightFlankBins) items, set to 0.0
                if (esDsAtNode[ithDim][i] < 0.0 || i >= ((nItems2Copy-1)*2)) {
                    esDsAtNode[ithDim][i] = 0.0;
                }
            }

            // passes tests in PropagatesQuaSSETest.java
//            for (int i=0; i<nXbins; ++i) {
//                // if negative value, set to 0.0
//                // if at last (nLeftFlankBins + nRightFlankBins) items, set to 0.0
//                if (esDsAtNode[ithDim][i] < 0.0 || i >= (nItems2Copy-1)) esDsAtNode[ithDim][i] = 0.0;
//            }
        }

        // Restore the left nRightFlankBins and right nLeftFlankBins saved before convolution.
        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int i=0, j=0; i < nRightFlankBins; ++i, j+=2) {
                esDsAtNode[ithDim][j] = scratchAtNode[ithDim][i];
            }
        }
        for (int ithDim=0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {
            for (int i=(nXbins - nPad - nLeftFlankBins), j=(nXbins - nPad - nLeftFlankBins)*2; i < (nXbins - nPad); ++i, j+=2) {
                    esDsAtNode[ithDim][j] = scratchAtNode[ithDim][i];
            }
        }

        // debugging
        // System.out.println("After scaling and zero-ing");
        // System.out.println("esDsAtNode (D's) = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode = " + Arrays.toString(scratchAtNode[1]));
    }

    /*
     * Carry out two FFTs, one on the Normal kernel, another on the 2D-array containing the E's and D's.
     * Then carry out inverse FFT on the resulting 2D-array, re-ordering elements so that all real elements occupy the first half.
     * Note: the input of FFT uses just the first half of all elements; the output interdigitates real and complex numbers;
     * the input of inverse-FFT is then the interdigitated array, and the output of inverse-FFT remains interdigitated.
     *
     * Function leaves results in scratch (after moving all real elements to the first half of each row of esDs)
     *
     * @param   esDsAtNode    2D-array containing E's followed by D's (e.g., for QuaSSE, esDs[0]=Es, esDs[1]=Ds)
     * @param   scratchAtNode 2D-array for storing/restoring flanking values and other math terms (from a node)
     * @param   fY  Normal kernel giving the density for changes in value of the quantitative trait
     * @param   nDimensionsE number of E equations (dimensions in plan) to solve for each quantitative ch
     * @param   nDimensionsD number of D equations (dimensions in plan) to solve for each quantitative ch
     * @param   fft instance of DoubleFFT_1D that will carry out FFT and inverse-FFT
     */
    public static void convolveInPlace(double[][] esDsAtNode, double[] fY, int nDimensionsE, int nDimensionsD, DoubleFFT_1D fft) {
        // doing E's and D's
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {

            // uncomment for testIntegrateOneBranchHiResOutsideClassJustX
            // System.out.println("Before FFT scratchAtNode[" + ithDim + "] = " + Arrays.toString(scratchAtNode[ithDim]));

            fft.realForwardFull(esDsAtNode[ithDim]); // FFT for each E and D dimension

            // uncomment for testIntegrateOneBranchHiResOutsideClassJustX
            // System.out.println("After FFT scratchAtNode[" + ithDim + "] = " + Arrays.toString(scratchAtNode[ithDim]));

            for (int i = 0; i < fY.length; i += 2) {
                // For (a + bi)(c + di), preserve a and b because the products overwrite
                // the input spectrum.
                double dataReal = esDsAtNode[ithDim][i];
                double dataImag = esDsAtNode[ithDim][i + 1];
                esDsAtNode[ithDim][i] = dataReal * fY[i] - dataImag * fY[i + 1];
                esDsAtNode[ithDim][i + 1] = dataReal * fY[i + 1] + dataImag * fY[i];
            }

            // uncomment for testIntegrateOneBranchHiResOutsideClassJustX
            // System.out.println("After FFT and * fY scratchAtNode[" + ithDim + "] = " + Arrays.toString(scratchAtNode[ithDim]));

            fft.complexInverse(esDsAtNode[ithDim], false); // inverse FFT for each E and D dimension

            // uncomment for testIntegrateOneBranchHiResOutsideClassJustX
            // System.out.println("After iFFT scratchAtNode[" + ithDim + "] = " + Arrays.toString(scratchAtNode[ithDim]));
        }

        // looking at things
        // System.out.println(Arrays.toString(scratchAtNode[0]));
        // System.out.println(Arrays.toString(scratchAtNode[1]));
    }

    /*
     * This convolution in place is done by calling methods of
     * RealArray (which invokes ModalFftService, a class that
     * initializes a JavaFftService instance)
     *
     * It works, but is a bit clunky. Below there is a cleaner
     * alternative that calls JavaFftService directly
     */
    public static void convolveInPlaceSSTModalFftService(double[][] esDsAtNode, ComplexArray fY, RealArray scratchRA, int nDimensionsE, int nDimensionsD) {
        int normalizingInverseFFTFactor = esDsAtNode[0].length;

        // doing E's and D's
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {

            int jthElem = 0;
            for (double v: esDsAtNode[ithDim]) {
                scratchRA.set(v, jthElem);

                jthElem++;
            }

            /*
             * We invoke the static ModalFftService (which implements interface FftService)
             * fft and ifft methods
             *
             * ModalFftService's fft and ifft methods can be accessed as members of ArrayBase,
             * called from within AbstractArray, which in turn is the superclass of AbstractComplexArray
             * -- created here with the .tocre() call
             *
             * Once the fft result is returned, we do element-wise multiplication with eMul,
             * inverse-fft, get the real part, and return the values
             */
            esDsAtNode[ithDim] = scratchRA.tocRe().fft().eMul(fY).ifft().torRe().values(); // already comes out normalized (unlike R version and JTransforms)
        }
    }

    /*
     * Retain JavaFftService callers while sharing the implementation-independent convolution below.
     */
    public static void convolveInPlaceSSTJavaFftService(double[][] esDsAtNode, double[][] fftEsDsAtNode, double[] fftFY, int nDimensionsE, int nDimensionsD, int[] nXbins, JavaFftService ffts) {
        // Compatibility entry point for service-based callers; remove when they use ComplexFFT.
        convolveInPlace(esDsAtNode, fftEsDsAtNode, fftFY, nDimensionsE, nDimensionsD,
                new SstFFT(nXbins, ffts));
    }

    // Transform each E/D array, multiply by the supplied kernel spectrum, then invert.
    // Arrays alternate real/imaginary entries; inverse() supplies 1/N, so do not scale again here.
    public static void convolveInPlace(double[][] esDsAtNode, double[][] fftEsDsAtNode, double[] fftFY,
            int nDimensionsE, int nDimensionsD, ComplexFFT fft) {
        // int normalizingInverseFFTFactor = esDsAtNode[0].length;

        // doing E's and D's
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {

            // fft-ing
            fft.forward(esDsAtNode[ithDim], fftEsDsAtNode[ithDim]); // real complex real complex...
            // kylieJNI.fft(nRE, input_array_REdoubles_with_twice_nRE, output_array_REcomplex_with_twice_nRE)

            // convolving
            for (int i=0; i < fftFY.length; i += 2) {
                // For (a + bi)(c + di), preserve a and b because the products overwrite
                // the input spectrum.
                double dataReal = fftEsDsAtNode[ithDim][i];
                double dataImag = fftEsDsAtNode[ithDim][i + 1];
                fftEsDsAtNode[ithDim][i] = dataReal * fftFY[i] - dataImag * fftFY[i + 1];
                fftEsDsAtNode[ithDim][i + 1] = dataReal * fftFY[i + 1] + dataImag * fftFY[i];
            }

            // ifft-ing
            fft.inverse(fftEsDsAtNode[ithDim], esDsAtNode[ithDim]);
            // kylieJNI.ifft(nRE, input_array_REcomplex_with_twice_nRE, output_array_REdoubles_with_twice_nRE)
        }
    }

    /*
     * Builds Normal distribution (yValues) where x is the CHANGE in quantitative trait values
     * (not the quantitative trait values themselves!). Under BM, for example, this
     * kernel is centered at 0.0
     *
     * Note: this kernel (yValues) is later FFT-ed in the likelihood class, and then used in the
     * convolution function, and for reasons I do not fully understand, this bell-shaped kernel\
     * needs to be cuts in two, with the left half being placed at the (right-)tail end of the kernel,
     * and the right half at the (left-)head of the kernel
     *
     * 'nLeftFlankBins' and 'nRightFlankBins' are also the number of bins on the left-
     * side (and right-side, after skipping (nLeftFlankBins + nRightFlankBins))
     * of E and D that are not updated by 'propagateEandDinXQuaLike'
     *
     * To be used with values distributed over bins as required by JTransforms or ModalFftService,
     * that is, elements should be placed consecutively
     *
     * @param   yValues (= fY) where the result is left; gives the probability density of a given change in quantitative trait value
     * @param   mean    (= changeInXNormalMean = diversitree's drift * -dt) is the mean expected change in quantitative trait value
     * @param   sd      (= changeInXNormalSd = squared root(diversitree's diffusion * dt)) is the standard deviation of the expected change in quantitative trait value
     * @param   nXbins  total number of bins resulting from discretizing quantitative trait-change normal kernel (fY and each row of esDs will have these many nXbins)
     * @param   nLeftFlankBins   how many bins on the right side of kernel are non-zero
     * @param   nRightFlankBins  how many bins on the left side of kernel are non-zero
     * @param   dx               size of each bin
     */
    public static void makeNormalKernelInPlace(double[] yValues, double mean, double sd, int nXbins, int nLeftFlankBins, int nRightFlankBins, double dx) {
        double total = 0.0;

        double x = 0.0;
        for (int i = 0; i <= nRightFlankBins; i++) {
            yValues[i] = getNormalDensity(x, mean, sd); // in diversitree's C code, this is further multiplied by dx (think this is unnecessary, b/c it doesn't change the result!)
            total += yValues[i];
            x += dx;
        }

        for (int i = nRightFlankBins + 1; i < nXbins - nLeftFlankBins; i++) {
            yValues[i] = 0;
        }

        x = -nLeftFlankBins * dx;
        for (int i = nXbins - nLeftFlankBins; i < nXbins; i++) {
            yValues[i] = getNormalDensity(x, mean, sd);
            total += yValues[i];
            x += dx;
        }

        if (!Double.isFinite(total) || total <= 0)
            throw new QuaSSEKernelException("Gaussian kernel normalization total must be positive and finite: " + total);
        for (int i = 0; i <= nRightFlankBins; i++) yValues[i] /= total;
        for (int i = (nXbins - nLeftFlankBins); i < nXbins; i++) yValues[i] /= total;
    }

    /*
     * Builds Normal distribution (yValues) where x is the CHANGE in quantitative trait values
     * (not the quantitative trait values themselves!). Under BM, for example, this
     * kernel is centered at 0.0
     *
     * Note: this kernel (yValues) is later FFT-ed in the likelihood class, and then used in the
     * convolution function, and for reasons I do not fully understand, this bell-shaped kernel\
     * needs to be cuts in two, with the left half being placed at the (right-)tail end of the kernel,
     * and the right half at the (left-)head of the kernel
     *
     * 'nLeftFlankBins' and 'nRightFlankBins' are also the number of bins on the left-
     * side (and right-side, after skipping (nLeftFlankBins + nRightFlankBins))
     * of E and D that are not updated by 'propagateEandDinXQuaLike'
     *
     * To be used with values distributed over bins as required by SST's JavaFftService,
     * that is, elements should be alternated, skipping even indices (that's why we multiply the
     * index integers in loops by two)
     *
     * @param   yValues (= fY) where the result is left; gives the probability density of a given change in quantitative trait value
     * @param   mean    (= changeInXNormalMean = diversitree's drift * -dt) is the mean expected change in quantitative trait value
     * @param   sd      (= changeInXNormalSd = squared root(diversitree's diffusion * dt)) is the standard deviation of the expected change in quantitative trait value
     * @param   nXbins  total number of bins resulting from discretizing quantitative trait-change normal kernel (fY and each row of esDs will have these many nXbins)
     * @param   nLeftFlankBins   how many bins on the right side of kernel are non-zero
     * @param   nRightFlankBins  how many bins on the left side of kernel are non-zero
     * @param   dx               size of each bin
     */
    public static void makeNormalKernelInPlaceSSTJavaFftService(double[] yValues, double mean, double sd, int nXbins, int nLeftFlankBins, int nRightFlankBins, double dx) {
        double total = 0.0;

        double x = 0.0;
        for (int i=0; i <= nRightFlankBins*2; i+=2) {
            yValues[i] = getNormalDensity(x, mean, sd); // in diversitree's C code, this is further multiplied by dx (think this is unnecessary, b/c it doesn't change the result!)
            total += yValues[i];
            x += dx;
        }

        for (int i=(nRightFlankBins+1)*2; i < (nXbins - nLeftFlankBins)*2; i+=2) {
            yValues[i] = 0;
        }

        x = -nLeftFlankBins * dx;
        for (int i=(nXbins - nLeftFlankBins)*2; i < nXbins*2; i+=2) {
            yValues[i] = getNormalDensity(x, mean, sd);
            total += yValues[i];
            x += dx;
        }

        if (!Double.isFinite(total) || total <= 0)
            throw new QuaSSEKernelException("Gaussian kernel normalization total must be positive and finite: " + total);
        for (int i=0; i <= nRightFlankBins*2; i+=2) yValues[i] /= total;
        for (int i=(nXbins - nLeftFlankBins)*2; i < nXbins*2; i+=2) yValues[i] /= total;
    }

    /*
     * (DEPRECATED: will be done in makeNormalKernelInPlace)
     */
    public static void normalizeArray(double[] doubleArray) {
        double theSum = 0.0;
        for (int i = 0; i < doubleArray.length; i++) {
            theSum += doubleArray[i];
        }

        for (int i = 0; i < doubleArray.length; i++) {
            doubleArray[i] = doubleArray[i] / theSum;
        }
    }

    /*
     * Return equivalent to R's dnorm(x, mean, sd)
     *
     * @param   x   random variable value
     * @param   mean    Mean of normal distribution
     * @param   sd  Standard deviation of normal distribution
     * @return  density of x under Normal distribution
     */
    public static double getNormalDensity(double x, double mean, double sd) {
        double x0 = x - mean;
        return Math.exp(-x0 * x0 / (2 * sd * sd)) / (sd * SQRT2PI);
    }

    /*
     * Move 'everyOther' items of an array to its head, skipping the first 'skipFirstN'
     * and the last 'skipLastN'
     *
     * @param   anArray source and destination array
     * @param   nXbins  total number of bins resulting from discretizing quantitative trait-change normal kernel (fY and each row of esDs will have these many nXbins)
     * @param   skipFirstN  number of discrete quantitative trait bins on the left-side of E and D that are not affected by 'propagateEandDinXQuaLike'
     * @param   skipLastN  number of discrete quantitative trait bins on the right-side of E and D that are not affected by 'propagateEandDinXQuaLike'
     * @param   everyOther  every other 'everyOther' elements (skipping the first skipFirstN and the last skipLastN) will be moved to the head of the array
     * @param   scaleBy will scale every other element by this
     */
    public static void everyOtherToHeadInPlace(double[] anArray, int nXbins, int skipFirstN, int skipLastN, int everyOther, double scaleBy) {
        //if (skipFirstN == 0 && skipLastN == 0) skipLastN = -1; // this should only happen in debugging, where there are no elements to skip at the start or end
        for (int i=skipFirstN * 2, j=skipFirstN; i <= (anArray.length-2); i+=everyOther, j++) {
            // first implementation if (j < (nXbins - skipLastN - (skipFirstN + skipLastN) - 1)) {
            if (j <= (nXbins - skipLastN - (skipFirstN + skipLastN) - 1)) {
                anArray[j] = anArray[i] * scaleBy;
            }
        }
    }

    /*
     * Simple function to take the first half of a double array and spread it
     * through the entire array, skipping every other element (setting
     * those to 0.0).
     */
    public static void everyOtherExpandInPlace(double[] anArray) {
        for (int i=anArray.length/2-1, j=anArray.length-2; i>=0; i--, j-=2) {
            anArray[j] = anArray[i];
            anArray[j +1 ] = 0.0;
        }
    }

    public static void hiToLoTransferInPlace(double[] fromArray, double[] toArray, int[] idxs4Transfer, boolean jtransforms) {
        int i = 0;
        for (int j: idxs4Transfer) {
            toArray[i] = fromArray[j];

            if (!jtransforms) i++; // add extra 1, so Ds are interdigitated for SST's JavaFftService
            i++;
        }
    }

    public static double calculateNormalizationFactorJTransforms(double[] dS, double binSize) {
        double normalizationFactorFromDs = 0.0;
        for (int i=0; i<dS.length/2; ++i) normalizationFactorFromDs += dS[i]; // dS is 2 * as large as it should be because it needs to hold complex parts resulting from FFTs

        // debugging
        // System.out.println("normalizationFactorFromDs inside calculateNormalizationFactor = " + normalizationFactorFromDs);
        // System.out.println("binSize inside calculateNormalizationFactor = " + binSize);

        normalizationFactorFromDs *= binSize;

        // debugging
        // System.out.println("normalizationFactorFromDs * binSize = " + normalizationFactorFromDs);

        return normalizationFactorFromDs;
    }

    public static double calculateNormalizationFactor(double[] dS, double binSize) {
        double normalizationFactorFromDs = 0.0;
        for (int i=0; i<dS.length; i += 2) normalizationFactorFromDs += dS[i]; // dS is 2 * as large as it should be because it needs to hold complex parts resulting from FFTs

        // debugging
        // System.out.println("normalizationFactorFromDs inside calculateNormalizationFactor = " + normalizationFactorFromDs);
        // System.out.println("binSize inside calculateNormalizationFactor = " + binSize);

        normalizationFactorFromDs *= binSize;

        // debugging
        // System.out.println("normalizationFactorFromDs * binSize = " + normalizationFactorFromDs);

        return normalizationFactorFromDs;
    }
}
