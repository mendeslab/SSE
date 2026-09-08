package SSE;

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

    /*
     * Propagates E's and D's, leaving the result in esDs.
     * Equivalent to Fitzjohn's fftR.propagate.t (R/model-quasse-fftR.R) and propagate_t (src/quasse-eqs-fftC.c)
     *
     * Note: only the first nUsefulTraitBins elements in each row of esDs are used (and in birthRate and deathRate);
     * This method is to be used whenever FFTing when propagating over x is carried out with JTransforms
     *
     * @param   esDs    2D-array containing E's followed by D's (e.g., for QuaSSE, esDs[0]=Es, esDs[1]=Ds)
     * @param   scratchAtNode 2D-array for storing/restoring flanking values and other math terms
     * @param   birthRate   birth rates (lambdas), one per quantitative ch useful bin
     * @param   deathRate   death rates (mus), one per quantitative ch useful bin
     * @param   dt  length of time slice over which to propagate
     * @param   nUsefulTraitBins number of useful bins (i.e., excluding the # of flanking bins = npad = nLeftFlankBins + nRightFlankbins) in discretized pdf over change in quantitative ch values
     * @param   nDimensionsD number of D equations (dimensions in plan) to solve for each quantitative ch
     * @return  no return, leaves result in 'esDs' and 'scratch'
     */
    public static void propagateEandDinTQuaSSEInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double[] birthRate, double[] deathRate, double dt, int nUsefulTraitBins, int nDimensionsD) {
        // iterating over all continuous character bins (total # = nx), to update E and D for each of those bins
        for (int i=0; i < nUsefulTraitBins; ++i) {
            double ithLambda = birthRate[i];
            double ithMu = deathRate[i];
            double netDivRate = ithLambda - ithMu;
            double ithZ = Math.exp(dt * netDivRate);
            double ithE = esDsAtNode[0][i]; // Ask Xia: note that because i starts at 0, we're grabbing E values in left-padding

            double tmp1 = ithMu - ithLambda * ithE;
            double tmp2 = ithZ * (ithE - 1);

            /* Updating E's */
            esDsAtNode[0][i] = (tmp1 + tmp2 * ithMu) / (tmp1 + tmp2 * ithLambda);

            // checking against R
            // System.out.println("ithMu = " + ithMu + " ithLambda = " + ithLambda + " ithE = " + ithE + " ithZ = " + ithZ);
            // System.out.println("esDsAtNode[0][i] = " + esDsAtNode[0][i]);

            // Saving values for updating D below
            tmp1 = (ithLambda - ithMu) / (ithZ * ithLambda - ithMu + (1 - ithZ) * ithLambda * ithE);

            // checking against R
            // System.out.println("i = " + i + " numerator (z * r * r) = " + (ithZ * Math.pow((ithLambda - ithMu), 2)));
            // System.out.println("i = " + i + " denominator (z * lambda - mu + (1-z) * lambda * e0) = " + (ithZ * ithLambda - ithMu + (1 - ithZ) * ithLambda * ithE));

            scratchAtNode[1][i] = ithZ * tmp1 * tmp1;

            // checking against R
            // System.out.println("i = " + i + " (numerator/denominator^2) = " + scratchAtNode[1][i]);
        }

        /* Updating D's */
        // iterating over dimensions in plan to transform (total # = nd = number of equations for D, say, 4 if A C G and T)
        // one quantitative ch: nd = 2 for QuaSSE (one eq for E, one for D), nd=5 for MoSSE (4 eqs for D, one for E)
        // ithDim starts at 1 because esDs[0] are the E's

        for (int ithDim=1; ithDim <= nDimensionsD; ithDim++) {
            // int ithDimStartIdx = nUsefulTraitBins * ithDim; // skipping E's

            // iterating over bins of this dimension
            for (int j=0; j < nUsefulTraitBins; j++) {
                if (esDsAtNode[ithDim][j] < 0) esDsAtNode[ithDim][j] = 0;
                else {
                    // checking against R
                    // System.out.println("Scr j=" + j + " scratchAtNode[0][j] = " + scratchAtNode[0][j] + " scratchAtNode[1][j] = " + scratchAtNode[1][j]);
                    // System.out.println("Bef j=" + j + " esDsAtNode[0][j] = " + esDsAtNode[0][j] + " esDsAtNode[1][j] = " + esDsAtNode[1][j]);
                    esDsAtNode[ithDim][j] *= scratchAtNode[1][j];

                    // checking against R
                    // System.out.println("Aft j=" + j + " esDsAtNode[0][j] = " + esDsAtNode[0][j] + " esDsAtNode[1][j] = " + esDsAtNode[1][j]);
                }
            }
        }
    }

    /*
     * Propagates E's and D's, leaving the result in esDs.
     * Equivalent to Fitzjohn's fftR.propagate.t (R/model-quasse-fftR.R) and propagate_t (src/quasse-eqs-fftC.c)
     *
     * Note: only the first nUsefulTraitBins elements in each row of esDs are used (and in birthRate and deathRate);
     * This method is to be used whenever FFTing when propagating over x is carried out with SST's library
     * class JavaFftService (which requires that elements are placed on every other index of the input array)
     *
     * @param   esDs    2D-array containing E's followed by D's (e.g., for QuaSSE, esDs[0]=Es, esDs[1]=Ds)
     * @param   scratchAtNode 2D-array for storing/restoring flanking values and other math terms
     * @param   birthRate   birth rates (lambdas), one per quantitative ch useful bin
     * @param   deathRate   death rates (mus), one per quantitative ch useful bin
     * @param   dt  length of time slice over which to propagate
     * @param   nUsefulTraitBins number of useful bins (i.e., excluding the # of flanking bins = npad = nLeftFlankBins + nRightFlankbins) in discretized pdf over change in quantitative ch values
     * @param   nDimensionsD number of D equations (dimensions in plan) to solve for each quantitative ch
     * @return  no return, leaves result in 'esDs' and 'scratch'
     */
    public static void propagateEandDinTQuaSSEInPlaceSSTJavaFftService(double[][] esDsAtNode, double[][] scratchAtNode, double[] birthRate, double[] deathRate, double dt, int nUsefulTraitBins, int nDimensionsD) {
        // iterating over all continuous character bins (total # = nx), to update E and D for each of those bins
        for (int i=0, j=0; i < nUsefulTraitBins; ++i, j+=2) {
            double ithLambda = birthRate[i];
            double ithMu = deathRate[i];
            double netDivRate = ithLambda - ithMu;
            double ithZ = Math.exp(dt * netDivRate);
            double ithE = esDsAtNode[0][j]; // Ask Xia: note that because i starts at 0, we're grabbing E values in left-padding

            double tmp1 = ithMu - ithLambda * ithE;
            double tmp2 = ithZ * (ithE - 1);

            /* Updating E's */
            esDsAtNode[0][j] = (tmp1 + tmp2 * ithMu) / (tmp1 + tmp2 * ithLambda);

            // checking against R
            // System.out.println("ithMu = " + ithMu + " ithLambda = " + ithLambda + " ithE = " + ithE + " ithZ = " + ithZ);
            // System.out.println("esDsAtNode[0][i] = " + esDsAtNode[0][i]);

            // Saving values for updating D below
            tmp1 = (ithLambda - ithMu) / (ithZ * ithLambda - ithMu + (1 - ithZ) * ithLambda * ithE);

            // checking against R
            // System.out.println("i = " + i + " numerator (z * r * r) = " + (ithZ * Math.pow((ithLambda - ithMu), 2)));
            // System.out.println("i = " + i + " denominator (z * lambda - mu + (1-z) * lambda * e0) = " + (ithZ * ithLambda - ithMu + (1 - ithZ) * ithLambda * ithE));

            // NOTE: we're not skipping even-numbered indices when just caching things here in scratch!
            // (so below, we use index i instead of j)
            scratchAtNode[1][i] = ithZ * tmp1 * tmp1;

            // checking against R
            // System.out.println("i = " + i + " (numerator/denominator^2) = " + scratchAtNode[1][i]);
        }

        /* Updating D's */
        // iterating over dimensions in plan to transform (total # = nd = number of equations for D, say, 4 if A C G and T)
        // one quantitative ch: nd = 2 for QuaSSE (one eq for E, one for D), nd=5 for MoSSE (4 eqs for D, one for E)
        // ithDim starts at 1 because esDs[0] are the E's

        for (int ithDim=1; ithDim <= nDimensionsD; ithDim++) {
            // int ithDimStartIdx = nUsefulTraitBins * ithDim; // skipping E's

            // iterating over bins of this dimension
            for (int i=0, j=0; i < nUsefulTraitBins; i++, j+=2) {
                if (esDsAtNode[ithDim][j] < 0) esDsAtNode[ithDim][j] = 0;
                else {
                    // checking against R
                    // System.out.println("Scr j=" + j + " scratchAtNode[0][j] = " + scratchAtNode[0][j] + " scratchAtNode[1][j] = " + scratchAtNode[1][j]);
                    // System.out.println("Bef j=" + j + " esDsAtNode[0][j] = " + esDsAtNode[0][j] + " esDsAtNode[1][j] = " + esDsAtNode[1][j]);
                    esDsAtNode[ithDim][j] *= scratchAtNode[1][i]; // have to use scratch's i-th element, not j-th!

                    // checking against R
                    // System.out.println("Aft j=" + j + " esDsAtNode[0][j] = " + esDsAtNode[0][j] + " esDsAtNode[1][j] = " + esDsAtNode[1][j]);
                }
            }
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

        int[] nDims = new int[] { nXbins };
        SSEUtils.convolveInPlaceSSTJavaFftService(esDsAtNode, fftEsDsAtNode, fftFY, nDimensionsE, nDimensionsD, nDims, ffts);

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
     * This convolution in place is done by directly calling methods
     * of the JavaFftService class
     */
    public static void convolveInPlaceSSTJavaFftService(double[][] esDsAtNode, double[][] fftEsDsAtNode, double[] fftFY, int nDimensionsE, int nDimensionsD, int[] nXbins, JavaFftService ffts) {
        // int normalizingInverseFFTFactor = esDsAtNode[0].length;

        // doing E's and D's
        for (int ithDim = 0; ithDim < (nDimensionsE + nDimensionsD); ithDim++) {

            // fft-ing
            ffts.fft(nXbins, esDsAtNode[ithDim], fftEsDsAtNode[ithDim]); // fftEsDsAtNode: real complex real complex...
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
            ffts.ifft(nXbins, fftEsDsAtNode[ithDim], esDsAtNode[ithDim]);
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
