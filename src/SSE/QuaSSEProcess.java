package SSE;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import org.jtransforms.fft.DoubleFFT_1D;
import org.shared.array.ComplexArray;
import org.shared.array.RealArray;
import org.shared.fft.JavaFftService;

import java.util.Arrays;

@Description("Specifies a quantitative trait(s) state-dependent speciation and" +
        "extinction birth-death process.")
public abstract class QuaSSEProcess extends Distribution {

    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);
    final public Input<RealParameter> dtMaxInput = new Input<>("dtMax", "Length of max time interval over which integration is carried out.", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> dynamicDtInput = new Input<>("dynDt", "If interval over which to carry out integration should be dynamically adjusted to maximize accuracy.", Input.Validate.REQUIRED);
    final public Input<RealParameter> tcInput = new Input<>("tc", "Time (backwards, i.e., present=0.0) when integration happens with discretization at low resolution.", Input.Validate.OPTIONAL);
    final public Input<IntegerParameter> nXbinsInput = new Input<>("nX", "Total number of quantitative trait bins after discretization at low resolution.", Input.Validate.REQUIRED);
    final public Input<RealParameter> dXBinInput = new Input<>("dX", "Width of quantitative trait bins at low resolution.", Input.Validate.REQUIRED);
    final public Input<RealParameter> xMidInput = new Input<>("xMid", "Midpoint to center the quantitative trait space.", Input.Validate.REQUIRED);
    final public Input<RealParameter> driftInput = new Input<>("drift", "Drift term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> diffusionInput = new Input<>("diffusion", "Diffusion term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> flankWidthScalerInput = new Input<>("flankWidthScaler", "Multiplier of normal standard deviation when determining number of flanking bins.", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> highLowRatioInput = new Input<>("hiLoRatio", "Scale nX by this when at high resolution.", Input.Validate.REQUIRED);
    final public Input<String> priorProbAtRootTypeInput = new Input<>("priorProbAtRootType", "Type of root prior probabilities for D's.", Input.Validate.REQUIRED);
    final public Input<RealParameter> priorProbsAtRootInput = new Input<>("givenPriorProbsAtRoot", "Root prior probabilities for D's at high resolution.", Input.Validate.XOR, priorProbAtRootTypeInput);

    protected Tree tree;
    protected RealParameter quTraits;

    // dealing with prior probability at root
    protected static final String FLAT = "Flat";
    protected static final String OBS = "Observed";
    protected double[] priorProbsAtRoot;
    protected String rootPriorType;

    // state for dimensioning things and setting up resolution of integration
    protected boolean got2LowRes; // if the pruning got to do low-resolution integrated (set in startRecursionAtRootNode())
    protected boolean providedPriorAtRoot = false;
    protected boolean dynamicallyAdjustDt;
    protected double dtMax, tc;
    protected double dXbin, flankWidthScaler, xMinLo, xMinHi, xMid;
    protected int nXbinsLo, nUsefulXbinsLo, nXbinsHi, nUsefulXbinsHi, hiLoRatio;
    protected int[] nXbinsLoSST, nXbinsHiSST; // for SST JavaFftService, same as nXbinsLo and nXbinsHi
    protected int[] nUsefulXbins, nLeftNRightFlanksHi, nLeftNRightFlanksLo;
    protected double[] xLo, xHi; // x rulers
    protected int[] hiLoIdxs4Transfer;

    // quantitative trait evolution
    protected double drift;
    protected double diffusion;
    protected double changeInXNormalMean; // (=diversitree's drift)
    protected double changeInXNormalSd; // (=diversitree's diffusion)
    protected double[] fYLo, fYHi, fftFYLo, fftFYHi;
    protected DoubleFFT_1D fftForEandDLo, fftForEandDHi;
    // using SST library
    JavaFftService jffts;

    @Override
    public void initAndValidate() {

        tree = treeInput.get();
        rootPriorType = priorProbAtRootTypeInput.get();

        if (priorProbsAtRootInput.get() != null) {
            throw new IllegalArgumentException("givenPriorProbsAtRoot is not supported: a supplied array "
                    + "has no defined mapping when the QuaSSE grid changes. Use Flat or Observed.");
        }
        if (!FLAT.equals(rootPriorType) && !OBS.equals(rootPriorType))
            throw new IllegalArgumentException("priorProbAtRootType must be Flat or Observed.");
        // quTraits = quTraitsInput.get();

        nLeftNRightFlanksLo = new int[2];
        nLeftNRightFlanksHi = new int[2];
        nUsefulXbins = new int[2];

        dynamicallyAdjustDt = dynamicDtInput.get().getValue();
        dtMax = dtMaxInput.get().getValue();
        dXbin = dXBinInput.get().getValue();
        tc = tcInput.get().getValue();

        hiLoRatio = highLowRatioInput.get().getValue();
        nXbinsLo = nXbinsInput.get().getValue();
        if (nXbinsLo <= 0 || (nXbinsLo & (nXbinsLo-1)) != 0)
            throw new IllegalArgumentException("Number of quantitative character bins must be a power of 2. It was " + nXbinsLo);
        if (hiLoRatio < 1 || 2L * nXbinsLo * hiLoRatio > Integer.MAX_VALUE)
            throw new IllegalArgumentException("QuaSSE hiLoRatio must be positive and FFT array sizes must fit in an int.");
        nXbinsHi = nXbinsLo * hiLoRatio;

        xMid = xMidInput.get().getValue();
        flankWidthScaler = flankWidthScalerInput.get().getValue();
        if (!Double.isFinite(dtMax) || dtMax <= 0 || !Double.isFinite(dXbin) || dXbin <= 0
                || !Double.isFinite(flankWidthScaler) || flankWidthScaler <= 0
                || !Double.isFinite(xMid) || !Double.isFinite(tc) || tc < 0)
            throw new IllegalArgumentException("QuaSSE requires finite positive dtMax, dX and flankWidthScaler, "
                    + "finite xMid and finite nonnegative tc.");
        drift = driftInput.get().getValue();
        diffusion = diffusionInput.get().getValue();
        changeInXNormalMean = drift * -dtMax;
        changeInXNormalSd = Math.sqrt(diffusion * dtMax);

        prepareDimensionsInPlace(); // in parent class
        prepareXRulers(); // in parent class

        // getting indices for transferring high- to low-res E's and D's during pruning
        hiLoIdxs4Transfer = new int[nUsefulXbinsLo];
        boolean jtransforms = false;
        populateIndicesHiLo(hiLoIdxs4Transfer, hiLoRatio, nUsefulXbinsLo, jtransforms); // populates hiLoIdxs4Transfer

        // debugging
        // System.out.println("hiLoIdxs4Transfer = " + Arrays.toString(hiLoIdxs4Transfer));

        // JTransforms version
        fftForEandDLo = new DoubleFFT_1D(nXbinsLo);
        fftForEandDHi = new DoubleFFT_1D(nXbinsHi);
//        fYLo = new double[2 * nXbinsLo]; // for real and complex part after FFT
//        fYHi = new double[2 * nXbinsHi]; // for real and complex part after FFT

        // SST
        fYLo = new double[nXbinsLo * 2]; // just real
        fYHi = new double[nXbinsHi * 2]; // just real
        fftFYLo = new double[nXbinsLo * 2]; // just real
        fftFYHi = new double[nXbinsHi * 2]; // just real
        nXbinsLoSST = new int[] { nXbinsLo };
        nXbinsHiSST = new int[] { nXbinsHi };
        jffts = new JavaFftService();

        // populatefY(dtMax, true, false, true, true); // force populate fY, and do FFT
        // populatefY(dtMax, true, false, true, false); // force populate fY, and do FFT
    }

    /*
     * Set up the number of left and right discrete bins at both low and high resolution
     * that will not contribute to... (fill this out later)
     */
    protected void prepareDimensionsInPlace() {
        int[] padding = calculatePadding(drift, diffusion);
        if (padding == null)
            throw new IllegalArgumentException("Invalid QuaSSE drift/diffusion or insufficient grid for kernel support.");
        nLeftNRightFlanksLo[0] = padding[0];
        nLeftNRightFlanksHi[0] = hiLoRatio * nLeftNRightFlanksLo[0];

        nLeftNRightFlanksLo[1] = padding[1];
        nLeftNRightFlanksHi[1] = hiLoRatio * nLeftNRightFlanksLo[1];

        // nXbinsHi = hiLoRatio * nXbinsLo;

        nUsefulXbinsLo = nXbinsLo - (nLeftNRightFlanksLo[0] + 1 + nLeftNRightFlanksLo[1]);
        nUsefulXbinsHi = nXbinsHi - (nLeftNRightFlanksHi[0] + 1 + nLeftNRightFlanksHi[1]);
        nUsefulXbins[0] = nUsefulXbinsLo;
        nUsefulXbins[1] = nUsefulXbinsHi;

        // debugging
        // System.out.println("\n\nSetting dimensions of QuaSSEProcess");
        // System.out.println("changeInXNormalMean = " + changeInXNormalMean + " changeInXNormalSd = " + changeInXNormalSd);
        // System.out.println("flankWidthScaler = " + flankWidthScaler + " dXbin = " + dXbin);
        // System.out.println("nLeftNRightFlanksLo[0] = " + nLeftNRightFlanksLo[0]);
        // System.out.println("nXbinsLo = " + nXbinsLo + " nXbinsHi = " + nXbinsHi);
        // System.out.println("nLeftFlankLo = " + nLeftNRightFlanksLo[0] + " nRightFlankLo = " + nLeftNRightFlanksLo[1]);
        // System.out.println("nUsefulXbinsLo = " + nUsefulXbinsLo);
        // System.out.println("nLeftFlankHi = " + nLeftNRightFlanksHi[0] + " nRightFlankHi = " + nLeftNRightFlanksHi[1]);
        // System.out.println("nUsefulXbinsHi = " + nUsefulXbinsHi);
    }

    // Calculate candidate support without mutating the live grid. Null denotes an invalid proposal.
    protected int[] calculatePadding(double candidateDrift, double candidateDiffusion) {
        if (!Double.isFinite(candidateDrift) || !Double.isFinite(candidateDiffusion) || candidateDiffusion <= 0)
            return null;
        // Backward support extents are ±v*t+a*sqrt(t), a=w*sqrt(diffusion). The growing side
        // peaks at T=dtMax; the opposing side peaks at min(T,(a/(2*abs(v)))^2).
        // Round only after maximising over all steps, then scale the coarse counts for the fine grid.
        double a = flankWidthScaler * Math.sqrt(candidateDiffusion);
        double speed = Math.abs(candidateDrift);
        double peakRoot = speed == 0 ? Math.sqrt(dtMax) : Math.min(Math.sqrt(dtMax), a / (2 * speed));
        double growing = speed * dtMax + a * Math.sqrt(dtMax);
        double opposing = peakRoot * (a - speed * peakRoot);
        double left = Math.ceil((candidateDrift >= 0 ? growing : opposing) / dXbin);
        double right = Math.ceil((candidateDrift >= 0 ? opposing : growing) / dXbin);
        // The preserved boundary strips must not overlap: useful > left+right.
        // This also bounds integer casts, scaled counts and the fine-to-coarse map before allocation.
        if (!Double.isFinite(left) || !Double.isFinite(right) || left < 1 || right < 1
                || 2 * (left + right) + 1 >= nXbinsLo)
            return null;
        int useful = nXbinsLo - ((int) left + (int) right + 1);
        double minLo = xMid - dXbin * Math.ceil((useful - 1.0) / 2.0);
        double minHi = minLo - dXbin * (1.0 - 1.0 / hiLoRatio);
        if (dXbin / hiLoRatio == 0 || !Double.isFinite(minLo) || !Double.isFinite(minHi)
                || !Double.isFinite(minLo + (useful - 1.0) * dXbin)
                || !Double.isFinite(minLo + (useful - 1.0 / hiLoRatio) * dXbin))
            return null;
        return new int[] {(int) left, (int) right};
    }

    /*
     *
     */
    protected void prepareXRulers() {

        xMinLo = xMid - dXbin * Math.ceil((nUsefulXbinsLo - 1.0) / 2.0);
        xMinHi = xMinLo - dXbin * (1.0 - 1.0 / hiLoRatio);

        // debugging
        // System.out.println("xMinLo = " + xMinLo + " xMinHi = " + xMinHi);

        // preparing x rulers
        xLo = new double[nUsefulXbinsLo];
        xLo[0] = xMinLo;
        for (int i = 1; i < nUsefulXbinsLo; i++) {
            // System.out.println("xLo[" + (i-1) + "]" + xLo[i-1]);
            xLo[i] = xLo[i-1] + dXbin;
            // xLo[i] = Math.round((xLo[i-1] + dXbin) * 1e4) / 1e4;
            // System.out.println("xLo[" + i + "] = " + xLo[i] + " dx = " + dXbin);
        }

        xHi = new double[nUsefulXbinsHi];
        xHi[0] = xMinHi;
        for (int i = 1; i< nUsefulXbinsHi; i++) {
            // System.out.println("xHi[" + (i-1) + "]" + xHi[i-1]);
            xHi[i] = xHi[i-1] + dXbin / hiLoRatio;
            // xHi[i] = Math.round((xHi[i-1] + (dXbin / hiLoRatio)) * 1e4) / 1e4;
            // System.out.println("xHi[" + i + "] = " + xHi[i] + " dx = " + dXbin);
        }
    }

    /*
     *
     */
    protected abstract void populateMacroevolParams(boolean ignoreRefresh);

    /*
     *
     */
    protected void populatefY(double aDt, boolean forceRecalcKernel, boolean dtChanged, boolean doFFT, boolean lowRes, boolean jtransforms) {
        // finding out if we need to recalculate fY or not
        boolean didRefreshNowMustRecalcKernel = false;
        if (dtChanged || forceRecalcKernel) {
            changeInXNormalMean = drift * -aDt;
            didRefreshNowMustRecalcKernel = true;
        }
        if (dtChanged || forceRecalcKernel) {
            changeInXNormalSd = Math.sqrt(diffusion * aDt);
            didRefreshNowMustRecalcKernel = true;
        }

        if (forceRecalcKernel || didRefreshNowMustRecalcKernel) {
        	
        	// System.out.println("changeInXNormalMean = " + changeInXNormalMean);
        	// System.out.println("changeInXNormalSd = " + changeInXNormalSd);
        	
            if (lowRes) {
                if (jtransforms) {
                    SSEUtils.makeNormalKernelInPlace(fYLo, changeInXNormalMean, changeInXNormalSd, nXbinsLo, nLeftNRightFlanksLo[0], nLeftNRightFlanksLo[1], dXbin); // normalizes inside already
                } else {
                    SSEUtils.makeNormalKernelInPlaceSSTJavaFftService(fYLo, changeInXNormalMean, changeInXNormalSd, nXbinsLo, nLeftNRightFlanksLo[0], nLeftNRightFlanksLo[1], dXbin);
                }

                // debugging
                // System.out.println("pre-FFT fYLo = " + Arrays.toString(fYLo));

                // FFTs normal kernel
                if (jtransforms && doFFT) fftForEandDLo.realForwardFull(fYLo); // with JTransforms
                else if (doFFT) {
                    jffts.fft(nXbinsLoSST, fYLo, fftFYLo); // result left in fftFYLo
                }
            }
            else {
                if (jtransforms) {
                    SSEUtils.makeNormalKernelInPlace(fYHi, changeInXNormalMean, changeInXNormalSd, nXbinsHi, nLeftNRightFlanksHi[0], nLeftNRightFlanksHi[1], dXbin / hiLoRatio); // normalizes inside already
                } else {
                    SSEUtils.makeNormalKernelInPlaceSSTJavaFftService(fYHi, changeInXNormalMean, changeInXNormalSd, nXbinsHi, nLeftNRightFlanksHi[0], nLeftNRightFlanksHi[1], dXbin / hiLoRatio); // normalizes inside already
                }

                // debugging
                // System.out.println("pre-FFT fYHi = " + Arrays.toString(fYHi));

                // FFTs normal kernel
                if (jtransforms && doFFT) fftForEandDHi.realForwardFull(fYHi); // with JTransforms
                else if (doFFT) {
                    jffts.fft(nXbinsHiSST, fYHi, fftFYHi); // result left in fftFYHi
                }
            }
        }
    }

    /*
     *
     */
    protected abstract void initializeEsDs(int nNodes, int nDimensionsFFT, int nXbinsLo, int nXbinsHi);

    protected void populateIndicesHiLo(int[] toArray, int ratioFromTo, int careAboutNTo, boolean jtransforms) {
        int multipl = 2; // SST's JavaFftService expects Es and Ds to be interdigitated
        if (jtransforms) multipl = 1;

        int nElementsCopied = 0;
        int i = ratioFromTo - 1; // offset because of indexing
        while (nElementsCopied < careAboutNTo) {
            toArray[nElementsCopied] = i * multipl;
            i += ratioFromTo;
            nElementsCopied++;
        }
    }

    /*
     *
     */
    protected abstract void populateTipsEsDs(int nDimensionsFFT, int nXbins, boolean ignoreRefresh, boolean jtransforms);

    /*
     * The D's at the root must be multiplied by a prior probability array.
     * This method populates it in place depending on what the user wants.
     */
    protected abstract void populatePriorProbAtRoot(double[] dsAtRoot, double dxAtRightRes, int nXbinAtRightRes, int nUsefulXbinAtRightRes, String rootPriorType, boolean jtransforms);

    /*
     *
     */
    public abstract void processBranch(Node aNode, boolean forceRecalcKernel, boolean jtransforms);

    /*
     *
     */
    protected abstract void integrateLength(int nodeIdx, double[][] esDsAtNode, double[][] scratchAtNode, double aLength, boolean dynamicallyAdjust, double maxDt, boolean lowRes, boolean forceRecalcKernel, boolean jtransforms);

    /*
     *
     */
    protected abstract double normalizeDs(int nodeIdx, boolean lowRes, boolean jtransforms);

    /*
     *
     */
    protected abstract void mergeChildrenNodes(Node aNode, boolean jtransforms);

    /*
     *
     */
    protected abstract void processInternalNode(Node aNode, boolean forceRecalcKernel, boolean jtransforms);

    /*
     *
     */
    protected abstract void startRecursionAtRootNode(Node rootNode, boolean forceRecalcKernel, boolean jtransforms);

    /*
     * Does integration in time and character space in place
     */
    protected abstract void doIntegrateInPlace(int nodeIdx, double dt, boolean lowRes);

    /*
     *
     */
    protected abstract void propagateTInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double dt, boolean lowRes, boolean jtransforms);

    /*
     *
     */
    // JTransforms version
    // protected abstract void propagateXInPlace(int nodeIdx, double[][] esDsAtNode, double[][] scratchAtNode, boolean lowRes);
    protected abstract void propagateXInPlace(double[][] esDsAtNode, double[][] fftBufferEsDsAtNode, double[][] scratchAtNode, boolean lowRes);

    /*
     * This method looks at the relevant objects in state,
     * computes the log-likelihood, and returns it
     */
    protected abstract double getLogPFromRelevantObjects(double[][] esDsAtRoot, double sumOfLogNormalizationFactors, double[] birthRates, double dXAtRightRes, boolean jtransforms);

    /*
     * Getters, setters and helper methods below
     */
    public int getnXbins(boolean lowRes) {
        if (lowRes) return nXbinsLo;
        else return nXbinsHi;
    }

    public int getNUsefulXbins(boolean lowRes) {
        if (lowRes) return nUsefulXbinsLo;
        else return nUsefulXbinsHi;
    }

    public int getNLeftFlanks(boolean lowRes) {
        if (lowRes) return nLeftNRightFlanksLo[0];
        else return nLeftNRightFlanksHi[0];
    }

    public int getNRightFlanks(boolean lowRes) {
        if (lowRes) return nLeftNRightFlanksLo[1];
        else return nLeftNRightFlanksHi[1];
    }

    public double getXMinLo() {
        return xMinLo;
    }

    public double getXMinHi() {
        return xMinHi;
    }

    public double[] getX(boolean lowRes) {
        if (lowRes) return xLo;
        else return xHi;
    }

    public double getdXbin() {
        return dXbin;
    }

    public int[] getHiLoIdxs4Transfer() {
        return hiLoIdxs4Transfer;
    }

    public double[] getPriorProbsAtRoot(String rootPriorType) {
        return priorProbsAtRoot;
    }

    // for debugging
    public void initializePriorProbAtRoot(int nxBinAtRightRes) {
        priorProbsAtRoot = new double[nxBinAtRightRes];
    }
}
