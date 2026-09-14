package SSE;

import SSE.fft.ComplexFFT;
import SSE.fft.SstFFT;
import SSE.fft.FftwFFT;

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

import java.util.Arrays;

@Description("Specifies a quantitative trait(s) state-dependent speciation and" +
        "extinction birth-death process.")
public abstract class QuaSSEProcess extends Distribution implements AutoCloseable {

    final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> dynamicDtInput = new Input<>("dynDt", "If interval over which to carry out integration should be dynamically adjusted to maximize accuracy.", Input.Validate.REQUIRED);
    final public Input<String> priorProbAtRootTypeInput = new Input<>("priorProbAtRootType", "Type of root prior probabilities for D's.", Input.Validate.REQUIRED);
    final public Input<RealParameter> priorProbsAtRootInput = new Input<>("givenPriorProbsAtRoot", "Root prior probabilities for D's at high resolution.", Input.Validate.XOR, priorProbAtRootTypeInput);
    final public Input<String> fftBackendInput = new Input<>("fftBackend",
            "Fourier transform implementation: sst (Java, default) or fftw (requires native library).", "sst");

    public final Input<QuaSSEGrid> gridInput = new Input<>("grid", "Fine/coarse grid; omit for observation-derived defaults.");
    protected QuaSSEGrid grid;
    private int gridRevision = -1;

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
    protected double dXbin, xMinLo, xMinHi;
    protected int nXbinsLo, nUsefulXbinsLo, nXbinsHi, nUsefulXbinsHi, hiLoRatio;
    protected int[] nLeftNRightFlanksHi, nLeftNRightFlanksLo;
    protected double[] xLo, xHi; // x rulers
    protected int[] hiLoIdxs4Transfer;

    // quantitative trait evolution
    protected double changeInXNormalMean; // (=diversitree's drift)
    protected double changeInXNormalSd; // (=diversitree's diffusion)
    protected double[] fYLo, fYHi, fftFYLo, fftFYHi;
    protected DoubleFFT_1D fftForEandDLo, fftForEandDHi;
    // One kernel per resolution (0=coarse, 1=fine). These arrays describe the contents of fY/fftFY,
    // not the accepted MCMC state; comparing values also handles clean restored parameters.
    private final boolean[] kernelValid = new boolean[2], kernelHasFFT = new boolean[2];
    private final boolean[] kernelJTransforms = new boolean[2];
    private final double[] kernelDt = new double[2], kernelDrift = new double[2];
    private final double[] kernelDiffusion = new double[2], kernelDx = new double[2];
    private final int[] kernelLeft = new int[2], kernelRight = new int[2];
    // One per resolution, shared by kernel and E/D transforms. Padding/parameter refresh reuses
    // them; reinitialization recreates them and invalidates kernels. Native scratch needs no snapshot.
    protected ComplexFFT fftLo, fftHi;

    @Override
    public void initAndValidate() {

        close();
        // Select once per initialization, not dynamically during likelihood evaluation.
        String backend = fftBackendInput.get();
        if (!"sst".equals(backend) && !"fftw".equals(backend))
            throw new IllegalArgumentException("fftBackend must be sst or fftw.");
        tree = treeInput.get();
        rootPriorType = priorProbAtRootTypeInput.get();

        if (priorProbsAtRootInput.get() != null) {
            throw new IllegalArgumentException("givenPriorProbsAtRoot is not supported: a supplied array "
                    + "has no defined mapping when the QuaSSE grid changes. Use Flat or Observed.");
        }
        if (!FLAT.equals(rootPriorType) && !OBS.equals(rootPriorType))
            throw new IllegalArgumentException("priorProbAtRootType must be Flat or Observed.");
        // quTraits = quTraitsInput.get();

        dynamicallyAdjustDt = dynamicDtInput.get().getValue();
        grid = gridInput.get();
        gridRevision = -1;
        if (!refreshGridGeometry()) throw new IllegalArgumentException("Invalid QuaSSE grid.");

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
        try {
            fftLo = "fftw".equals(backend) ? new FftwFFT(nXbinsLo) : new SstFFT(nXbinsLo);
            fftHi = "fftw".equals(backend) ? new FftwFFT(nXbinsHi) : new SstFFT(nXbinsHi);
        } catch (RuntimeException | Error failure) {
            close();
            throw failure;
        }
        Arrays.fill(kernelValid, false);

        // populatefY(dtMax, true, true, true, false); // force populate fY, and do FFT (coarse)
        // populatefY(dtMax, true, true, false, false); // force populate fY, and do FFT (fine)
    }

    // Release both resolutions on explicit closure or before reinitializing the distribution.
    @Override
    public void close() {
        if (fftLo != null) fftLo.close();
        if (fftHi != null) fftHi.close();
    }

    // Copy geometry only after a layout change; integration never mutates the grid node's arrays.
    // Kernels read current drift/diffusion directly from the grid, even without a layout change.
    protected boolean refreshGridGeometry() {
        if (!grid.update()) return false;
        if (gridRevision == grid.getRevision()) return true;
        gridRevision = grid.getRevision();
        dtMax = grid.getDtMax();
        tc = grid.getTc();
        dXbin = grid.getDx();
        hiLoRatio = grid.getRatio();
        nXbinsLo = grid.getBins(true);
        nXbinsHi = grid.getBins(false);
        nLeftNRightFlanksLo = grid.copyPadding(true);
        nLeftNRightFlanksHi = grid.copyPadding(false);
        xLo = grid.copyX(true);
        xHi = grid.copyX(false);
        xMinLo = xLo[0];
        xMinHi = xHi[0];
        nUsefulXbinsLo = xLo.length;
        nUsefulXbinsHi = xHi.length;
        hiLoIdxs4Transfer = grid.copyTransfer();
        return true;
    }

    protected abstract void populateMacroevolParams();

    /*
     *
     */
    // Reuse only a kernel whose parameters, resolution and representation match this request.
    protected void populatefY(double aDt, boolean forceRecalcKernel,
                              boolean doFFT, boolean lowRes, boolean jtransforms) {
        int resolution = lowRes ? 0 : 1;
        int size = lowRes ? nXbinsLo : nXbinsHi;
        int[] padding = lowRes ? nLeftNRightFlanksLo : nLeftNRightFlanksHi;
        double dx = lowRes ? dXbin : dXbin / hiLoRatio;
        double currentDrift = grid.getDrift();
        double currentDiffusion = grid.getDiffusion();
        if (!Double.isFinite(aDt) || aDt <= 0 || !Double.isFinite(currentDrift)
                || !Double.isFinite(currentDiffusion) || currentDiffusion <= 0)
            throw new IllegalArgumentException("A QuaSSE kernel requires finite drift and positive finite dt/diffusion.");
        // The actual per-resolution key decides validity, including changes to the timestep.
        // Require the same raw/FFT request: JTransforms overwrites its raw array with the spectrum.
        if (!forceRecalcKernel && kernelValid[resolution] && kernelDt[resolution] == aDt
                && kernelDrift[resolution] == currentDrift && kernelDiffusion[resolution] == currentDiffusion
                && kernelDx[resolution] == dx && kernelLeft[resolution] == padding[0]
                && kernelRight[resolution] == padding[1] && kernelJTransforms[resolution] == jtransforms
                && kernelHasFFT[resolution] == doFFT)
            return;

        kernelValid[resolution] = false;
        double[] raw = lowRes ? fYLo : fYHi;
        double[] spectrum = lowRes ? fftFYLo : fftFYHi;
        Arrays.fill(raw, 0.0);
        Arrays.fill(spectrum, 0.0);
        changeInXNormalMean = -currentDrift * aDt;
        changeInXNormalSd = Math.sqrt(currentDiffusion * aDt);
        if (jtransforms) {
            SSEUtils.makeNormalKernelInPlace(raw, changeInXNormalMean, changeInXNormalSd,
                    size, padding[0], padding[1], dx);
            if (doFFT) (lowRes ? fftForEandDLo : fftForEandDHi).realForwardFull(raw);
        } else {
            SSEUtils.makeNormalKernelInPlaceSSTJavaFftService(raw, changeInXNormalMean, changeInXNormalSd,
                    size, padding[0], padding[1], dx);
            if (doFFT) (lowRes ? fftLo : fftHi).forward(raw, spectrum);
        }
        kernelDt[resolution] = aDt;
        kernelDrift[resolution] = currentDrift;
        kernelDiffusion[resolution] = currentDiffusion;
        kernelDx[resolution] = dx;
        kernelLeft[resolution] = padding[0];
        kernelRight[resolution] = padding[1];
        kernelJTransforms[resolution] = jtransforms;
        kernelHasFFT[resolution] = doFFT;
        kernelValid[resolution] = true;
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
