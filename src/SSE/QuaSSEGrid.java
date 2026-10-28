package SSE;

import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;
import java.util.Arrays;

/**
 * Fine/coarse grid geometry with fixed bin counts and defaults derived from current parameters.
 * BEAST notifications and restore invalidate derived data. Callers must use BEAST's lifecycle.
 */
public class QuaSSEGrid extends CalculationNode {
    public final Input<Tree> treeInput = new Input<>("tree", "Tree supplying the root age for root conditioning.");
    public final Input<RealParameter> originInput = new Input<>("origin",
            "Origin age for origin conditioning.", Input.Validate.XOR, treeInput);
    public final Input<RealParameter> traitsInput = new Input<>("traits", "Observed traits for spatial defaults.", Input.Validate.REQUIRED);
    public final Input<RealParameter> driftInput = new Input<>("drift", "Trait drift.", Input.Validate.REQUIRED);
    public final Input<RealParameter> diffusionInput = new Input<>("diffusion", "Trait diffusion variance rate.", Input.Validate.REQUIRED);
    public final Input<Integer> nXInput = new Input<>("nX", "Coarse bins; defaults to 1024.");
    public final Input<Double> dXInput = new Input<>("dX", "Coarse spacing; defaults cover observations, diffusion and padding.");
    public final Input<Double> xMidInput = new Input<>("xMid", "Midpoint; defaults to the center of the observation/diffusion interval.");
    public final Input<Double> rangeMultiplierInput = new Input<>("rangeMultiplier", "Observed range multiplier.", 5.0);
    public final Input<Integer> hiLoRatioInput = new Input<>("hiLoRatio", "Fine/coarse resolution ratio.", 4);
    public final Input<Double> flankWidthScalerInput = new Input<>("flankWidthScaler", "Kernel support width in SDs.", 5.0);
    public final Input<Double> dtMaxInput = new Input<>("dtMax", "Maximum step; defaults to current origin or root age / 250.");
    public final Input<Double> tcInput = new Input<>("tc", "Resolution switch; defaults to current origin or root age / 10.");

    private double dtMax, tc, dXbin, xMid, flankWidthScaler, drift, diffusion;
    private int nXbinsLo, nXbinsHi, hiLoRatio, revision;
    private int[] paddingLo, paddingHi, transfer;
    private double[] xLo, xHi;
    private boolean needsUpdate = true, valid;

    // Resolve dimensions once. All omitted physical controls are derived by update(), including
    // before initial validation, so a singleton can use its origin rather than a zero root age.
    @Override
    public void initAndValidate() {
        if ((treeInput.get() == null) == (originInput.get() == null))
            throw new IllegalArgumentException("QuaSSEGrid requires exactly one of tree and origin.");
        nXbinsLo = nXInput.get() == null ? 1024 : nXInput.get();
        // Preserve spacing-only configuration: choose dimensions once from the observed range.
        if (nXInput.get() == null && dXInput.get() != null) {
            double[] values = traitsInput.get().getDoubleValues();
            double width = (Arrays.stream(values).max().orElse(Double.NaN)
                    - Arrays.stream(values).min().orElse(Double.NaN)) * rangeMultiplierInput.get();
            double count = Math.pow(2, Math.ceil(Math.log(width / dXInput.get()) / Math.log(2)));
            if (!Double.isFinite(count) || count < 1 || count > Integer.MAX_VALUE)
                throw new IllegalArgumentException("Cannot derive nX from trait range and dX; supply nX explicitly.");
            nXbinsLo = (int) count;
        }
        hiLoRatio = hiLoRatioInput.get();
        flankWidthScaler = flankWidthScalerInput.get();
        if (nXbinsLo <= 0 || (nXbinsLo & (nXbinsLo - 1)) != 0)
            throw new IllegalArgumentException("Number of quantitative character bins must be a power of 2. It was " + nXbinsLo);
        if (hiLoRatio < 1 || 2L * nXbinsLo * hiLoRatio > Integer.MAX_VALUE)
            throw new IllegalArgumentException("QuaSSE hiLoRatio must be positive and FFT array sizes must fit in an int.");
        if (!Double.isFinite(flankWidthScaler) || flankWidthScaler <= 0
                || !Double.isFinite(rangeMultiplierInput.get()) || rangeMultiplierInput.get() <= 0)
            throw new IllegalArgumentException("QuaSSE rangeMultiplier and flankWidthScaler must be positive and finite.");
        nXbinsHi = nXbinsLo * hiLoRatio;
        paddingLo = null;
        needsUpdate = true;
        if (!update()) throw new IllegalArgumentException("Invalid QuaSSE grid controls or defaults; "
                + "check origin/root age, drift, diffusion, and kernel support; supply explicit dX/dtMax if needed.");
    }

    // BEAST calls this only along affected dependency paths; geometry is rebuilt on demand.
    @Override
    protected boolean requiresRecalculation() {
        needsUpdate = true;
        return true;
    }

    // The arrays may describe a rejected proposal even though BEAST resets its own dirty flag.
    @Override
    protected void restore() {
        needsUpdate = true;
        super.restore();
    }

    // Construct candidate geometry entirely in locals: invalid proposals must not alter the last
    // usable grid. A revision includes time controls, even when no spatial coordinates change.
    public boolean update() {
        if (!needsUpdate) return valid;
        needsUpdate = false;
        valid = false;
        // Read current parameters and validate the physical time controls.
        double candidateDrift = driftInput.get().getValue(), candidateDiffusion = diffusionInput.get().getValue();
        RealParameter origin = originInput.get();
        if (origin != null && origin.getDimension() != 1) return false;
        double age = origin == null ? treeInput.get().getRoot().getHeight() : origin.getValue();
        double candidateDt = dtMaxInput.get() == null ? age / 250 : dtMaxInput.get();
        double candidateTc = tcInput.get() == null ? age / 10 : tcInput.get();
        if (!Double.isFinite(age) || age < 0 || !Double.isFinite(candidateDrift)
                || !Double.isFinite(candidateDiffusion) || candidateDiffusion <= 0
                || !Double.isFinite(candidateDt) || candidateDt <= 0
                || !Double.isFinite(candidateTc) || candidateTc < 0) return false;

        // Determine the kernel's reach and the trait interval the active grid must cover.
        double[] support = calculateKernelSupportDistances(candidateDrift, candidateDiffusion, candidateDt);
        double leftDistance = support[0], rightDistance = support[1];
        double[] bounds = calculateTraitBounds(candidateDrift, candidateDiffusion, age);
        double lower = bounds[0], upper = bounds[1];

        // Apply explicit midpoint and spacing overrides where supplied.
        double candidateMid = xMidInput.get() == null ? lower / 2 + upper / 2 : xMidInput.get();
        boolean automaticSpacing = dXInput.get() == null;
        double candidateDx = calculateSpacing(upper - lower, leftDistance, rightDistance);
        if (!Double.isFinite(candidateDx) || candidateDx <= 0 || !Double.isFinite(candidateMid)
                || (automaticSpacing && (!(upper > lower) || nXbinsLo <= 6))) return false;

        // Round support into padding bins, then place and validate both resolutions.
        int[] padding = calculatePadding(leftDistance, rightDistance, candidateDx);
        if (padding == null) return false;
        int[] high = {padding[0] * hiLoRatio, padding[1] * hiLoRatio};
        int usefulLo = nXbinsLo - padding[0] - padding[1] - 1;
        int usefulHi = nXbinsHi - high[0] - high[1] - 1;
        // Convolution preserves pR on the left and pL on the right. Center the active span;
        // explicit spacing retains the historical placement, including its bin rounding.
        double minLo = candidateMid - candidateDx * (automaticSpacing
                ? (usefulLo - 1.0 + padding[1] - padding[0]) / 2 : Math.ceil((usefulLo - 1.0) / 2));
        double minHi = minLo - candidateDx * (1.0 - 1.0 / hiLoRatio);
        if (candidateDx / hiLoRatio == 0 || !Double.isFinite(minLo) || !Double.isFinite(minHi)
                || !Double.isFinite(minLo + (usefulLo - 1.0) * candidateDx)
                || !Double.isFinite(minLo + (usefulLo - 1.0 / hiLoRatio) * candidateDx)) return false;

        // Reuse unchanged geometry; otherwise construct arrays before publishing the new grid.
        valid = true;
        drift = candidateDrift;
        diffusion = candidateDiffusion;
        if (Arrays.equals(padding, paddingLo) && candidateDx == dXbin && candidateMid == xMid
                && candidateDt == dtMax && candidateTc == tc) return true;
        double[] lowX = new double[usefulLo], highX = new double[usefulHi];
        lowX[0] = minLo;
        highX[0] = minHi;
        // Preserve accumulated-addition rounding from the existing discretization.
        for (int i = 1; i < usefulLo; ++i) lowX[i] = lowX[i - 1] + candidateDx;
        for (int i = 1; i < usefulHi; ++i) highX[i] = highX[i - 1] + candidateDx / hiLoRatio;
        int[] indices = new int[usefulLo];
        for (int i = 0; i < usefulLo; ++i) indices[i] = 2 * (hiLoRatio - 1 + i * hiLoRatio);
        dXbin = candidateDx;
        xMid = candidateMid;
        dtMax = candidateDt;
        tc = candidateTc;
        paddingLo = padding;
        paddingHi = high;
        xLo = lowX;
        xHi = highX;
        transfer = indices;
        ++revision;
        return true;
    }

    // Return {left, right} physical support distances covering all steps up to candidateDt.
    private double[] calculateKernelSupportDistances(double candidateDrift, double candidateDiffusion,
            double candidateDt) {
        // Backward support extents are ±v*t+a*sqrt(t). The growing side peaks at candidateDt;
        // the opposing side peaks at min(candidateDt,(a/(2*abs(v)))²). These physical distances
        // bound every integration step before they are rounded outward into padding bins.
        double a = flankWidthScaler * Math.sqrt(candidateDiffusion), speed = Math.abs(candidateDrift);
        double peakRoot = speed == 0 ? Math.sqrt(candidateDt) : Math.min(Math.sqrt(candidateDt), a / (2 * speed));
        double growing = speed * candidateDt + a * Math.sqrt(candidateDt);
        double opposing = peakRoot * (a - speed * peakRoot);
        double leftDistance = candidateDrift >= 0 ? growing : opposing;
        double rightDistance = candidateDrift >= 0 ? opposing : growing;
        return new double[] {leftDistance, rightDistance};
    }

    // Return {lower, upper}: the union of the observation interval and the backward diffusion
    // interval. Six sigma is the total diffusion width, centered at the midpoint minus drift*age.
    private double[] calculateTraitBounds(double candidateDrift, double candidateDiffusion, double age) {
        double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
        for (double value : traitsInput.get().getDoubleValues()) {
            min = Math.min(min, value);
            max = Math.max(max, value);
        }
        double middle = min / 2 + max / 2;
        double halfObserved = (max - min) * rangeMultiplierInput.get() / 2;
        double halfDiffusion = 3 * Math.sqrt(candidateDiffusion * age);
        double lower = Math.min(middle - halfObserved, middle - candidateDrift * age - halfDiffusion);
        double upper = Math.max(middle + halfObserved, middle - candidateDrift * age + halfDiffusion);
        return new double[] {lower, upper};
    }

    // Use explicit spacing or allow enough width for the requested interval and kernel padding.
    private double calculateSpacing(double requestedWidth, double leftDistance, double rightDistance) {
        // N bins reserve P+1 for FFT padding and P for fixed boundary strips, P=pL+pR.
        // The active center span is (N-2P-2)*dx. Since ceil(aL/dx)+ceil(aR/dx)
        // is at most (aL+aR)/dx+2, dx=(W+2(aL+aR))/(N-6) leaves span >= W.
        return dXInput.get() == null
                ? (requestedWidth + 2 * (leftDistance + rightDistance)) / (nXbinsLo - 6) : dXInput.get();
    }

    // Round physical support distances outward using the candidate spacing, without changing
    // the live grid. Null denotes nonfinite padding or overlapping preserved boundary strips.
    private int[] calculatePadding(double leftDistance, double rightDistance, double candidateDx) {
        double left = Math.ceil(leftDistance / candidateDx), right = Math.ceil(rightDistance / candidateDx);
        // The preserved strips cannot overlap. This also bounds integer casts and the fine map.
        if (!Double.isFinite(left) || !Double.isFinite(right) || left < 1 || right < 1
                || 2 * (left + right) + 1 >= nXbinsLo) return null;
        return new int[] {(int) left, (int) right};
    }

    // Call update() successfully before reading derived values. Copies protect owned geometry.
    public double[] copyX(boolean low) { return (low ? xLo : xHi).clone(); }
    public int[] copyPadding(boolean low) { return (low ? paddingLo : paddingHi).clone(); }
    public int[] copyTransfer() { return transfer.clone(); }
    public int getRevision() { return revision; }
    public int getBins(boolean low) { return low ? nXbinsLo : nXbinsHi; }
    public int getRatio() { return hiLoRatio; }
    public double getDx() { return dXbin; }
    public double getMidpoint() { return xMid; }
    public double getDtMax() { return dtMax; }
    public double getTc() { return tc; }
    public double getDrift() { return drift; }
    public double getDiffusion() { return diffusion; }
}
