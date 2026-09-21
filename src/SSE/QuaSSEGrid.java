package SSE;

import beast.base.core.Input;
import beast.base.inference.CalculationNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;
import java.util.Arrays;

/**
 * Fine/coarse grid geometry. Defaults are resolved once; only drift/diffusion change its padding.
 * BEAST notifications and restore invalidate derived data. Callers must use BEAST's lifecycle.
 */
public class QuaSSEGrid extends CalculationNode {
    public final Input<Tree> treeInput = new Input<>("tree", "Initialization tree for time defaults.", Input.Validate.REQUIRED);
    public final Input<RealParameter> traitsInput = new Input<>("traits", "Observed traits for spatial defaults.", Input.Validate.REQUIRED);
    public final Input<RealParameter> driftInput = new Input<>("drift", "Trait drift.", Input.Validate.REQUIRED);
    public final Input<RealParameter> diffusionInput = new Input<>("diffusion", "Trait diffusion variance rate.", Input.Validate.REQUIRED);
    public final Input<Integer> nXInput = new Input<>("nX", "Coarse bins; defaults to 1024.");
    public final Input<Double> dXInput = new Input<>("dX", "Coarse spacing; derived from the trait range.");
    public final Input<Double> xMidInput = new Input<>("xMid", "Midpoint; defaults to the observed range midpoint.");
    public final Input<Double> rangeMultiplierInput = new Input<>("rangeMultiplier", "Observed range multiplier.", 5.0);
    public final Input<Integer> hiLoRatioInput = new Input<>("hiLoRatio", "Fine/coarse resolution ratio.", 4);
    public final Input<Double> flankWidthScalerInput = new Input<>("flankWidthScaler", "Kernel support width in SDs.", 5.0);
    public final Input<Double> dtMaxInput = new Input<>("dtMax", "Maximum step; defaults to initial tree height / 250.");
    public final Input<Double> tcInput = new Input<>("tc", "Resolution switch; defaults to initial tree height / 10.");

    private double dtMax, tc, dXbin, xMid, flankWidthScaler, drift, diffusion;
    private int nXbinsLo, nXbinsHi, hiLoRatio, revision;
    private int[] paddingLo, paddingHi, transfer;
    private double[] xLo, xHi;
    private boolean needsUpdate = true, valid;

    // Resolve configuration without changing it during proposals; explicit reinitialization resolves anew.
    @Override
    public void initAndValidate() {
        double min = Double.POSITIVE_INFINITY, max = Double.NEGATIVE_INFINITY;
        for (double value : traitsInput.get().getDoubleValues()) {
            min = Math.min(min, value);
            max = Math.max(max, value);
        }
        double width = (max - min) * rangeMultiplierInput.get();
        double height = treeInput.get().getRoot().getHeight();
        nXbinsLo = nXInput.get() == null ? 1024 : nXInput.get();
        // A spacing-only override retains at least the requested range by rounding bins up to a power of two.
        if (nXInput.get() == null && dXInput.get() != null) {
            double count = Math.pow(2, Math.ceil(Math.log(width / dXInput.get()) / Math.log(2)));
            if (!Double.isFinite(count) || count < 1 || count > Integer.MAX_VALUE)
                throw new IllegalArgumentException("Cannot derive nX from trait range and dX; supply nX explicitly.");
            nXbinsLo = (int) count;
        }
        dXbin = dXInput.get() == null ? width / nXbinsLo : dXInput.get();
        xMid = xMidInput.get() == null ? (min + max) / 2 : xMidInput.get();
        dtMax = dtMaxInput.get() == null ? height / 250 : dtMaxInput.get();
        tc = tcInput.get() == null ? height / 10 : tcInput.get();
        hiLoRatio = hiLoRatioInput.get();
        flankWidthScaler = flankWidthScalerInput.get();
        if (nXbinsLo <= 0 || (nXbinsLo & (nXbinsLo - 1)) != 0)
            throw new IllegalArgumentException("Number of quantitative character bins must be a power of 2. It was " + nXbinsLo);
        if (hiLoRatio < 1 || 2L * nXbinsLo * hiLoRatio > Integer.MAX_VALUE)
            throw new IllegalArgumentException("QuaSSE hiLoRatio must be positive and FFT array sizes must fit in an int.");
        if (!Double.isFinite(dXbin) || dXbin <= 0 || !Double.isFinite(xMid)
                || !Double.isFinite(dtMax) || dtMax <= 0 || !Double.isFinite(tc) || tc < 0
                || !Double.isFinite(flankWidthScaler) || flankWidthScaler <= 0
                || !Double.isFinite(rangeMultiplierInput.get()) || rangeMultiplierInput.get() <= 0)
            throw new IllegalArgumentException("Invalid QuaSSE grid controls or defaults; supply positive dX/dtMax, "
                    + "finite xMid, nonnegative tc, and positive finite rangeMultiplier/flankWidthScaler.");
        nXbinsHi = nXbinsLo * hiLoRatio;
        paddingLo = null;
        needsUpdate = true;
        if (!update()) throw new IllegalArgumentException("Invalid QuaSSE drift/diffusion or insufficient grid for kernel support.");
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

    // Build all geometry before publishing; invalid proposals leave the last usable arrays intact.
    public boolean update() {
        if (!needsUpdate) return valid;
        needsUpdate = false;
        drift = driftInput.get().getValue();
        diffusion = diffusionInput.get().getValue();
        int[] padding = calculatePadding(drift, diffusion);
        valid = padding != null;
        if (!valid || Arrays.equals(padding, paddingLo)) return valid;
        int[] high = {padding[0] * hiLoRatio, padding[1] * hiLoRatio};
        int usefulLo = nXbinsLo - padding[0] - padding[1] - 1;
        int usefulHi = nXbinsHi - high[0] - high[1] - 1;
        double minLo = xMid - dXbin * Math.ceil((usefulLo - 1.0) / 2.0);
        double minHi = minLo - dXbin * (1.0 - 1.0 / hiLoRatio);
        double[] lowX = new double[usefulLo], highX = new double[usefulHi];
        lowX[0] = minLo;
        highX[0] = minHi;
        // Preserve accumulated-addition rounding from the existing discretization.
        for (int i = 1; i < usefulLo; ++i) lowX[i] = lowX[i - 1] + dXbin;
        for (int i = 1; i < usefulHi; ++i) highX[i] = highX[i - 1] + dXbin / hiLoRatio;
        int[] indices = new int[usefulLo];
        for (int i = 0; i < usefulLo; ++i) indices[i] = 2 * (hiLoRatio - 1 + i * hiLoRatio);
        paddingLo = padding;
        paddingHi = high;
        xLo = lowX;
        xHi = highX;
        transfer = indices;
        ++revision;
        return true;
    }

    // Calculate candidate support without mutating the live grid. Null denotes an invalid proposal.
    private int[] calculatePadding(double candidateDrift, double candidateDiffusion) {
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
