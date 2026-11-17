package SSE;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;

@Description("Gaussian rate peak or trough between nonnegative center and tail rates.")
public class GaussianLinkFn extends CalculationNode implements LinkFn {
    public final Input<Function> baselineInput = new Input<>("baseline", "Rate in the tails.", Input.Validate.REQUIRED);
    public final Input<Function> centerRateInput = new Input<>("center_rate", "Rate at the center.", Input.Validate.REQUIRED);
    public final Input<Function> centerInput = new Input<>("center", "Center trait value.", Input.Validate.REQUIRED);
    public final Input<Function> widthInput = new Input<>("width", "Positive Gaussian width (not variance).",
            Input.Validate.REQUIRED);

    private double baseline, centerRate, center, width;

    @Override
    public void initAndValidate() {
        refreshParams();
    }

    // Validate all inputs before updating the cache; value comparisons also detect restored parameters.
    @Override
    public boolean refreshParams() {
        Function a = baselineInput.get(), b = centerRateInput.get(), c = centerInput.get(), w = widthInput.get();
        if (a.getDimension() != 1 || b.getDimension() != 1 || c.getDimension() != 1 || w.getDimension() != 1)
            throw new IllegalArgumentException("Gaussian rate parameters must be scalars.");
        double nextBaseline = a.getArrayValue(), nextCenterRate = b.getArrayValue();
        double nextCenter = c.getArrayValue(), nextWidth = w.getArrayValue();
        if (!Double.isFinite(nextBaseline) || nextBaseline < 0
                || !Double.isFinite(nextCenterRate) || nextCenterRate < 0
                || !Double.isFinite(nextCenter) || !Double.isFinite(nextWidth) || nextWidth <= 0)
            throw new IllegalArgumentException(
                    "Gaussian rates must be finite and nonnegative; center finite; width positive and finite.");
        boolean changed = baseline != nextBaseline || centerRate != nextCenterRate
                || center != nextCenter || width != nextWidth;
        baseline = nextBaseline;
        centerRate = nextCenterRate;
        center = nextCenter;
        width = nextWidth;
        return changed;
    }

    // Evaluate each grid point as a convex combination of tail and center rates. With t = -z²/2,
    // -expm1(t) retains the small tail weight near the center; exp(t) retains the center weight
    // far away. Equal levels give a constant, and the center itself has exactly the center rate.
    @Override
    public double[] getY(double[] x, double[] y, boolean ignoreRefresh) {
        boolean changed = !ignoreRefresh && refreshParams();
        if (ignoreRefresh || changed) {
            if (x.length != y.length) throw new IllegalArgumentException("Trait and rate arrays must have equal lengths.");
            for (int i = 0; i < x.length; i++) {
                if (baseline == centerRate) {
                    y[i] = baseline;
                    continue;
                }
                double delta = x[i] - center;
                // For finite traits, rescale an overflowed subtraction before taking the difference.
                // Opposite signs caused that overflow, so this cannot introduce cancellation.
                double z = Double.isInfinite(delta) ? x[i] / width - center / width : delta / width;
                double t = -0.5 * z * z;
                y[i] = baseline * -Math.expm1(t) + centerRate * Math.exp(t);
            }
        }
        return y;
    }

    @Override
    public String getLinkFnName() {
        return "gaussian";
    }
}
