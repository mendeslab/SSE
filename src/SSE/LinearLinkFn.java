package SSE;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.CalculationNode;

@Description("Nonnegative linear rate with optional constant tails at fixed trait bounds.")
public class LinearLinkFn extends CalculationNode implements LinkFn {
    public final Input<Function> interceptInput = new Input<>("intercept", "Linear intercept.", Input.Validate.REQUIRED);
    public final Input<Function> slopeInput = new Input<>("slope", "Linear slope.", Input.Validate.REQUIRED);
    public final Input<Double> lowerTraitInput = new Input<>("lower_trait", "Optional fixed lower trait bound.");
    public final Input<Double> upperTraitInput = new Input<>("upper_trait", "Optional fixed upper trait bound.");

    private double intercept, slope;
    private Double lowerTrait, upperTrait;

    // Bounds are fixed settings; only the scalar coefficients participate in parameter refresh.
    @Override
    public void initAndValidate() {
        lowerTrait = lowerTraitInput.get();
        upperTrait = upperTraitInput.get();
        if ((lowerTrait != null && !Double.isFinite(lowerTrait))
                || (upperTrait != null && !Double.isFinite(upperTrait))
                || (lowerTrait != null && upperTrait != null && lowerTrait > upperTrait)) {
            throw new IllegalArgumentException("Trait bounds must be finite and lower_trait <= upper_trait.");
        }
        refreshParams();
    }

    // Compare values rather than dirty flags: restored BEAST parameters can already be clean.
    @Override
    public boolean refreshParams() {
        Function a = interceptInput.get(), b = slopeInput.get();
        double currentIntercept = a.getArrayValue(), currentSlope = b.getArrayValue();
        if (a.getDimension() != 1 || b.getDimension() != 1
                || !Double.isFinite(currentIntercept) || !Double.isFinite(currentSlope)) {
            throw new IllegalArgumentException("Linear intercept and slope must be finite scalars.");
        }
        boolean changed = intercept != currentIntercept || slope != currentSlope;
        intercept = currentIntercept;
        slope = currentSlope;
        return changed;
    }

    // Clamp the trait, evaluate the line, then replace negative rates with zero (diversitree order).
    @Override
    public double[] getY(double[] x, double[] y, boolean ignoreRefresh) {
        boolean changed = !ignoreRefresh && refreshParams();
        if (ignoreRefresh || changed) {
            if (x.length != y.length) throw new IllegalArgumentException("Trait and rate arrays must have equal lengths.");
            for (int i = 0; i < x.length; i++) {
                double trait = x[i];
                if (lowerTrait != null) trait = Math.max(lowerTrait, trait);
                if (upperTrait != null) trait = Math.min(upperTrait, trait);
                y[i] = Math.max(0.0, intercept + slope * trait);
            }
        }
        return y;
    }

    @Override
    public String getLinkFnName() {
        return "linear";
    }
}
