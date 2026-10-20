package SSE;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.distribution.ParametricDistribution;
import org.apache.commons.math.distribution.CauchyDistributionImpl;
import org.apache.commons.math.distribution.ContinuousDistribution;

@Description("Cauchy distribution with a location (median) and positive scale. Mean and variance are undefined.")
public class Cauchy extends ParametricDistribution {
    public final Input<Function> locationInput = new Input<>("location", "Location (median)", Input.Validate.REQUIRED);
    public final Input<Function> scaleInput = new Input<>("scale", "Positive scale, not standard deviation", Input.Validate.REQUIRED);
    private final CauchyDistributionImpl distribution = new CauchyDistributionImpl();

    @Override
    public void initAndValidate() {
        getDistribution();
    }

    // Refresh from BEAST parameters on every access, including after proposal restoration.
    @Override
    @SuppressWarnings("deprecation")
    public ContinuousDistribution getDistribution() {
        Function location = locationInput.get();
        Function scale = scaleInput.get();
        if (location == null || scale == null || location.getDimension() != 1 || scale.getDimension() != 1)
            throw new IllegalArgumentException("Cauchy location and scale must be scalar parameters.");
        double m = location.getArrayValue(), s = scale.getArrayValue();
        if (!Double.isFinite(m) || !Double.isFinite(s) || s <= 0)
            throw new IllegalArgumentException("Cauchy requires finite location and finite positive scale.");
        distribution.setMedian(m);
        distribution.setScale(s);
        return distribution;
    }

    @Override
    public double getMeanWithoutOffset() {
        return Double.NaN;
    }
}
