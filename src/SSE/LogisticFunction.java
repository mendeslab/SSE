package SSE;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Function;

/**
 * This class implements a type of link function that is used
 * by QuaSSE to convert an x value (continuous-trait value)
 * into a y value (macroevolutionary parameter, e.g., birth rate)
 *
 * @author Fabio K. Mendes
 */

@Description("Logistic link function for converting x into y, " +
        "where x is a continuous trait and y is a macroevolutionary" +
        "parameter.")
public class LogisticFunction extends CalculationNode implements LinkFn {

    final public Input<Function> curveYBaseValueInput = new Input<>("curveYBaseValue", "Curve y base value.", Input.Validate.REQUIRED);
    final public Input<Function> curveMaxYInput = new Input<>("curveMaxY", "Curve maximum y value.", Input.Validate.REQUIRED);
    final public Input<Function> sigmoidMidpointInput = new Input<>("sigmoidMidpoint", "Midpoint of sigmoid curve.", Input.Validate.REQUIRED);
    final public Input<Function> logisticGrowthRateInput = new Input<>("logisticGrowthRate", "Growth rate of logistic curve.", Input.Validate.REQUIRED);

    private double y0, y1, x0, r, curveMaxMinusBaseValue;
    private static final String LINKFUNCTION = "logistic";

    @Override
    public void initAndValidate() {
        y0 = curveYBaseValueInput.get().getArrayValue();
        y1 = curveMaxYInput.get().getArrayValue();
        r = logisticGrowthRateInput.get().getArrayValue();
        x0 = sigmoidMidpointInput.get().getArrayValue();
    }

    @Override
    public boolean refreshParams() {

        boolean refreshedSomething = false;

        // Compare values: rejection restores parameters, but not these cached scalars.
        double currentY0 = curveYBaseValueInput.get().getArrayValue();
        if (y0 != currentY0) {
            y0 = currentY0;
            refreshedSomething = true;
        }

        double currentY1 = curveMaxYInput.get().getArrayValue();
        if (y1 != currentY1) {
            y1 = currentY1;
            refreshedSomething = true;
        }

        double currentR = logisticGrowthRateInput.get().getArrayValue();
        if (r != currentR) {
            r = currentR;
            refreshedSomething = true;
        }

        double currentX0 = sigmoidMidpointInput.get().getArrayValue();
        if (x0 != currentX0) {
            x0 = currentX0;
            refreshedSomething = true;
        }

        return refreshedSomething;
    }

    @Override
    public double[] getY(double[] x, double[] y, boolean ignoreRefresh) {
        boolean refreshedSomething = false;
        if (!ignoreRefresh) refreshedSomething = refreshParams(); // if something changed in deterministic function parameters, we need to repopulate macroevol arrays

        /*
         * if either we don't care about refreshing, or we do and something was refreshed,
         * we repopulate y
         */
        if (ignoreRefresh || refreshedSomething) {
            if (x.length != y.length) throw new RuntimeException("Sizes of x (qu trait) and y (macroevol param) differ. Exiting...");

            for (int i = 0; i < x.length; i++) {
                curveMaxMinusBaseValue = y1 - y0;
                y[i] = y0 + curveMaxMinusBaseValue / (1.0 + Math.exp(r * (x0 - x[i])));
            }
        }

        return y;
    }

    @Override
    public String getLinkFnName() {
        return LINKFUNCTION;
    }

}
