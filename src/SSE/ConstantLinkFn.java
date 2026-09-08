package SSE;

import beast.base.inference.CalculationNode;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

/*
 * Applies the same y (macroevol param) value to all x (qu trait) bins.
 * If you put a prior on the y value, it's the same as assuming y is
 * distributed according to that prior, and independent of x
 */
public class ConstantLinkFn extends CalculationNode implements LinkFn {

    final public Input<RealParameter> yValueInput = new Input<>("yV", "Constant value of dependent variable (quantitative trait).", Input.Validate.REQUIRED);

    double yValue;
    private static final String LINKFUNCTION = "constant";

    @Override
    public void initAndValidate() {
        yValue = yValueInput.get().getValue();
    }

    @Override
    public boolean refreshParams() {

        boolean refreshedSomething = false;

        // Restored parameters may be clean while this unsnapshotted scalar still holds a proposal.
        if (yValue != yValueInput.get().getValue()) {
            yValue = yValueInput.get().getValue();
            refreshedSomething = true;
        }

        return refreshedSomething;
    }

    @Override
    public double[] getY(double[] x, double[] y, boolean ignoreRefresh) {
        boolean refreshedSomething = false;
        if (!ignoreRefresh) refreshedSomething = refreshParams();

        /*
         * if either we don't care about refreshing, or we do and something was refreshed,
         * we repopulate y
         */
        if (ignoreRefresh || refreshedSomething) {
            if (x.length != y.length) throw new RuntimeException("Sizes of x (qu trait) and y (macroevol param) differ. Exiting...");

            for (int i=0; i<x.length; i++) {
                y[i] = yValue;
            }
        }

        return y;
    }

    @Override
    public String getLinkFnName() {
        return LINKFUNCTION;
    }

}
