package SSE;

import beast.base.inference.CalculationNode;
import beast.base.core.Input;
import beast.base.inference.util.InputUtil;
import beast.base.inference.parameter.RealParameter;

import java.util.Arrays;

/*
 * Used to initialize D's when quantitative trait is observed
 */
public class NormalCenteredAtObservedLinkFn extends CalculationNode implements LinkFn {

    final public Input<RealParameter> quTraitsInput = new Input<>("quTraits", "Quantitative trait values observed at tips", Input.Validate.REQUIRED);
    // final public Input<RealParameter> dtInput = new Input<>("dt", "Length of time interval over which integration is carried out.", Input.Validate.REQUIRED);
    // final public Input<RealParameter> diffusionInput = new Input<>("diffusion", "Diffusion term of quantitative trait diffusion process.", Input.Validate.REQUIRED);
    final public Input<RealParameter> sdNormalQuTrValueInput = new Input<>("sdNormalQuTrValue", "User-specified standard deviation of normal distribution for trait values (not trait value change!).", Input.Validate.REQUIRED);

    protected RealParameter quTraits;
    private double sdNormalQuTrValue; // not a parameter in the model
    private static final String LINKFUNCTION = "normalcenteredatdata";

    @Override
    public void initAndValidate() {
        quTraits = quTraitsInput.get();
        sdNormalQuTrValue = sdNormalQuTrValueInput.get().getValue();
    }

    @Override
    public boolean refreshParams() {

        boolean refreshedSomething = false;

        if (InputUtil.isDirty(quTraitsInput)) {
            quTraits = quTraitsInput.get();
            refreshedSomething = true;
        }

        // This scalar is not snapshotted; a clean restored parameter must also refresh it.
        if (sdNormalQuTrValue != sdNormalQuTrValueInput.get().getValue()) {
            sdNormalQuTrValue = sdNormalQuTrValueInput.get().getValue();
            refreshedSomething = true;
        }

        return refreshedSomething;
    }

    /*
     *
     */
    @Override
    public double[] getY(double[] x, double[] y, int[] nLeftNRightFlanksHi, String spName, boolean ignoreRefresh) {
        boolean refreshedSomething = false;
        if (!ignoreRefresh) refreshedSomething = refreshParams();
        int nLeftFlanks = nLeftNRightFlanksHi[0];
        int nRightFlanks = nLeftNRightFlanksHi[1];

        /*
         * if either we don't care about refreshing, or we do and something was refreshed,
         * we repopulate y
         */
        if (ignoreRefresh || refreshedSomething) {
            /*
             * for my unit test, I had to do +1 below, because my number of useful bins results from
             * subtracts an additional bin on top of the left and right flanks
             *
             * we divide y/2 because esDs array will be twice the number of useful bins to store
             * real and complex parts
             */
            // if ((x.length + nLeftFlanks + nRightFlanks + 1) != y.length/2) throw new RuntimeException("Sizes of x (qu trait) and y (esDs) differ. Exiting...");

            for (int i=0; i<x.length; i++) {
                y[i] = SSEUtils.getNormalDensity(x[i], quTraits.getValue(spName), sdNormalQuTrValue);
                // System.out.println("x["+ i + "] = " + x[i] + " : " + y[i]);
            }
        }

        return y;
    }

    @Override
    public String getLinkFnName() {
        return LINKFUNCTION;
    }

}
