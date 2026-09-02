package mosse;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeDistribution;

import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.lang.UnsupportedOperationException;


/**
 * @author Kylie Chen
 */

@Description("Mosse tree model")
public class MosseDistribution extends TreeDistribution {

    final public Input<IntegerParameter> nxInput = new Input<>("nx", "number of bins for substitution rate", new IntegerParameter("1024"));
    final public Input<RealParameter> dxInput = new Input<>("dx", "distance between xs", new RealParameter("0.0001"));
    final public Input<List<IntegerParameter>> ndInput = new Input<>("nd", "plan dimensions for fftw3", Arrays.asList(new IntegerParameter("5")));
    final public Input<RealParameter> driftInput = new Input<>("drift", "drift parameter", new RealParameter("0.0"));
    final public Input<RealParameter> diffusionInput = new Input<>("diffusion", "diffusion parameter", new RealParameter("0.001"));
    final public Input<RealParameter> dtInput = new Input<>("dt", "time interval dt", new RealParameter("0.01"));
    final public Input<IntegerParameter> padLeftInput = new Input<>("padleft", "number of elements to pad left of kernel", new IntegerParameter("40"));
    final public Input<IntegerParameter> padRightInput = new Input<>("padright", "number of elements to pad right of kernel", new IntegerParameter("40"));

    final public int FLAG_FFTW3_DEFAULT = 0;
    final public int NUCLEOTIDE_STATES_DEFAULT = 4; // we can get this value from datatype

    static {
        System.loadLibrary("test");
    }

    /**
     * initialize mosse C object
     * @param nx number of bins for substitution rate
     * @param dx distance between xs
     * @param array_nd plan dimensions for FFTW3 integration
     * @param flags flags for FFTW3 integration
     * @return mosse object pointer
     */
    private native long makeMosseFFT(int nx, double dx, int[] array_nd, int flags);

    /**
     * destroy mosse C object
     * @param obj_ptr mosse object pointer
     */
    private native void mosseFinalize(long obj_ptr);

    private native double[] doIntegrateMosse(long obj_ptr, double[] vars, double[] lambda, double[] mu,
            double drift, double diffusion, double[] Q, int nt, double dt, int pad_left, int pad_right);

    private double[][] eProbs;
    private double[][] dProbs;

    public void initAndValidate() {

    }

    public double calculateLogP(Tree tree) {
        double logP = 0.0;
        return logP;
    }

    public double calculateBranchLogP(double[] vars, double[] lambda, double[] mu, double[] Q) {
        double logP = 0.0;
        // getting parameter values
        int nx = nxInput.get().getValue();
        double dx = dxInput.get().getValue();
        int[] nd = new int[ndInput.get().size()];
        ndInput.get().toArray();
        double drift = driftInput.get().getValue();
        double diffusion = diffusionInput.get().getValue();
        double dt = dtInput.get().getValue();
        int padLeft = padLeftInput.get().getValue();
        int padRight = padRightInput.get().getValue();

        double[] result = doIntegration(nx, dx, nd, FLAG_FFTW3_DEFAULT,
            vars, lambda, mu,
            drift, diffusion,
            Q, NUCLEOTIDE_STATES_DEFAULT, dt,
            padLeft, padRight);
        int ncol = vars.length / nx;
        double[][] ans = new double[nx][ncol];
        logP = calculateBranchLogP(result, nx, ncol, dx, ans);
        return logP;
    }

    /**
     * calculate the log probability on a single branch
     * @param array input non logged partials from doIntegration()
     * @param nx number of bins for substitution rate
     * @param ncol number of columns
     * @param dx distance between xs
     * @param ans array for storing logged result
     * @return log probability
     */
    public double calculateBranchLogP(double[] array, int nx, int ncol, double dx, double[][] ans) {
        double logP = logCompensation(nx, ncol, dx, array, ans);
        return logP; // return log compensated result
    }

    /**
     * perform integration along a single branch
     * @param nx number of bins for substitution rate
     * @param dx distance between xs
     * @param nd plan dimensions for FFT3 integration
     * @param flags flags for FFTW3 integration
     * @param vars array of tips or partials
     * @param lambda array of birth-rates
     * @param mu array of death-rates
     * @param drift drift parameter
     * @param diffusion diffusion parameter
     * @param Q exponentiated form of the Q matrix
     * @param nt number of nucleotide types (default = 4)
     * @param dt_max maximum time step
     * @param pad_left padding size left of the kernel (zero padding)
     * @param pad_right padding size right of the kernel (zero padding)
     * @return partial probabilities
     */
    public double[] doIntegration(int nx, double dx, int[] nd, int flags,
                              double[] vars, double[] lambda, double[] mu,
                              double drift, double diffusion,
                              double[] Q, int nt, double dt_max,
                              int pad_left, int pad_right) {
        // make mosse fft object pointer
        long ptr = makeMosseFFT(nx, dx, nd, flags);
        // integrate using C propagate x and propagate t
        double[] result = doIntegrateMosse(ptr, vars, lambda, mu, drift, diffusion, Q, nt, dt_max, pad_left, pad_right);
        mosseFinalize(ptr); // destroy obj pointer
        return result; // return non logged results
    }

    /**
     *
     * @param nrow number of rows
     * @param ncol number of columns
     * @param dx distance between xs
     * @param result non-logged partial probabilities
     * @param ans logged partial probabilities
     * @return logged partial probabilities
     */
    public double logCompensation(int nrow, int ncol, double dx, double[] result, double[][] ans) {
        double logP = 0.0;
        int count = 0;
        for (int j = 0; j < ncol; j++) {
            for (int i = 0; i < nrow; i++) {
                ans[i][j] = result[count];
                count++;
            }
        }
        if (ncol > 1 ) {
            double sum = 0.0;
            // sum except first col
            for (int i = 0; i < nrow; i++) {
                for (int j = 1; j < ncol; j++) {
                    sum += ans[i][j] * dx;
                }
            }
            double q = sum;
            // update ans except first col
            for (int i = 0; i < nrow; i++) {
                for (int j = 1; j < ncol; j++) {
                    ans[i][j] = ans[i][j] / q;
                }
            }
            logP = Math.log(q);
        } else {
            logP = 0.0;
        }
        return logP;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
        throw new UnsupportedOperationException();
    }

}
