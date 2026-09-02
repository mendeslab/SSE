package mosse;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.distribution.NormalDistribution;


/**
 * @author Kylie Chen
 */

@Description("Mosse tip likelihood using regression model of tip rate given traits and parameters beta and epsilon")
public class MosseTipLikelihood extends CalculationNode {

	final public Input<RealParameter> betaInput = new Input<>("beta", "beta coefficients for each trait", Input.Validate.REQUIRED);
	final public Input<RealParameter> epsilonInput = new Input<>("epsilon", "error term of regression model", Input.Validate.REQUIRED);
	final public Input<BooleanParameter> logScaleInput = new Input<>("logscale", "whether to use log scale for substitution rate (defaults to true)", Input.Validate.OPTIONAL);

	private RealParameter beta;
	private RealParameter epsilon;
	private boolean logScale;

	public MosseTipLikelihood() {

	}

	/**
	 * Returns P(tip rate | beta0, beta1, epsilon, trait0, trait1) ~
	 * Gaussian(mean = beta0 * trait0 + beta1 * trait1 + ... + epsilon, sd = epsilon)
	 * within the tip rate interval (a,b)
	 * @param a start value of tip rate interval
	 * @param b end value of tip rate interval
	 * @param traits trait values
	 * @return tip likelihood between intervals a and b
	 */
	public double getTipLikelihood(double a, double b, double[] traits) {
		double mean = epsilon.getValue();
		for (int i = 0; i < traits.length; i++) {
			int numBetas = beta.getDimension();
			if (numBetas != traits.length) {
				throw new IllegalArgumentException("beta dimension not equal to trait dimension!");
			}
			mean += beta.getValue(i) * traits[i];
		}
		double sd = epsilon.getValue();
		// Gaussian distribution
		NormalDistribution normalDist = new NormalDistribution(mean, sd);
		double startProb = normalDist.cumulativeProbability(a);
		double endProb = normalDist.cumulativeProbability(b);
		return endProb - startProb;
	}

	/**
	 *
	 * @param traits array of trait values
	 * @param numBins number of bins for substitution rate discretization
	 * @param startSubsRate substitution rate lower bound
	 * @param endSubsRate substitution rate upper bound
	 * @return array of tip likelihoods
	 */
	public double[] getTipLikelihoods(double[] traits, int numBins, double startSubsRate, double endSubsRate) {
		double subsInterval = (endSubsRate - startSubsRate) / numBins;
		double[] tipLikelihoods = new double[numBins];
		for (int i = 0; i < numBins; i++) {
			double a = startSubsRate + i * subsInterval;
			double b = startSubsRate + (i + 1) * subsInterval;
			if (logScale) {
				a = Math.log(a);
				b = Math.log(b);
			}
			tipLikelihoods[i] = getTipLikelihood(a, b, traits);
		}
		return tipLikelihoods;
	}

	@Override
	public void initAndValidate() {
		beta = betaInput.get();
		epsilon = epsilonInput.get();
		if (logScaleInput.get() == null) {
			logScale = true;
		} else {
			logScale = logScaleInput.get().getValue();
		}
	}
}