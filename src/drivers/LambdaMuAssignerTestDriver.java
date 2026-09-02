package drivers;

import java.util.Arrays;

import SSE.LambdaMuAssigner;
import beast.base.inference.parameter.RealParameter;

public class LambdaMuAssignerTestDriver {

	public static void main(String[] args) {
		String lambdasToStatesString = "0,0,1,1";
		String musToStatesString = "0,1,0,1";
		
		Double lambda1 = 0.2108637; // 0A and 1A
		Double lambda2 = 0.0688838; // 0B and 1B
		Double[] lambdas = { lambda1, lambda2 };
		RealParameter lambda = new RealParameter(lambdas);
		
		Double mu1 = 0.04285281; // 0A and 1A
		Double mu2 = 1.419801E-10; // 0B and 1B
		Double[] mus = { mu1, mu2 };
		RealParameter mu = new RealParameter(mus);
		
		LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
		lambdaMuAssigner.initByName("totalNumberOfStates", 4, "nDistinctLambdas", 2, "nDistinctMus", 2, "lambdasToStates", lambdasToStatesString, "lambda", lambda, "musToStates", musToStatesString, "mu", mu);
		
		Double[] lambdaAssigned = lambdaMuAssigner.getLambdas();
		Double[] muAssigned = lambdaMuAssigner.getMus();
		
		System.out.println("Assigned lambdas: " + Arrays.toString(lambdaAssigned));
		System.out.println("Assigned mus: " + Arrays.toString(muAssigned));
	}
}
