package drivers;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.math3.ode.*;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import SSE.InstantaneousRateMatrix;
import SSE.SSEODE;
import SSE.TraitStash;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;

public class SSEODETestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		int numberOfStates = 2; // BiSSE
		String[] spNames = new String[] { "Human" };
		Double[] lambda = new Double[] {0.222222222, 0.222222222};
		Double[] mu = new Double[] {0.0, 0.0}; // pure birth
		
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", "0.9 0.9");
		
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", "Human=2");
		double[] y = traitStash.getSpLks("Human");
		// double[] y = new double[] { 0.0, 0.0, 0.0, 1.0 };
		
		// initializing integrator and ODE
		FirstOrderIntegrator dp853 = new 
				DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		
		SSEODE ode = new SSEODE(mu, irm, 1.0, false, true, false); // incorporateCladogenesis = false (BiSSE, or MuSSE), backwardTime = true, extinctionOnly = false
		ode.setSpeciationRates(lambda);
		
		// running
		System.out.println("Initial conditions: " + Arrays.toString(y));
		for (double i = 0.0; i < 1.0; i += 0.02) {
			double end = i + 0.02;
			dp853.integrate(ode, i, y, end, y);
			System.out.println("Conditions at time " + end + ": " + Arrays.toString(y));
		}
	}
}
