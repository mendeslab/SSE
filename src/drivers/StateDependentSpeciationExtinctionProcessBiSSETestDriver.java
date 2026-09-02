package drivers;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;

public class StateDependentSpeciationExtinctionProcessBiSSETestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		int numberOfStates = 2; // BiSSE
		String[] spNames = new String[] { "Human", "Chimp", "Gorilla" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", "Human=2,Chimp=2,Gorilla=2");
		traitStash.printLksMap();
		
		Double birthRate = 0.222222222;
		Double deathRate = 0.1;
		Double[] mus = { deathRate, deathRate };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);
		
		Double[] lambdas = new Double[numberOfStates];
		Arrays.fill(lambdas, birthRate);
		RealParameter lambda = new RealParameter(lambdas);
		
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", "0.1 0.1");
		irm.printMatrix();
		
		Double[] piEs = new Double[numberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numberOfStates];
		Arrays.fill(piDs, (1.0/numberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.5, 0.5
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
                
        StateDependentSpeciationExtinctionProcess sdsep = new StateDependentSpeciationExtinctionProcess();
        sdsep.initByName(
        		"tree", myTree,
        		"traitStash", traitStash,
        		"instantaneousRateMatrix", irm,
        		"lambda", lambda,
        		"mu", mu,
        		"pi", pi,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	System.out.println(sdsep.calculateLogP());
	}
}
