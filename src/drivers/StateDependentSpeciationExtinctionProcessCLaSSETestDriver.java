package drivers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;

import org.apache.commons.lang3.ArrayUtils;

import SSE.CladoTriplet;
import SSE.CladogeneticSpeciationRateStash;
import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import SSE.CladoTriplet.speciationType;

public class StateDependentSpeciationExtinctionProcessCLaSSETestDriver {
	
	public static void main(String[] args) {
		
		// initializing parameter values
		int numStates = 4; // state = 1 (index 0) is 'null state'
		String[] spNames = new String[] { "Human", "Chimp", "Gorilla", "Orang" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numStates, "taxa", taxonSet, "value", "Human=2,Chimp=2,Gorilla=2,Orang=3");
		traitStash.printLksMap();
		
		double sympProb = 1.0; // DEC-like
		double subsympProb = 1.0 / 6.0;
		double vicProb = 1.0 / 6.0;
		double jProb = 0.0; // no jump dispersal
		
		double birthRate = 0.32222224;
		double deathRate = 0.1; // DEC-like
		
		Double[] mus = { deathRate, deathRate, deathRate, deathRate };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);
		
		Double[] sSpeciationRate = {sympProb * birthRate};
		Double[] ssSpeciationRate = {subsympProb * birthRate};
		Double[] vSpeciationRate = {vicProb * birthRate};
		Double[] jSpeciationRate = {jProb * birthRate}; // 0.0
		RealParameter sympatricSpeciationRate = new RealParameter(sSpeciationRate);
		RealParameter subSympatricSpeciationRate = new RealParameter(ssSpeciationRate);
		RealParameter vicariantSpeciationRate = new RealParameter(vSpeciationRate);
		RealParameter jumpSpeciationRate = new RealParameter(jSpeciationRate);
		
		CladoTriplet nullTriplet = new CladoTriplet();
		nullTriplet.initByName("parentState", 1,
				"leftChildState", 1,
				"rightChildState", 1,
				"speciationType", speciationType.SYMPATRY);
		
		CladoTriplet sTriplet1 = new CladoTriplet();
		sTriplet1.initByName("parentState", 2,
				"leftChildState", 2,
				"rightChildState", 2,
				"speciationType", speciationType.SYMPATRY);
		
		CladoTriplet sTriplet2 = new CladoTriplet();
		sTriplet2.initByName("parentState", 3,
				"leftChildState", 3,
				"rightChildState", 3,
				"speciationType", speciationType.SYMPATRY);
		
		CladoTriplet jTriplet1 = new CladoTriplet();
		jTriplet1.initByName("parentState", 2,
				"leftChildState", 2,
				"rightChildState", 3,
				"speciationType", speciationType.JUMPDISPERSAL);
		
		CladoTriplet jTriplet2 = new CladoTriplet();
		jTriplet2.initByName("parentState", 3,
				"leftChildState", 2,
				"rightChildState", 3,
				"speciationType", speciationType.JUMPDISPERSAL);
		
		CladoTriplet vTriplet1 = new CladoTriplet();
		vTriplet1.initByName("parentState", 4,
				"leftChildState", 2,
				"rightChildState", 3,
				"speciationType", speciationType.VICARIANCE);
		
		CladoTriplet ssTriplet1 = new CladoTriplet();
		ssTriplet1.initByName("parentState", 4,
				"leftChildState", 2,
				"rightChildState", 4,
				"speciationType", speciationType.SUBSYMPATRY);
		
		CladoTriplet ssTriplet2 = new CladoTriplet();
		ssTriplet2.initByName("parentState", 4,
				"leftChildState", 3,
				"rightChildState", 4,
				"speciationType", speciationType.SUBSYMPATRY);
		
		List<CladoTriplet> cladoTripletList = new ArrayList<CladoTriplet>();
		Collections.addAll(cladoTripletList, nullTriplet, sTriplet1, sTriplet2, jTriplet1, jTriplet2, vTriplet1, ssTriplet1, ssTriplet2);
		
		CladogeneticSpeciationRateStash csrt = new CladogeneticSpeciationRateStash();
		csrt.initByName("cladoTriplets", cladoTripletList,
				"sympatricRate", sympatricSpeciationRate,
				"subsympatricRate", subSympatricSpeciationRate,
				"vicariantRate", vicariantSpeciationRate,
				"jumpRate", jumpSpeciationRate);
		csrt.printEventMap();
		
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		String FlatQMatrixString = "0.0 0.0 0.0 0.01 0.0 0.01 0.01 0.0 0.01 0.0 0.01 0.01"; 
		                           
		irm.initByName("numberOfStates", numStates, "flatQMatrix", FlatQMatrixString);
		irm.printMatrix();
		
		Double[] piEs = new Double[numStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numStates];
		Arrays.fill(piDs, 1.0/(numStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "(((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0):1.0,Orang:3.0);";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        // Node root = myTree.getRoot();
        
        boolean incorporateCladogenesis = true;
        
        StateDependentSpeciationExtinctionProcess sdsep = new StateDependentSpeciationExtinctionProcess();
        sdsep.initByName(
        		"tree", myTree,
        		"traitStash", traitStash,
        		"instantaneousRateMatrix", irm,
        		"cladogeneticStash", csrt,
        		"mu", mu,
        		"pi", pi,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);

        System.out.println("Log-likelihood = " + Double.toString(sdsep.calculateLogP()));
	}
}
