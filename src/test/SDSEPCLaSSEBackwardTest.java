/*
Given a tree of 4 tips and a CLaSSE model, we check
1. The backwards pass has a final tree likelihood matching diversitree
2. The backwards pass has a final tree likelihood matching diversitree while sampling in the branches
3. After the backwards pass, each node has the same likelihood with and without sampling in the branches
 */

package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.CladoTriplet;
import SSE.CladogeneticSpeciationRateStash;
import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import SSE.CladoTriplet.speciationType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;

public class SDSEPCLaSSEBackwardTest {
	final static double EPSILON = 1e-10;
	private StateDependentSpeciationExtinctionProcess sdsep;

	@Before
	public void setUp() throws Exception {
		// initializing states
		int numberOfStates = 4; // state = 1 (index 0) is 'null state'
		String[] spNames = new String[] { "Human", "Chimp", "Gorilla", "Orang" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", "Human=2,Chimp=2,Gorilla=2,Orang=3");
		traitStash.printLksMap();
				
		// initializing birth-death parameters
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
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", FlatQMatrixString);
		irm.printMatrix();
		
		Double[] piEs = new Double[numberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numberOfStates];
		Arrays.fill(piDs, 1.0/(numberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "(((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0):1.0,Orang:3.0);";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
		
        boolean incorporateCladogenesis = true;
        
        sdsep = new StateDependentSpeciationExtinctionProcess();
        sdsep.initByName(
        		"tree", myTree,
        		"traitStash", traitStash,
        		"instantaneousRateMatrix", irm,
        		"cladogeneticStash", csrt,
        		"mu", mu,
        		"pi", pi,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
        
        System.out.println(sdsep.calculateLogP());
	}

	private void assert2DArrayEquals(double[][] arr1, double[][] arr2) {
		Assert.assertEquals(arr1.length, arr2.length);

		for (int i = 0; i < arr1.length; i++) {
			for (int j = 0; j < arr1[i].length; j++) {
				Assert.assertEquals(arr1[i][j], arr2[i][j], 1e-7);
			}
		}
	}

	private double[][] deep2DArrayCopy(double[][] arr) {
		double[][] ret = new double[arr.length][arr[0].length];

		for (int i = 0; i < arr.length; i++) {
			System.arraycopy(arr[i], 0, ret[i], 0, arr[i].length);
		}
		return ret;
	}

	@Test
	public void test() {
		// Assert.assertEquals(-10.59346884351, sdsep.calculateLogP(), EPSILON); // Used in original version with fixed-step size ODE solver
		Assert.assertEquals(-10.59346882658, sdsep.calculateLogP(), EPSILON);
		double[][] nodePartialsNoCharHist = deep2DArrayCopy(sdsep.getNodePartialScaledLksPostOde());
		System.out.println("Passed likelihood test!");

		sdsep.setSampleCharacterHistory(true);
		Assert.assertEquals(-10.59346882658, sdsep.calculateLogP(), 1e-6);
		double[][] nodePartialsWCharHist = sdsep.getNodePartialScaledLksPostOde();
		System.out.println("Passed likelihood with sample character history test!");

		assert2DArrayEquals(nodePartialsNoCharHist, nodePartialsWCharHist);
		System.out.println("Passed matching likelihood test!");

		int[] jointNodes = sdsep.drawJointAncestralStatesAtRoot();
		int[] stocNodes = sdsep.drawStochasticAncestralStateAtRoot();
		System.out.println(Arrays.toString(jointNodes));
		System.out.println(Arrays.toString(stocNodes));

		double[][] post = sdsep.sampleAndSummarizeCLaSSE(100, false);

	}
}
