//
/*
Compare posterior probabilities at all internal nodes for drawJoint and drawStochastic
 */
package test;

import SSE.CladoTriplet;
import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import SSE.CladogeneticSpeciationRateStash;
import SSE.CladoTriplet.speciationType;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.*;

public class SDSEPSampleCLaSSEComparisonTest {
	final static double EPSILON = 1e-2;
	final static int numTrials = 10000;
	private StateDependentSpeciationExtinctionProcess sdsep;

	public void runExperiment(String treeStr, String spAttr, String[] spNames, String expName,
							  int numTimeSlices) throws Exception {
		// initializing states
		int numberOfStates = 4; // CLaSSE
		int numSpecies = spNames.length;
		int numNodes = numSpecies * 2 - 1;
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", spAttr);
		traitStash.printLksMap();

		// initializing birth-death parameters
		double sympProb = 1.0; // DEC-like
		double subsympProb = 1.0 / 6.0;
		double vicProb = 1.0 / 6.0;
		double jProb = 0.0; // no jump dispersal

		double birthRate = 0.32222224;
		double deathRate = 0.1; // DEC-like

		Double[] mus = {deathRate, deathRate, deathRate, deathRate };
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
		Arrays.fill(piDs, (1.0/numberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		// fixed equilibrium frequency
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);

		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
		// tree is given

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

		// Sample many times with drawStochasticChar and calculate the posterior
		sdsep.setNumTimeSlices(numTimeSlices);
		double[][] posteriorStoc = sdsep.sampleAndSummarizeCLaSSE(numTrials, false);
		double[][] posteriorJoint = sdsep.sampleAndSummarizeCLaSSE(numTrials, true);

//		double[] posteriorStoc1D = new double[numNodes];
//		double[] posteriorJoint1D = new double[numNodes];

		for (int i = 0; i < numNodes; i++) {
			System.out.println("Stoch for node " + i + Arrays.toString(posteriorStoc[i]));
			System.out.println("Joint for node " + i + Arrays.toString(posteriorJoint[i]));

//			posteriorStoc1D[i] = posteriorStoc[i][0];
//			posteriorJoint1D[i] = posteriorJoint[i][0];
		}

		TestHelper.compareNestedArr(posteriorJoint, posteriorStoc);

//		TestHelper.prepareAndWriteToCSV(posteriorStoc1D, "beast_results/" + expName + "-stoc", sdsep);
//		TestHelper.prepareAndWriteToCSV(posteriorJoint1D, "beast_results/" + expName + "-joint", sdsep);
	}

	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void testSmall() throws Exception {
		String treeStr = "(((Human:1.0,Chimp:1.0):1.0,Gorilla:2.0):1.0,Orang:3.0);";
		String[] spNames = new String[] { "Human", "Chimp", "Gorilla", "Orang" };
		String spAttr = "Human=3,Chimp=1,Gorilla=4,Orang=3";
		runExperiment(treeStr, spAttr, spNames, "CLaSSE", 50);
	}

	@Test
	public void testBig() throws Exception {
		String treeStr = "(((sp15:0.5701922606,(sp22:0.1174274481,sp23:0.1174274481)nd22:0.4527648125)nd6:5.46955786,((sp4:2.913008462,(sp16:0.4790358056,sp17:0.4790358056)nd11:2.433972656)nd9:1.72680138,sp2:4.639809842)nd7:1.399940278)nd2:8.039087646,((sp1:5.262858931,((((sp10:1.936988093,sp11:1.936988093)nd15:0.8700699862,((sp20:0.1813602217,sp21:0.1813602217)nd17:2.59756285,sp6:2.778923072)nd16:0.02813500652)nd12:0.1038009358,(sp14:1.103215563,(sp18:0.2976700868,sp19:0.2976700868)nd21:0.805545476)nd13:1.807643452)nd10:0.5229591127,sp3:3.433818127)nd8:1.829040804)nd4:1.760591904,((((sp8:1.951198056,sp9:1.951198056)nd20:0.153294648,sp7:2.104492704)nd18:0.5588707339,sp12:2.663363438)nd14:0.2401874525,sp5:2.90355089)nd5:4.119899945)nd3:7.055386931)nd1;";
		String[] spNames = new String[]{"sp1","sp2","sp3","sp4","sp5","sp6","sp7","sp8","sp9","sp10","sp11","sp12","sp14","sp15","sp16","sp17","sp18","sp19","sp20","sp21","sp22","sp23"};
		String spAttr = "sp1=1,sp2=1,sp3=1,sp4=2,sp5=1,sp6=1,sp7=2,sp8=1,sp9=1,sp10=1,sp11=1,sp12=1,sp14=1,sp15=2,sp16=2,sp17=1,sp18=1,sp19=1,sp20=2,sp21=2,sp22=1,sp23=1";
		runExperiment(treeStr, spAttr, spNames, "CLaSSELarge", 50);
	}

}
