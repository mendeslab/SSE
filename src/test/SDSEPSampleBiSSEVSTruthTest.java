/*
Under BiSSE model, sample the tree many times with drawJointConditionalAncestralState and drawStochasticCharacterMapping.
For each node, determine the most frequently visited state. Compare this with the ground truth from diversitree simulation
Report accuracy - What proportion of the time do we sample the ground truth state
IMPORTANT: Make sure the div data here is TRUE STATES not POSTERIORS
 */

package test;

import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;
import org.apache.commons.lang3.ArrayUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

public class SDSEPSampleBiSSEVSTruthTest {
	final static double EPSILON = 1e-2;
	final static int numTrials = 10000;
	private StateDependentSpeciationExtinctionProcess sdsep;

	public void runExperiment(String treeStr, String spAttr, String[] spNames, String expName,
							  Double[] lambdas, Double[] mus, String q, int numTimeSlices,
							  String[] divLbls, String[] divStates, double divAcc) throws Exception {
		// initializing states
		int numberOfStates = 2; // BiSSE
		int numSpecies = spNames.length;
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", spAttr);
		traitStash.printLksMap();

		// initializing birth-death parameters
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);  // if matrix gets flattened, minor dimension is the # of rows/cols

		RealParameter lambda = new RealParameter(lambdas);

		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", q);
		irm.printMatrix();

		Double[] piEs = new Double[numberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numberOfStates];
		Arrays.fill(piDs, (1.0/numberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.5, 0.5
		// fixed equilibrium frequency
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);

		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
		// tree is given

		boolean incorporateCladogenesis = false;

		sdsep = new StateDependentSpeciationExtinctionProcess();
		sdsep.initByName(
				"tree", myTree,
				"traitStash", traitStash,
				"instantaneousRateMatrix", irm,
				"lambda", lambda,
				"mu", mu,
				"pi", pi,
				"incorporateCladogenesis", incorporateCladogenesis
		);

		// Sample many times with drawStochasticChar and calculate the posterior
		sdsep.setNumTimeSlices(numTimeSlices);
		double[] posteriorStoc = sdsep.sampleAndSummarizeBiSSE(numTrials, false);
		double[] posteriorJoint = sdsep.sampleAndSummarizeBiSSE(numTrials, true);

		// Write only the ancestral states to csv
		TestHelper.prepareAndWriteToCSV(posteriorJoint, expName + "-joint", sdsep);
		TestHelper.prepareAndWriteToCSV(posteriorStoc, expName + "-stoc", sdsep);

		int[] posteriorStocMostCommon = TestHelper.parsimony(posteriorStoc);
		int[] posteriorJointMostCommon = TestHelper.parsimony(posteriorJoint);

		posteriorStocMostCommon = TestHelper.trimTipsInt(posteriorStocMostCommon);
		posteriorJointMostCommon = TestHelper.trimTipsInt(posteriorJointMostCommon);


		HashMap<String, Double> divMap = TestHelper.getDivMap(divLbls, divStates);
		String[] idxLabelMapper = sdsep.getNodeIndexNameMapper();
		double accStoc = TestHelper.compareDivTruth(divMap, idxLabelMapper, posteriorStocMostCommon);
        double accJoint = TestHelper.compareDivTruth(divMap, idxLabelMapper, posteriorJointMostCommon);
		System.out.println(Arrays.toString(idxLabelMapper));
        System.out.println("Stoc accuracy: " + accStoc);
		System.out.println("Joint accuracy: " + accJoint);
		System.out.println("Div accuracy: " + divAcc);
		Assert.assertTrue(accStoc + EPSILON > divAcc);
        System.out.println("Draw stochastic mapping gets accuracy close to diversitree asr");
		Assert.assertTrue(accJoint + EPSILON > divAcc);
		System.out.println("Draw joint gets accuracy close to diversitree asr");
	}
	@Before
	public void setUp() throws Exception {
	}

	@Test
	public void test() throws Exception {
		// Small
		String treeStr = "(((sp5:1.753921149,sp6:1.753921149)nd5:10.54206596,sp2:12.2959871)nd3:5.60266132,(sp3:6.029064844,sp4:6.029064844)nd4:11.86958358)nd2;";
		String spAttr = "sp2=1,sp3=1,sp4=1,sp5=2,sp6=2";
		String[] spNames = new String[] {"sp2","sp3","sp4","sp5","sp6"};
		Double[] lambdas = {0.08, 0.08};
		Double[] mus = { 0.01, 0.01 };
		String q = "0.01 0.01";
		String[] divLbls = {"nd2", "nd3", "nd5", "nd4"};
		String[] divStates = {"0", "0", "1", "0"};
		double divAcc = 1;
		runExperiment(treeStr, spAttr, spNames, "beast_small", lambdas, mus, q, 500, divLbls, divStates, divAcc);
		Assert.assertEquals(-17.97332, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
		sdsep.setSampleCharacterHistory(true);
		Assert.assertEquals(-17.97332, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
        System.out.println("Small test passed");


        // RB
		treeStr = "(((sp15:0.5701922606,(sp22:0.1174274481,sp23:0.1174274481)nd22:0.4527648125)nd6:5.46955786,((sp4:2.913008462,(sp16:0.4790358056,sp17:0.4790358056)nd11:2.433972656)nd9:1.72680138,sp2:4.639809842)nd7:1.399940278)nd2:8.039087646,((sp1:5.262858931,((((sp10:1.936988093,sp11:1.936988093)nd15:0.8700699862,((sp20:0.1813602217,sp21:0.1813602217)nd17:2.59756285,sp6:2.778923072)nd16:0.02813500652)nd12:0.1038009358,(sp14:1.103215563,(sp18:0.2976700868,sp19:0.2976700868)nd21:0.805545476)nd13:1.807643452)nd10:0.5229591127,sp3:3.433818127)nd8:1.829040804)nd4:1.760591904,((((sp8:1.951198056,sp9:1.951198056)nd20:0.153294648,sp7:2.104492704)nd18:0.5588707339,sp12:2.663363438)nd14:0.2401874525,sp5:2.90355089)nd5:4.119899945)nd3:7.055386931)nd1;";
		spAttr = "sp1=1,sp2=1,sp3=1,sp4=2,sp5=1,sp6=1,sp7=2,sp8=1,sp9=1,sp10=1,sp11=1,sp12=1,sp14=1,sp15=2,sp16=2,sp17=1,sp18=1,sp19=1,sp20=2,sp21=2,sp22=1,sp23=1";
		spNames = new String[]{"sp1","sp2","sp3","sp4","sp5","sp6","sp7","sp8","sp9","sp10","sp11","sp12","sp14","sp15","sp16","sp17","sp18","sp19","sp20","sp21","sp22","sp23"};
		lambdas = new Double[] {0.2, 0.4};
		mus = new Double[] {0.01, 0.1};
		q = "0.1 0.4";
		divLbls = new String[] {"nd1","nd2","nd6","nd22","nd7","nd9","nd11","nd3","nd4","nd8","nd10","nd12","nd15","nd16","nd17","nd13","nd21","nd5","nd14","nd18","nd20"};
		divStates = new String[] {"0","1","1","0","1","1","1","0","0","0","0","0","0","0","1","0","0","0","0","0","0"};
		divAcc = 0.7619048;
		runExperiment(treeStr, spAttr, spNames, "beast_rb", lambdas, mus, q, 500, divLbls, divStates, divAcc);
		Assert.assertEquals(-63.0014, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
		sdsep.setSampleCharacterHistory(true);
		Assert.assertEquals(-63.0014, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver


        // Test2
		treeStr = "((((sp13:1.091879977,sp14:1.091879977)nd18:0.6759211435,sp5:1.76780112)nd6:2.772232755,(((sp22:0.05611638081,sp23:0.05611638081)nd22:0.2713390715,sp19:0.3274554523)nd10:2.797675452,sp1:3.125130904)nd7:1.414902971)nd2:4.225862134,(((sp15:0.6802057993,sp16:0.6802057993)nd8:3.035233042,(((sp7:1.683835836,sp8:1.683835836)nd19:0.09043471409,(sp17:0.4689228608,sp18:0.4689228608)nd16:1.30534769)nd12:0.5383307309,(sp20:0.07351713161,sp21:0.07351713161)nd13:2.23908415)nd9:1.40283756)nd4:1.522763703,(sp2:2.860986153,(sp3:2.074989102,((((sp11:1.102714652,sp12:1.102714652)nd21:0.1185607974,sp10:1.221275449)nd20:0.3858338004,sp9:1.60710925)nd17:0.1656769653,sp4:1.772786215)nd14:0.3022028868)nd11:0.7859970509)nd5:2.377216391)nd3:3.527693465)nd1;";
		spAttr = "sp1=2,sp2=2,sp3=2,sp4=2,sp5=2,sp7=2,sp8=2,sp9=2,sp10=2,sp11=2,sp12=2,sp13=2,sp14=2,sp15=1,sp16=1,sp17=2,sp18=2,sp19=2,sp20=2,sp21=2,sp22=2,sp23=2";
		spNames = new String[]{"sp1","sp2","sp3","sp4","sp5","sp7","sp8","sp9","sp10","sp11","sp12","sp13","sp14","sp15","sp16","sp17","sp18","sp19","sp20","sp21","sp22","sp23"};
		lambdas = new Double[] {0.5, 0.4};
		mus = new Double[] {0.02, 0.1};
		q = "0.1 0.02";
		divLbls = new String[] {"nd1","nd2","nd6","nd18","nd7","nd10","nd22","nd3","nd4","nd8","nd9","nd12","nd19","nd16","nd13","nd5","nd11","nd14","nd17","nd20","nd21"};
		divStates = new String[] {"0","1","1","1","1","1","1","0","0","0","1","1","1","1","1","1","1","1","1","1","1"};
		divAcc = 0.8571429;
		runExperiment(treeStr, spAttr, spNames, "test2", lambdas, mus, q, 500, divLbls, divStates, divAcc);

		Assert.assertEquals(-46.09716, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
		sdsep.setSampleCharacterHistory(true);
		Assert.assertEquals(-46.09716, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver


		//Asym
		treeStr = "((((sp8:11.27251117,((sp10:8.014869105,(sp24:0.3241612551,sp25:0.3241612551)nd19:7.69070785)nd18:0.7317424213,sp9:8.746611526)nd17:2.525899642)nd13:1.137652071,sp6:12.41016324)nd6:16.19398372,sp2:28.60414696)nd2:12.13673702,((sp14:3.316933757,sp15:3.316933757)nd7:14.26641326,((((((sp18:0.6285711163,sp19:0.6285711163)nd24:0.7208918168,sp17:1.349462933)nd14:10.58048027,((sp12:7.949237475,(sp13:4.175964768,(sp22:0.4419647162,sp23:0.4419647162)nd22:3.734000052)nd21:3.773272707)nd16:3.847572034,sp7:11.79680951)nd15:0.1331336907)nd12:0.5189043765,sp5:12.44884758)nd10:0.0122921505,(sp16:1.796173193,(sp20:0.5167068077,sp21:0.5167068077)nd23:1.279466385)nd11:10.66496653)nd9:2.7121409,sp4:15.17328063)nd8:2.410066389)nd4:23.15753696)nd1;";
		spAttr = "sp2=1,sp4=2,sp5=2,sp6=2,sp7=2,sp8=1,sp9=2,sp10=2,sp12=1,sp13=2,sp14=2,sp15=2,sp16=2,sp17=2,sp18=2,sp19=2,sp20=2,sp21=2,sp22=1,sp23=1,sp24=2,sp25=2";
		spNames = new String[]{"sp2","sp4","sp5","sp6","sp7","sp8","sp9","sp10","sp12","sp13","sp14","sp15","sp16","sp17","sp18","sp19","sp20","sp21","sp22","sp23","sp24","sp25"};
		lambdas = new Double[] {0.04, 0.08};
		mus = new Double[] {0.01, 0.02};
		q = "0.04 0.01";
		divLbls = new String[] {"nd1","nd2","nd6","nd13","nd17","nd18","nd19","nd4","nd7","nd8","nd9","nd10","nd12","nd14","nd24","nd15","nd16","nd21","nd22","nd11","nd23"};
		divStates = new String[] {"1","0","1","1","1","1","1","1","1","1","1","1","1","1","1","1","0","0","0","1","1"};
		divAcc = 0.8571429;
		runExperiment(treeStr, spAttr, spNames, "beast_asym", lambdas, mus, q, 500, divLbls, divStates, divAcc);
		Assert.assertEquals(-84.5913, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
		sdsep.setSampleCharacterHistory(true);
		Assert.assertEquals(-84.5913, sdsep.calculateLogP(), 1e-3); // Used in original version with fixed-step size ODE solver
	}

}
