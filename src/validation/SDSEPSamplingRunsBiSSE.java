package validation;

/*
Sample the tree with both drawJointConditionalAncestralState and drawStochasticCharacterMapping.
Write the posteriors to csv's
 */

//package validation;

import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;
import test.TestHelper;

import org.apache.commons.lang3.ArrayUtils;

import java.util.Arrays;
import java.util.List;

public class SDSEPSamplingRunsBiSSE {
	final static double EPSILON = 1e-2;
	final static int numTrials = 10000;
	int numTimeSlices = 100;
	private StateDependentSpeciationExtinctionProcess sdsep;

	// Builds a BiSSE model and writes joint and stochastic ancestral-state summaries.
	public void runExperiment(String treeStr, String spAttr, String[] spNames, String expName,
							  Double[] lambdas, Double[] mus, String q, int numTimeSlices) throws Exception {
		// initializing states
		int numberOfStates = 2; // BiSSE
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
		System.out.println(Arrays.toString(sdsep.getNodeIndexNameMapper()));

		// Write only the ancestral states to csv
		TestHelper.prepareAndWriteToCSV(posteriorJoint, expName + "-joint", sdsep);
		TestHelper.prepareAndWriteToCSV(posteriorStoc, expName + "-stoc", sdsep);

//		String[] divLbls = {"nd2","nd3","nd5","nd4"};
//		String[] divLks = {"0.87785842562904", "0.829891574906887", "0.00242853993720819", "0.995986779987027"};
//		String[] indexNameMapper = sdsep.getNodeIndexNameMapper();
//		TestHelper.compareDiv(divLbls, divLks, indexNameMapper, posterior);
	}

	// Runs the original small symmetric-rate validation scenario.
	public void runSmall() throws Exception {
		String treeStr = "(((sp5:1.753921149,sp6:1.753921149)nd5:10.54206596,sp2:12.2959871)nd3:5.60266132,(sp3:6.029064844,sp4:6.029064844)nd4:11.86958358)nd2;";
		String spAttr = "sp2=1,sp3=1,sp4=1,sp5=2,sp6=2";
		String[] spNames = new String[] {"sp2","sp3","sp4","sp5","sp6"};
		Double[] lambdas = {0.08, 0.08};
		Double[] mus = { 0.01, 0.01 };
		String q = "0.01 0.01";
		runExperiment(treeStr, spAttr, spNames, "beast_results/beast_small", lambdas, mus, q, numTimeSlices);
	}

	// Runs the original larger symmetric-rate validation scenario.
	public void runLarge() throws Exception {
		String treeStr = "(((((((sp23:0.1790943213,sp24:0.1790943213)nd24:1.092201813,(sp25:0.07202419731,sp26:0.07202419731)nd25:1.199271937)nd13:11.36119223,sp4:12.63248837)nd11:1.446039834,(((sp9:10.42249592,sp10:10.42249592)nd15:0.404458011,sp6:10.82695394)nd14:0.008741084797,sp5:10.83569502)nd12:3.24283318)nd8:5.258830888,(sp7:10.79794068,sp8:10.79794068)nd9:8.539418405)nd5:10.54206596,((((sp13:6.968607434,(sp21:0.8635028776,sp22:0.8635028776)nd22:6.105104556)nd19:0.7298402864,(sp15:7.029426078,sp12:7.029426078)nd20:0.6690216423)nd18:2.616431804,(sp19:1.17454796,sp20:1.17454796)nd23:9.140331564)nd10:7.268558415,sp3:17.58343794)nd6:12.2959871)nd3:5.60266132,((sp17:2.68800466,sp18:2.68800466)nd7:20.92449812,sp2:23.61250278)nd4:11.86958358)nd2;";
		String spAttr = "sp2=1,sp3=2,sp4=2,sp5=2,sp6=2,sp7=2,sp8=2,sp9=2,sp10=2,sp12=1,sp13=1,sp15=1,sp17=1,sp18=1,sp19=1,sp20=1,sp21=1,sp22=1,sp23=2,sp24=2,sp25=2,sp26=2";
		String[] spNames = new String[] {"sp2","sp3","sp4","sp5","sp6","sp7","sp8","sp9","sp10","sp12","sp13","sp15","sp17","sp18","sp19","sp20","sp21","sp22","sp23","sp24","sp25","sp26"};
		Double[] lambdas = {0.08, 0.08};
		Double[] mus = { 0.01, 0.01 };
		String q = "0.01 0.01";
		runExperiment(treeStr, spAttr, spNames, "beast_results/beast_large", lambdas, mus, q, numTimeSlices);
	}

	// Runs the first original asymmetric-rate validation scenario.
	public void runAsym() throws Exception {
		String treeStr = "((((sp8:11.27251117,((sp10:8.014869105,(sp24:0.3241612551,sp25:0.3241612551)nd19:7.69070785)nd18:0.7317424213,sp9:8.746611526)nd17:2.525899642)nd13:1.137652071,sp6:12.41016324)nd6:16.19398372,sp2:28.60414696)nd2:12.13673702,((sp14:3.316933757,sp15:3.316933757)nd7:14.26641326,((((((sp18:0.6285711163,sp19:0.6285711163)nd24:0.7208918168,sp17:1.349462933)nd14:10.58048027,((sp12:7.949237475,(sp13:4.175964768,(sp22:0.4419647162,sp23:0.4419647162)nd22:3.734000052)nd21:3.773272707)nd16:3.847572034,sp7:11.79680951)nd15:0.1331336907)nd12:0.5189043765,sp5:12.44884758)nd10:0.0122921505,(sp16:1.796173193,(sp20:0.5167068077,sp21:0.5167068077)nd23:1.279466385)nd11:10.66496653)nd9:2.7121409,sp4:15.17328063)nd8:2.410066389)nd4:23.15753696)nd1;";
		String spAttr = "sp2=1,sp4=2,sp5=2,sp6=2,sp7=2,sp8=1,sp9=2,sp10=2,sp12=1,sp13=2,sp14=2,sp15=2,sp16=2,sp17=2,sp18=2,sp19=2,sp20=2,sp21=2,sp22=1,sp23=1,sp24=2,sp25=2";
		String[] spNames = new String[]{"sp2","sp4","sp5","sp6","sp7","sp8","sp9","sp10","sp12","sp13","sp14","sp15","sp16","sp17","sp18","sp19","sp20","sp21","sp22","sp23","sp24","sp25"};
		Double[] lambdas = new Double[] {0.04, 0.08};
		Double[] mus = new Double[] {0.01, 0.02};
		String q = "0.04 0.01";
		runExperiment(treeStr, spAttr, spNames, "beast_results/beast_asym", lambdas, mus, q, numTimeSlices);
	}

	// Runs the second original asymmetric-rate validation scenario.
	public void runAsym2() throws Exception {
		String treeStr = "((((sp4:12.52481592,sp5:12.52481592)nd9:4.999507452,(sp2:17.41136713,sp3:17.41136713)nd10:0.1129562496)nd5:11.44349255,sp7:28.96781593)nd3:6.156770681,(((sp22:0.4087368897,sp23:0.4087368897)nd11:15.39647534,(((sp19:2.600802432,(sp24:0.3120853668,sp25:0.3120853668)nd24:2.288717065)nd17:7.108483465,(((sp13:6.074866948,sp14:6.074866948)nd23:0.6432900407,sp12:6.718156989)nd22:2.091402571,sp11:8.809559559)nd18:0.8997263372)nd13:2.962050882,((((sp17:3.769029668,sp18:3.769029668)nd21:5.455359347,sp10:9.224389015)nd20:0.1013655259,sp9:9.325754541)nd15:1.963129902,(sp6:9.701426779,(sp15:6.015896058,sp16:6.015896058)nd19:3.685530721)nd16:1.587457663)nd14:1.382452335)nd12:3.133875451)nd7:6.27587594,(sp20:1.817290631,sp21:1.817290631)nd8:20.26379754)nd4:13.04349844)nd2;";
		String spAttr = "sp2=1,sp3=1,sp4=1,sp5=1,sp6=1,sp7=1,sp9=1,sp10=1,sp11=1,sp12=1,sp13=1,sp14=1,sp15=1,sp16=1,sp17=1,sp18=1,sp19=1,sp20=1,sp21=1,sp22=1,sp23=1,sp24=1,sp25=1";
		String[] spNames = new String[]{"sp2","sp3","sp4","sp5","sp6","sp7","sp9","sp10","sp11","sp12","sp13","sp14","sp15","sp16","sp17","sp18","sp19","sp20","sp21","sp22","sp23","sp24","sp25"};
		Double[] lambdas = new Double[] {0.08, 0.04};
		Double[] mus = new Double[] {0.01, 0.02};
		String q = "0.001 0.01";
		runExperiment(treeStr, spAttr, spNames, "beast_results/beast_asym2", lambdas, mus, q, numTimeSlices);
	}

	// Runs the fixed validation scenario used for the diversitree comparison.
	public void runRb(String pathToWriteIn) throws Exception {
		String treeStr = "(((sp15:0.5701922606,(sp22:0.1174274481,sp23:0.1174274481)nd22:0.4527648125)nd6:5.46955786,((sp4:2.913008462,(sp16:0.4790358056,sp17:0.4790358056)nd11:2.433972656)nd9:1.72680138,sp2:4.639809842)nd7:1.399940278)nd2:8.039087646,((sp1:5.262858931,((((sp10:1.936988093,sp11:1.936988093)nd15:0.8700699862,((sp20:0.1813602217,sp21:0.1813602217)nd17:2.59756285,sp6:2.778923072)nd16:0.02813500652)nd12:0.1038009358,(sp14:1.103215563,(sp18:0.2976700868,sp19:0.2976700868)nd21:0.805545476)nd13:1.807643452)nd10:0.5229591127,sp3:3.433818127)nd8:1.829040804)nd4:1.760591904,((((sp8:1.951198056,sp9:1.951198056)nd20:0.153294648,sp7:2.104492704)nd18:0.5588707339,sp12:2.663363438)nd14:0.2401874525,sp5:2.90355089)nd5:4.119899945)nd3:7.055386931)nd1;";
		String spAttr = "sp1=1,sp2=1,sp3=1,sp4=2,sp5=1,sp6=1,sp7=2,sp8=1,sp9=1,sp10=1,sp11=1,sp12=1,sp14=1,sp15=2,sp16=2,sp17=1,sp18=1,sp19=1,sp20=2,sp21=2,sp22=1,sp23=1";
		String[] spNames = new String[]{"sp1","sp2","sp3","sp4","sp5","sp6","sp7","sp8","sp9","sp10","sp11","sp12","sp14","sp15","sp16","sp17","sp18","sp19","sp20","sp21","sp22","sp23"};
		Double[] lambdas = new Double[] {0.2, 0.4};
		Double[] mus = new Double[] {0.01, 0.1};
		String q = "0.1 0.4";
		runExperiment(treeStr, spAttr, spNames, pathToWriteIn, lambdas, mus, q, numTimeSlices);
	}

	// Runs every built-in validation scenario in sequence.
	public void run() throws Exception {
		runSmall();
		runLarge();
		runAsym();
		runAsym2();
//		runRb();
	}

	// Loads a tree and tip states from files and runs them with the fixed BiSSE rates.
	public void runUserInput(String testFile, String expName) throws Exception {
		Double[] lambdas = new Double[] {0.2, 0.4};
		Double[] mus = new Double[] {0.01, 0.1};
		String q = "0.1 0.4";

		String treeSuffix = ".tree";
		String spAttrSuffix = "-beast_str.txt";
		String tipSuffix = "-tips.csv";

		String treeStr = TestHelper.readFirstLine(testFile + treeSuffix);
		String spAttr = TestHelper.readFirstLine(testFile + spAttrSuffix);
		String tips = TestHelper.readFirstLine(testFile + tipSuffix);
		String[] spNames = tips.split(",");

		runExperiment(treeStr, spAttr, spNames, "beast_results/" + expName, lambdas, mus, q, numTimeSlices);
	}

	// Writes the fixed diversitree-comparison scenario to the requested output prefix.
	public static void main(String[] args) throws Exception {
		SDSEPSamplingRunsBiSSE runs = new SDSEPSamplingRunsBiSSE();
		runs.runRb(args[0]);
//		runs.runUserInput("validation/diversitree_data/rb", "user_input");
	}

}
