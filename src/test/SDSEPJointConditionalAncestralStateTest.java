/*
Under BiSSE model, sample the tree many times with drawJointConditionalAncestralState.
Then, check that the posterior probability of a node being in state 0 matches that given by diversitree for all nodes
 */

package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;
import test.TestHelper;

public class SDSEPJointConditionalAncestralStateTest {
	final static double EPSILON = 1e-5;
	private StateDependentSpeciationExtinctionProcess sdsep;

	@Before
	public void setUp() throws Exception {
		// initializing states
		int numberOfStates = 2; // BiSSE

		String[] spNames = new String[] {"sp1","sp2","sp3","sp4","sp5","sp6","sp7","sp8","sp9","sp10","sp11","sp12","sp14","sp15","sp16","sp17","sp18","sp19","sp20","sp21","sp22","sp23" };
		int numSpecies = spNames.length;
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", "sp1=1,sp2=1,sp3=1,sp4=2,sp5=1,sp6=1,sp7=2,sp8=1,sp9=1,sp10=1,sp11=1,sp12=1,sp14=1,sp15=2,sp16=2,sp17=1,sp18=1,sp19=1,sp20=2,sp21=2,sp22=1,sp23=1");
		traitStash.printLksMap();

		// initializing birth-death parameters
		Double[] mus = { 0.001, 0.1 };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);  // if matrix gets flattened, minor dimension is the # of rows/cols

		Double[] lambdas = {0.2, 0.4};
		RealParameter lambda = new RealParameter(lambdas);

		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", "0.1 0.4");
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

		String treeStr = "(((sp15:0.5781568251,(sp22:0.1191601926,sp23:0.1191601926)nd22:0.4589966325)nd6:5.546658275,((sp4:2.962542849,(sp16:0.485813654,sp17:0.485813654)nd11:2.476729195)nd9:1.75035681,sp2:4.712899659)nd7:1.411915441)nd2:8.269293764,((sp1:5.342102079,((((sp10:1.968953288,sp11:1.968953288)nd15:0.8858381301,((sp20:0.1840811929,sp21:0.1840811929)nd17:2.642069799,sp6:2.826150992)nd16:0.0286404258)nd12:0.1055697963,(sp14:1.120276815,(sp18:0.302320397,sp19:0.302320397)nd21:0.8179564183)nd13:1.840084399)nd10:0.5333843982,sp3:3.493745612)nd8:1.848356466)nd4:1.785661442,((((sp8:1.983406387,sp9:1.983406387)nd20:0.155819381,sp7:2.139225768)nd18:0.5692033509,sp12:2.708429119)nd14:0.2445061279,sp5:2.952935247)nd5:4.174828273)nd3:7.266345344)nd1;";
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

		// Sample many times with drawJointConditional and calculate the posterior
		sdsep.setNumTimeSlices(500);
		double[] posterior = sdsep.sampleAndSummarizeBiSSE(10000, true);
		posterior = TestHelper.trimTips(posterior);

		// Write only the ancestral states to csv
//		String dir = "/Users/jeff/Documents/Research/Phylogenetics/SSE/validation/beast_results";
//		TestHelper.writeToCSV(dir, "rb.csv", posterior, sdsep);

		String[] divLbls = {"nd1","nd2","nd6","nd22","nd7","nd9","nd11","nd3","nd4","nd8","nd10","nd12","nd15","nd16","nd17","nd13","nd21","nd5","nd14","nd18","nd20"};
		String[] divLks = {"0.504778971188059","0.761371584451414","0.546427028721335","0.990929751746494","0.700705453757888","0.520712315675997","0.400707333658542","0.798875159009061","0.773218277930007","0.693681343224611","0.651616182147769","0.641510203760138","0.801524137771315","0.639132464851372","0.00104369318998754","0.939501019441591","0.99509883391553","0.725987011544158","0.706256765576992","0.654024903986533","0.693527170575002"};
		HashMap<String, Double> divMap = TestHelper.getDivMap(divLbls, divLks);
		String[] indexNameMapper = sdsep.getNodeIndexNameMapper();
		TestHelper.compareDivPosterior(divMap, indexNameMapper, posterior);
		System.out.println("Passed drawJointCond test - posteriors calculated from sampling is the same as diversitree!");
	}


	@Test
	public void test() {
		Assert.assertEquals(-63.26608, sdsep.calculateLogP(), EPSILON); // Used in original version with fixed-step size ODE solver
		System.out.println("Passed logP test - backwards pass is good!");

//		double[] lks = {0, 0, 0.4, 0.6};
//		int numTrials = 10000;
//		int count = 0;
//		for (int i = 0; i < numTrials; i++) {
//			int foo = sdsep.sampleLksArray(lks);
//			if (foo == 0) {
//				count++;
//			}
//		}
//		Assert.assertEquals(0.4, 1.0 * count / numTrials, 0.001);

	}

}
