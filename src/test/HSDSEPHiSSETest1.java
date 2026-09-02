package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;
import SSE.HiddenStateDependentSpeciationExtinctionProcess;
import SSE.HiddenTraitStash;
import SSE.LambdaMuAssigner;

import java.util.Arrays;
import java.util.List;
import org.apache.commons.lang3.ArrayUtils;

import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;

public class HSDSEPHiSSETest1 {
	final static double EPSILON = 1E-2;
	private double negLnl;

	@Before
	public void setUp() throws Exception {
		int numberOfStates = 2;
		int numberOfHiddenStates = 2;
		int totalNumberOfStates = numberOfStates + numberOfHiddenStates;
		String[] spNames = new String[] { "sp6", "sp10", "sp11", "sp12", "sp14", "sp15", "sp16", "sp17", "sp18", "sp19", "sp21", "sp22", "sp23", "sp24", "sp25", "sp26", "sp27", "sp28", "sp29", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42", "sp43", "sp44", "sp45", "sp46", "sp47", "sp48", "sp49", "sp50", "sp51", "sp52", "sp53", "sp54", "sp55", "sp56", "sp57", "sp58", "sp59", "sp60" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		
		String hiddenStatesString = "0,1"; // observed state 1 and 2 will transition to hidden states 1 and 2, respectively
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		
		HiddenTraitStash hiddenTraitStash = new HiddenTraitStash();
		hiddenTraitStash.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet, "hiddenObsStateMapper", stateMapper, "value", "sp6=1,sp10=1,sp11=1,sp12=1,sp14=2,sp15=1,sp16=2,sp17=2,sp18=1,sp19=1,sp21=1,sp22=1,sp23=1,sp24=1,sp25=1,sp26=2,sp27=2,sp28=1,sp29=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=1,sp36=1,sp37=2,sp38=2,sp39=1,sp40=1,sp41=1,sp42=2,sp43=2,sp44=2,sp45=2,sp46=2,sp47=2,sp48=2,sp49=1,sp50=1,sp51=2,sp52=2,sp53=1,sp54=1,sp55=2,sp56=2,sp57=1, sp58=1,sp59=1,sp60=1");
		hiddenTraitStash.printLksMap();
		
		String lambdasToStatesString = "0,1,2,3";
		Double lambda1 = 6.440928E-10; // 0A
		Double lambda2 = 0.1128975; // 1A
		Double lambda3 = 0.1308483; // 0B
		Double lambda4 = 85.03492; // 1B
		Double[] lambdas = { lambda1, lambda2, lambda3, lambda4 };
		RealParameter lambda = new RealParameter(lambdas);
		
		String musToStatesString = "0,1,2,3";
		Double mu1 = 1.417061E-09; // 0A
		Double mu2 = 0.04597049; // 1A
		Double mu3 = 0.03687651; // 0B
		Double mu4 = 1.266694E-06; // 1B
		Double[] mus = { mu1, mu2, mu3, mu4 };
		RealParameter mu = new RealParameter(mus);

		Double[] piEs = new Double[totalNumberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[totalNumberOfStates];
		Arrays.fill(piDs, (1.0/totalNumberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
				
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		
		LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
		lambdaMuAssigner.initByName("totalNumberOfStates", 4, "nDistinctLambdas", 4, "nDistinctMus", 4, "lambdasToStates", lambdasToStatesString, "lambda", lambda, "musToStates", musToStatesString, "mu", mu, "pi", pi);
		System.out.println("Lambdas: " + Arrays.toString(lambdaMuAssigner.getLambdas()));
		System.out.println("Mus: " + Arrays.toString(lambdaMuAssigner.getMus()));
		System.out.println("Pis: " + Arrays.toString(lambdaMuAssigner.getPis()));
		
		boolean disallowDoubleTransitions = true;
		int symmetrifyAcrossDiagonal = -1; // not used
		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString = "2.11773E-09 2.061154E-09 6.432635e-03 2.479683E-09 2.355085E-09 1.565876E-02 100 100";
		hirm.initByName("numberOfStates", 2, "numberOfHiddenStates", 2, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper);
		
		System.out.println("Qs:");
		hirm.printMatrix();
		
		String treeStr = "(((sp15:10.27880339,(sp57:0.4327353378,sp58:0.4327353378):9.846068053):21.30935137,((((sp49:1.322566942,sp50:1.322566942):6.531246386,(((((sp42:1.618558172,sp43:1.618558172):1.249323508,sp37:2.86788168):0.4105311845,sp36:3.278412865):1.110829025,sp28:4.38924189):2.453996398,((sp53:0.6765630317,sp54:0.6765630317):5.834067793,sp21:6.510630824):0.3326074635):1.01057504):6.546385565,sp12:14.40019889):3.891878236,((((sp18:8.595427361,((sp19:6.988162304,((sp39:1.941330272,(sp59:0.4256083779,sp60:0.4256083779):1.515721894):1.374985348,sp35:3.31631562):3.671846684):1.028692949,(sp24:5.527011086,(sp25:5.478875203,(sp40:1.898502308,sp41:1.898502308):3.580372894):0.04813588287):2.489844168):0.5785721075):0.8605508177,((sp47:1.324188282,sp48:1.324188282):1.210143714,sp38:2.534331996):6.921646183):1.848794077,(sp22:6.144323416,sp23:6.144323416):5.160448839):4.752352041,sp10:16.0571243):2.234952832):13.29607763):8.9940146,(sp6:33.80408947,(((sp29:4.271294196,sp30:4.271294196):3.963360008,(sp46:1.515605972,(sp51:0.6842469553,sp52:0.6842469553):0.8313590168):6.719048232):21.69107479,((((sp44:1.517683119,sp45:1.517683119):13.83340518,((sp33:3.451233406,sp34:3.451233406):7.318030201,sp14:10.76926361):4.581824694):2.3268441,((sp31:3.988873926,sp32:3.988873926):13.39833,(sp26:5.46221229,sp27:5.46221229):11.92499164):0.2907284735):12.10203097,((sp16:9.676541191,sp17:9.676541191):11.55054389,(sp11:16.00734921,(sp55:0.6152478573,sp56:0.6152478573):15.39210136):5.219735869):8.552878292):0.1457656227):3.878360468):6.778079891):0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
                
        HiddenStateDependentSpeciationExtinctionProcess hsdsep = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep.initByName(
        		"tree", myTree,
        		"hiddenTraitStash", hiddenTraitStash,
        		"lambdaMuAssigner", lambdaMuAssigner,
        		"hiddenInstantaneousRateMatrix", hirm,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	negLnl = hsdsep.calculateLogP();
    	System.out.println(negLnl); // -172.02167619855217
	}

	@Test
	public void againstDiversitreeHiSSE() {
		Assert.assertEquals(-172.02, negLnl, EPSILON); // difference due to precision and rounding
	}
}
