package drivers;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;
import SSE.HiddenStateDependentSpeciationExtinctionProcess;
import SSE.HiddenTraitStash;
import SSE.LambdaMuAssigner;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;

public class HiddenStateDependentSpeciationExtinctionProcessHiSSETestDriver {

	public static void main(String[] args) {
		
		// initializing parameter values
		int numberOfStates = 2;
		int numberOfHiddenStates1 = 2;
		int numberOfHiddenStates2 = 1;
		int totalNumberOfStates1 = numberOfStates + numberOfHiddenStates1;
		int totalNumberOfStates2 = numberOfStates + numberOfHiddenStates2;
		String[] spNames1 = new String[] { "sp6", "sp10", "sp11", "sp12", "sp14", "sp15", "sp16", "sp17", "sp18", "sp19", "sp21", "sp22", "sp23", "sp24", "sp25", "sp26", "sp27", "sp28", "sp29", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42", "sp43", "sp44", "sp45", "sp46", "sp47", "sp48", "sp49", "sp50", "sp51", "sp52", "sp53", "sp54", "sp55", "sp56", "sp57", "sp58", "sp59", "sp60" }; // test 1
		String[] spNames2 = new String[] { "sp1", "sp10", "sp11", "sp12", "sp14", "sp15", "sp16", "sp18", "sp19", "sp20", "sp22", "sp24", "sp25", "sp26", "sp27", "sp28", "sp29", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42" }; // test 2
		List<Taxon> taxaList1 = Taxon.createTaxonList(Arrays.asList(spNames1));
		List<Taxon> taxaList2 = Taxon.createTaxonList(Arrays.asList(spNames2));
		TaxonSet taxonSet1 = new TaxonSet(taxaList1);
		TaxonSet taxonSet2 = new TaxonSet(taxaList2);
		
		String hiddenStatesString1 = "0,1"; // observed state 1 and 2 will transition to hidden states 1 and 2, respectively
		HiddenObservedStateMapper stateMapper1 = new HiddenObservedStateMapper();
		stateMapper1.initByName("hiddenStates", hiddenStatesString1);
		stateMapper1.makeMaps();

		HiddenTraitStash hiddenTraitStash1 = new HiddenTraitStash();
		hiddenTraitStash1.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates1, "taxa", taxonSet1, "hiddenObsStateMapper", stateMapper1, "value", "sp6=1,sp10=1,sp11=1,sp12=1,sp14=2,sp15=1,sp16=2,sp17=2,sp18=1,sp19=1,sp21=1,sp22=1,sp23=1,sp24=1,sp25=1,sp26=2,sp27=2,sp28=1,sp29=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=1,sp36=1,sp37=2,sp38=2,sp39=1,sp40=1,sp41=1,sp42=2,sp43=2,sp44=2,sp45=2,sp46=2,sp47=2,sp48=2,sp49=1,sp50=1,sp51=2,sp52=2,sp53=1,sp54=1,sp55=2,sp56=2,sp57=1, sp58=1,sp59=1,sp60=1");
		hiddenTraitStash1.printLksMap();
		
		String hiddenStatesString2 = "-1,0"; // observed state 1 and 2 will transition to no hidden state, and hidden states 1, respectively
		HiddenObservedStateMapper stateMapper2 = new HiddenObservedStateMapper();
		stateMapper2.initByName("hiddenStates", hiddenStatesString2);
		stateMapper2.makeMaps();
		
		HiddenTraitStash hiddenTraitStash2 = new HiddenTraitStash();
		hiddenTraitStash2.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates2, "taxa", taxonSet2, "hiddenObsStateMapper", stateMapper2, "value", "sp1=2,sp10=1,sp11=1,sp12=1,sp14=1,sp15=2,sp16=1,sp18=1,sp19=1,sp20=2,sp22=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp29=2,sp30=1,sp31=1,sp32=2,sp33=1,sp34=1,sp35=2,sp36=2,sp37=2,sp38=2,sp39=2,sp40=2,sp41=2,sp42=2");
		hiddenTraitStash2.printLksMap();
		
		String lambdasToStatesString1 = "0,1,2,3"; // first lambda to first state, second lambda to second state, and so no...
		Double lambda1 = 6.440928E-10;
		Double lambda2 = 0.1128975;
		Double lambda3 = 0.1308483;
		Double lambda4 = 85.03492;
		Double[] lambdas1 = { lambda1, lambda2, lambda3, lambda4 }; // 0A, 1A, 0B, 1B
		System.out.println("Lambdas test 1: " + Arrays.toString(lambdas1));
		
		String lambdasToStatesString2 = "0,1,2";
		lambda1 = 0.08885618;
		lambda2 = 0.1509081;
		lambda3 = 3.597251;
		Double[] lambdas2 = { lambda1, lambda2, lambda3 }; // 0A, 1A, 1B
		System.out.println("Lambdas test 2: " + Arrays.toString(lambdas2));
		
		String musToStatesString1 = "0,1,2,3"; // first mu to first state, second mu to second state
		Double mu1 = 1.417061E-09;
		Double mu2 = 0.04597049;
		Double mu3 = 0.03687651;
		Double mu4 = 1.266694E-06;
		Double[] mus1 = { mu1, mu2, mu3, mu4 };		
		System.out.println("Mus test 1: " + Arrays.toString(mus1));
		
		String musToStatesString2 = "0,1,2";
		mu1 = 0.0224435;
		mu2 = 0.05379272;
		mu3 = 1.809666e-08;
		Double[] mus2 = { mu1, mu2, mu3 };
		System.out.println("Mus test 2: " + Arrays.toString(mus2));
		
		RealParameter muTest1 = new RealParameter(mus1);
		RealParameter muTest2 = new RealParameter(mus2);
		muTest1.initByName("minordimension", 1);
		muTest2.initByName("minordimension", 1);
		
		RealParameter lambdaTest1 = new RealParameter(lambdas1);
		RealParameter lambdaTest2 = new RealParameter(lambdas2);

		boolean disallowDoubleTransitions = true;
		HiddenInstantaneousRateMatrix hirm1 = new HiddenInstantaneousRateMatrix();
		HiddenInstantaneousRateMatrix hirm2 = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString1 = "2.11773E-09 2.061154E-09 6.432635e-03 2.479683E-09 2.355085E-09 1.565876E-02 100 100"; // test 1
		String flatQMatrixString2 = "0.02597671 0.04404249 2.061154E-09 4.59653214"; // test 2
		
		hirm1.initByName("numberOfStates", 2, "numberOfHiddenStates", 2, "flatQMatrix", flatQMatrixString1, "disallowDoubleTransitions", disallowDoubleTransitions, "hiddenObsStateMapper", stateMapper1, "symmetrifyAcrossDiagonal", -1);
		hirm1.printMatrix();
		hirm2.initByName("numberOfStates", 2, "numberOfHiddenStates", 1, "flatQMatrix", flatQMatrixString2, "disallowDoubleTransitions", disallowDoubleTransitions, "hiddenObsStateMapper", stateMapper2, "symmetrifyAcrossDiagonal", -1); 
		hirm2.printMatrix();
		
		Double[] piEs1 = new Double[totalNumberOfStates1];
		Arrays.fill(piEs1, 0.0);
		Double[] piDs1 = new Double[totalNumberOfStates1];
		Arrays.fill(piDs1, (1.0/totalNumberOfStates1));
		Double[] pis1 = ArrayUtils.addAll(piEs1, piDs1); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		
		Double[] piEs2 = new Double[totalNumberOfStates2];
		Arrays.fill(piEs2, 0.0);
		Double[] piDs2 = new Double[totalNumberOfStates2];
		Arrays.fill(piDs2, (1.0/totalNumberOfStates2));
		Double[] pis2 = ArrayUtils.addAll(piEs2, piDs2); // 0.0, 0.0, 0.0, 0.333, 0.333, 0.333
		
		System.out.println("Pi test 1 is: " + Arrays.toString(pis1));
		RealParameter pi1 = new RealParameter(pis1);
		pi1.initByName("minordimension", 1);
		
		System.out.println("Pi test 2 is: " + Arrays.toString(pis2));
		RealParameter pi2 = new RealParameter(pis2);
		pi2.initByName("minordimension", 1);
		
		LambdaMuAssigner lambdaMuAssigner1 = new LambdaMuAssigner();
		lambdaMuAssigner1.initByName("totalNumberOfStates", 4, "nDistinctLambdas", 4, "nDistinctMus", 4, "lambdasToStates", lambdasToStatesString1, "lambda", lambdaTest1, "musToStates", musToStatesString1, "mu", muTest1, "pi", pi1);

		LambdaMuAssigner lambdaMuAssigner2 = new LambdaMuAssigner();
		lambdaMuAssigner2.initByName("totalNumberOfStates", 3, "nDistinctLambdas", 3, "nDistinctMus", 3, "lambdasToStates", lambdasToStatesString2, "lambda", lambdaTest2, "musToStates", musToStatesString2, "mu", muTest2, "pi", pi2);
		
		String treeStr1 = "(((sp15:10.27880339,(sp57:0.4327353378,sp58:0.4327353378):9.846068053):21.30935137,((((sp49:1.322566942,sp50:1.322566942):6.531246386,(((((sp42:1.618558172,sp43:1.618558172):1.249323508,sp37:2.86788168):0.4105311845,sp36:3.278412865):1.110829025,sp28:4.38924189):2.453996398,((sp53:0.6765630317,sp54:0.6765630317):5.834067793,sp21:6.510630824):0.3326074635):1.01057504):6.546385565,sp12:14.40019889):3.891878236,((((sp18:8.595427361,((sp19:6.988162304,((sp39:1.941330272,(sp59:0.4256083779,sp60:0.4256083779):1.515721894):1.374985348,sp35:3.31631562):3.671846684):1.028692949,(sp24:5.527011086,(sp25:5.478875203,(sp40:1.898502308,sp41:1.898502308):3.580372894):0.04813588287):2.489844168):0.5785721075):0.8605508177,((sp47:1.324188282,sp48:1.324188282):1.210143714,sp38:2.534331996):6.921646183):1.848794077,(sp22:6.144323416,sp23:6.144323416):5.160448839):4.752352041,sp10:16.0571243):2.234952832):13.29607763):8.9940146,(sp6:33.80408947,(((sp29:4.271294196,sp30:4.271294196):3.963360008,(sp46:1.515605972,(sp51:0.6842469553,sp52:0.6842469553):0.8313590168):6.719048232):21.69107479,((((sp44:1.517683119,sp45:1.517683119):13.83340518,((sp33:3.451233406,sp34:3.451233406):7.318030201,sp14:10.76926361):4.581824694):2.3268441,((sp31:3.988873926,sp32:3.988873926):13.39833,(sp26:5.46221229,sp27:5.46221229):11.92499164):0.2907284735):12.10203097,((sp16:9.676541191,sp17:9.676541191):11.55054389,(sp11:16.00734921,(sp55:0.6152478573,sp56:0.6152478573):15.39210136):5.219735869):8.552878292):0.1457656227):3.878360468):6.778079891):0.0;"; // test 1
		String treeStr2 = "(((sp11:7.520431298,sp12:7.520431298):14.58559324,(((sp15:7.127781869,((sp35:0.4791031718,sp36:0.4791031718):0.6918706338,sp29:1.170973806):5.956808064):0.9713601967,sp14:8.099142066):13.81179333,(((sp10:11.2968662,sp16:11.2968662):4.233766343,(((sp27:1.270878296,sp28:1.270878296):3.103802529,sp25:4.374680825):1.378878877,(sp41:0.2553128505,sp42:0.2553128505):5.498246851):9.777072843):4.554919652,((sp26:3.972596005,(sp39:0.294787963,sp40:0.294787963):3.677808042):15.28504624,sp1:19.25764224):0.8279099527):1.825383202):0.1950891369):0.2486838349,((sp18:6.128125769,sp19:6.128125769):6.076105996,(sp20:10.59925438,((((sp30:1.085875942,sp31:1.085875942):4.632366157,sp22:5.718242099):3.38586537,(sp24:4.67364838,(sp33:0.5364403626,sp34:0.5364403626):4.137208017):4.43045909):0.07736434735,((sp37:0.4032131745,sp38:0.4032131745):0.3340409803,sp32:0.7372541548):8.444217662):1.417782565):1.604977382):10.15047661):0.0;"; // test 2
        TreeParser myTree1 = new TreeParser(treeStr1, false, false, true, 0); // true b/c species are labelled, offset=0
        TreeParser myTree2 = new TreeParser(treeStr2, false, false, true, 0);
        
        boolean incorporateCladogenesis = false;
                
        HiddenStateDependentSpeciationExtinctionProcess hsdsep1 = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep1.initByName(
        		"tree", myTree1,
        		"hiddenTraitStash", hiddenTraitStash1,
        		"hiddenInstantaneousRateMatrix", hirm1,
        		"lambdaMuAssigner", lambdaMuAssigner1,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
        
        HiddenStateDependentSpeciationExtinctionProcess hsdsep2 = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep2.initByName(
        		"tree", myTree2,
        		"hiddenTraitStash", hiddenTraitStash2,
        		"hiddenInstantaneousRateMatrix", hirm2,
        		"lambdaMuAssigner", lambdaMuAssigner2,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
        System.out.println(hsdsep1.calculateLogP());
        
//    	System.out.println(hsdsep2.calculateLogP());
	}
}
