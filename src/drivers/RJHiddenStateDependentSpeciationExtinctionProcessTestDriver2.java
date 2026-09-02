package drivers;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;
import SSE.HiddenStateDependentSpeciationExtinctionProcess;
import SSE.HiddenTraitStash;
import SSE.LambdaMuAssigner;
import SSE.MasqueradeBall;
import SSE.BSSVSStateDependentSpeciationExtinctionProcess;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;

public class RJHiddenStateDependentSpeciationExtinctionProcessTestDriver2 {

	public static void main(String[] args) {
		String[] spNames = new String[] { "sp4", "sp6", "sp10", "sp11", "sp12", "sp15", "sp16", "sp17", "sp18", "sp20", "sp21", "sp23", "sp24", "sp25", "sp26", "sp27", "sp28", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		
		String hiddenStatesString1 = "0,1,2,3"; // observed state 1,2,3 and 4 will transition to hidden states 1,2,3 and 4, respectively
		HiddenObservedStateMapper stateMapper1 = new HiddenObservedStateMapper();
		stateMapper1.initByName("hiddenStates", hiddenStatesString1);
		
		HiddenTraitStash hiddenTraitStash1 = new HiddenTraitStash();
		hiddenTraitStash1.initByName("numberOfStates", 4, "numberOfHiddenStates", 4, "taxa", taxonSet, "hiddenObsStateMapper", stateMapper1, "value", "sp4=1,sp6=4,sp10=2,sp11=3,sp12=1,sp15=1,sp16=4,sp17=1,sp18=4,sp20=2,sp21=3,sp23=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=4,sp36=2,sp37=2,sp38=2,sp39=1,sp40=1,sp41=2,sp42=2");
		// hiddenTraitStash.printLksMap();
		
		String hiddenStatesString2 = "-1,0,-1,1"; // observed state 1,2,3 and 4 will transition to hidden states 1,2,3 and 4, respectively
		HiddenObservedStateMapper stateMapper2 = new HiddenObservedStateMapper();
		stateMapper2.initByName("hiddenStates", hiddenStatesString2);
		
		HiddenTraitStash hiddenTraitStash2 = new HiddenTraitStash();
		hiddenTraitStash2.initByName("numberOfStates", 4, "numberOfHiddenStates", 2, "taxa", taxonSet, "hiddenObsStateMapper", stateMapper2, "value", "sp4=1,sp6=4,sp10=2,sp11=3,sp12=1,sp15=1,sp16=4,sp17=1,sp18=4,sp20=2,sp21=3,sp23=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=4,sp36=2,sp37=2,sp38=2,sp39=1,sp40=1,sp41=2,sp42=2");
		// hiddenTraitStash.printLksMap();
		
		// note that here there are two hidden traits, each with 4 states, but the lambdas and mus are shared across traits
		// also, transitions only occur across states of the same hidden trait, or within the same state across traits
		// i.e., this test uses the HiSSE machinery, but boils down to MuSSE
		String lambdasToStatesString1 = "0,1,2,3,4,5,6,7";
		String lambdasToStatesString2 = "0,1,2,3,4,5";
		Double lambda1 = 0.1; // 0A
		Double lambda2 = 0.15; // 1A
		Double lambda3 = 0.2; // 2A
		Double lambda4 = 0.1; // 3A
		Double lambda5 = 0.1; // 0B
		Double lambda6 = 0.15; // 1B
		Double lambda7 = 0.2; // 2B
		Double lambda8 = 0.1; // 3B
		Double[] lambdas1 = { lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8 };
		Double[] lambdas2 = { lambda1, lambda2, lambda3, lambda4, lambda6, lambda8 };
		RealParameter lambdaOne = new RealParameter(lambdas1);
		RealParameter lambdaTwo = new RealParameter(lambdas2);
				
		String musToStatesString1 = "0,1,2,3,4,5,6,7";
		String musToStatesString2 = "0,1,2,3,4,5";
		Double mu1 = 0.03;
		Double mu2 = 0.045;
		Double mu3 = 0.06;
		Double mu4 = 0.03;
		Double mu5 = 0.03;
		Double mu6 = 0.045;
		Double mu7 = 0.06;
		Double mu8 = 0.03;
		Double[] mus1 = { mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8 };
		Double[] mus2 = { mu1, mu2, mu3, mu4, mu6, mu8 };
		RealParameter muOne = new RealParameter(mus1);
		RealParameter muTwo = new RealParameter(mus2);
		
		Double[] pis1 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0 };
		RealParameter pi1 = new RealParameter(pis1);
		
		Double[] pis2 = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0};
		RealParameter pi2 = new RealParameter(pis2);

		LambdaMuAssigner lambdaMuAssigner1 = new LambdaMuAssigner();
		lambdaMuAssigner1.initByName("totalNumberOfStates", 8, 
				"nDistinctLambdas", 8, 
				"nDistinctMus", 8, 
				"lambdasToStates", lambdasToStatesString1, 
				"lambda", lambdaOne, 
				"musToStates", musToStatesString1, 
				"mu", muOne, 
				"pi", pi1);

		LambdaMuAssigner lambdaMuAssigner2 = new LambdaMuAssigner();
		lambdaMuAssigner2.initByName("totalNumberOfStates", 6, 
				"nDistinctLambdas", 6, 
				"nDistinctMus", 6, 
				"lambdasToStates", lambdasToStatesString2, 
				"lambda", lambdaTwo, 
				"musToStates", musToStatesString2, 
				"mu", muTwo, 
				"pi", pi2);
		
		boolean disallowDoubleTransitions = true; // not used
		int symmetrifyAcrossDiagonal = -1;
		
		HiddenInstantaneousRateMatrix hirm1 = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString1 = "0.05 0.05 0.0 1.0 0.05 0.0 0.05 1.0 0.05 0.0 0.05 1.0 0.0 0.05 0.05 1.0 1.0 0.05 0.05 0.0 1.0 0.05 0.0 0.05 1.0 0.05 0.0 0.05 1.0 0.0 0.05 0.05";
		hirm1.initByName("numberOfStates", 4, "numberOfHiddenStates", 4, "flatQMatrix", flatQMatrixString1, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper1);	
		// System.out.println("hirm1 Qs:");
		// hirm1.printMatrix();
		
		HiddenInstantaneousRateMatrix hirm2 = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString2 = "0.05 0.05 0.0 0.05 0.0 0.05 1.0 0.05 0.0 0.05 0.0 0.05 0.05 1.0 1.0 0.05 1.0 0.05";
		// String flatQMatrixString2 = "0.05 0.05 0.05 0.0 0.0 0.05 0.05 0.05 0.05 0.05 0.05 0.05 0.0 0.0 0.05 0.05 0.05 0.05 0.0 0.05 0.0 0.0 0.0 0.0 0.05 0.0";
		hirm2.initByName("numberOfStates", 4, "numberOfHiddenStates", 2, "flatQMatrix", flatQMatrixString2, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper2);
		// System.out.println("hirm2 Qs:");
		// hirm2.printMatrix();
				
//		Double[] mask1 = new Double[] { 2.0, 2.0, 2.0, 2.0, 0.0 }; // last one means CID is off
//		RealParameter modelMask1 = new RealParameter(mask1);
		
		Integer[] aStatesMaskPart = new Integer[] { 2, 2, 2, 2 };
		IntegerParameter stateMask1 = new IntegerParameter(aStatesMaskPart);
		Integer[] aCIDMaskPart = new Integer[] { 0 };
		IntegerParameter cidMask1 = new IntegerParameter(aCIDMaskPart);
		
		MasqueradeBall masqueradeBall = new MasqueradeBall();
		masqueradeBall.initByName(
				"stateMask", stateMask1,
				"cidMask", cidMask1,
				"hiddenTraitStash", hiddenTraitStash1,
				"hiddenInstantaneousRateMatrix", hirm1, 
				"lambdaMuAssigner", lambdaMuAssigner1);
		
		String treeStr = "(((sp39:0.518912972,sp40:0.518912972):19.54206195,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
        
        HiddenStateDependentSpeciationExtinctionProcess hsdsep1 = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep1.initByName(
        		"tree", myTree,
        		"hiddenTraitStash", hiddenTraitStash1,
        		"lambdaMuAssigner", lambdaMuAssigner1,
        		"hiddenInstantaneousRateMatrix", hirm1,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
        
        HiddenStateDependentSpeciationExtinctionProcess hsdsep2 = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep2.initByName(
        		"tree", myTree,
        		"hiddenTraitStash", hiddenTraitStash2,
        		"lambdaMuAssigner", lambdaMuAssigner2,
        		"hiddenInstantaneousRateMatrix", hirm2,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
        
        BSSVSStateDependentSpeciationExtinctionProcess rjhsdsep = new BSSVSStateDependentSpeciationExtinctionProcess();
        rjhsdsep.initByName(
        		"tree", myTree,
        		"hiddenTraitStash", hiddenTraitStash1,
        		"masqueradeBall", masqueradeBall,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	double negLnl1 = hsdsep1.calculateLogP();
    	System.out.println("hsdsep1: " + negLnl1); // -122.88014179920914
    	
    	double negLnl2 = hsdsep2.calculateLogP();
    	System.out.println("hsdsep2: " + negLnl2); // -122.84037671593603
    	
    	double negLnl3 = rjhsdsep.calculateLogP();
    	System.out.println("rjhsdsep1: " + negLnl3); // -122.88014179920914
    	
//    	Double[] mask2 = new Double[] { 0.0, 2.0, 0.0, 2.0, 0.0 }; // applying new mask
    	Integer[] aStatesMaskPart2 = new Integer[] { 0, 2, 0, 2 };
		Integer[] aCIDMaskPart2 = new Integer[] { 0 };
    	rjhsdsep.setMask(aStatesMaskPart2, aCIDMaskPart2);
		double negLnl4 = rjhsdsep.calculateLogP();
    	System.out.println("rjhsdsep2: " + negLnl4); // -122.84037671593603
	}
}
