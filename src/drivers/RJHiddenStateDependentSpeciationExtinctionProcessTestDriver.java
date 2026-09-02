package drivers;

import java.util.Arrays;
import java.util.List;

import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;
import SSE.HiddenTraitStash;
import SSE.LambdaMuAssigner;
import SSE.MasqueradeBall;
import SSE.BSSVSStateDependentSpeciationExtinctionProcess;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;

public class RJHiddenStateDependentSpeciationExtinctionProcessTestDriver {

	public static void main(String[] args) {
		int numberOfStates = 2;
		int numberOfHiddenStates = 2; // this is equivalent to fig. 1 in HiSSE paper (just one hidden state, linked to one of the observed states)
		
		String[] spNames = new String[] { "sp1", "sp10", "sp11", "sp12", "sp14", "sp15", "sp16", "sp18", "sp19", "sp20", "sp22", "sp24", "sp25", "sp26", "sp27", "sp28", "sp29", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		
		String hiddenStatesString = "0,1"; // observed state 1 and 2 will transition to no hidden state, and hidden states 1, respectively
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		
		HiddenTraitStash hiddenTraitStash = new HiddenTraitStash();
		hiddenTraitStash.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet, "hiddenObsStateMapper", stateMapper, "value", "sp1=2,sp10=1,sp11=1,sp12=1,sp14=1,sp15=2,sp16=1,sp18=1,sp19=1,sp20=2,sp22=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp29=2,sp30=1,sp31=1,sp32=2,sp33=1,sp34=1,sp35=2,sp36=2,sp37=2,sp38=2,sp39=2,sp40=2,sp41=2,sp42=2");
		// hiddenTraitStash.printLksMap();
		
		Double[] pis = new Double[] { 0.0, 0.0, 0.0, 0.0, .25, .25, .25, .25 };
		RealParameter pi = new RealParameter(pis);
		
		String lambdasToStatesString = "0,1,2,3";
		Double lambda1 = 0.08885618; // 0A
		Double lambda2 = 0.1509081; // 1A
		Double lambda3 = 0.1; // 0B
		Double lambda4 = 3.597251; // 1B
		Double[] lambdas = { lambda1, lambda2, lambda3, lambda4 };
		RealParameter lambda = new RealParameter(lambdas);
		
		String musToStatesString = "0,1,2,3";
		Double mu1 = 0.0224435; // 0A
		Double mu2 = 0.05379272; // 1A
		Double mu3 = 0.01; // 0B
		Double mu4 = 1.809666E-08; // 1B
		Double[] mus = { mu1, mu2, mu3, mu4 };
		RealParameter mu = new RealParameter(mus);
		
		LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
		lambdaMuAssigner.initByName("totalNumberOfStates", 4, 
				"nDistinctLambdas", 4, 
				"nDistinctMus", 4, 
				"lambdasToStates", lambdasToStatesString, 
				"lambda", lambda, 
				"musToStates", musToStatesString, 
				"mu", mu, 
				"pi", pi);
		
		boolean disallowDoubleTransitions = true;
		int symmetrifyAcrossDiagonal = -1;
		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString = "0.02597671 2.061154E-09 0.04404249 2.061154E-09 2.355085E-09 1.565876E-02 4.59653214 100";
		
		hirm.initByName("numberOfStates", 2, 
				"numberOfHiddenStates", 2, 
				"flatQMatrix", flatQMatrixString, 
				"disallowDoubleTransitions", disallowDoubleTransitions, 
				"hiddenObsStateMapper", stateMapper, 
				"symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal); 
		System.out.println("\nFull Q Matrix");
		hirm.printMatrix();
		System.out.println();
		
		// Up to this point, we have initialized everything to two hidden and two observed states
		// Below, we apply a mask and force just one hidden state -- it should match HSDSEPHiSSETest2
		
		Double[] mask1 = new Double[] { 0.0, 2.0, 0.0 }; // last one means CID is off
		RealParameter modelMask1 = new RealParameter(mask1);
		
		MasqueradeBall masqueradeBall = new MasqueradeBall();
		masqueradeBall.initByName("modelMask", modelMask1, 
				"hiddenTraitStash", hiddenTraitStash,
				"hiddenInstantaneousRateMatrix", hirm, 
				"lambdaMuAssigner", lambdaMuAssigner);
		
		String treeStr = "(((sp11:7.520431298,sp12:7.520431298):14.58559324,(((sp15:7.127781869,((sp35:0.4791031718,sp36:0.4791031718):0.6918706338,sp29:1.170973806):5.956808064):0.9713601967,sp14:8.099142066):13.81179333,(((sp10:11.2968662,sp16:11.2968662):4.233766343,(((sp27:1.270878296,sp28:1.270878296):3.103802529,sp25:4.374680825):1.378878877,(sp41:0.2553128505,sp42:0.2553128505):5.498246851):9.777072843):4.554919652,((sp26:3.972596005,(sp39:0.294787963,sp40:0.294787963):3.677808042):15.28504624,sp1:19.25764224):0.8279099527):1.825383202):0.1950891369):0.2486838349,((sp18:6.128125769,sp19:6.128125769):6.076105996,(sp20:10.59925438,((((sp30:1.085875942,sp31:1.085875942):4.632366157,sp22:5.718242099):3.38586537,(sp24:4.67364838,(sp33:0.5364403626,sp34:0.5364403626):4.137208017):4.43045909):0.07736434735,((sp37:0.4032131745,sp38:0.4032131745):0.3340409803,sp32:0.7372541548):8.444217662):1.417782565):1.604977382):10.15047661):0.0;";
		TreeParser myTree = new TreeParser(treeStr, false, false, true, 0);
		
		boolean incorporateCladogenesis = false;
		
		BSSVSStateDependentSpeciationExtinctionProcess rjhsdsep = new BSSVSStateDependentSpeciationExtinctionProcess();
        rjhsdsep.initByName(
        		"tree", myTree,
        		"hiddenTraitStash", hiddenTraitStash,
        		"masqueradeBall", masqueradeBall, 
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
        
        System.out.println(rjhsdsep.calculateLogP()); // -100.67039644731487
	}
}
