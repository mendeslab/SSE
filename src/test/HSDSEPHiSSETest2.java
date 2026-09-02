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

public class HSDSEPHiSSETest2 {
	final static double EPSILON = 1E-3;
	private double negLnl;

	@Before
	public void setUp() throws Exception {
		int numberOfStates = 2;
		int numberOfHiddenStates = 1; // this is equivalent to fig. 1 in HiSSE paper (just one hidden state, linked to one of the observed states)
		int totalNumberOfStates = numberOfStates + numberOfHiddenStates;
		String[] spNames = new String[] { "sp1", "sp10", "sp11", "sp12", "sp14", "sp15", "sp16", "sp18", "sp19", "sp20", "sp22", "sp24", "sp25", "sp26", "sp27", "sp28", "sp29", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		
		String hiddenStatesString = "-1,0"; // observed state 2 will transition to hidden states 1
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
		stateMapper.initByName("hiddenStates", hiddenStatesString);
		
		HiddenTraitStash hiddenTraitStash = new HiddenTraitStash();
		hiddenTraitStash.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet, "hiddenObsStateMapper", stateMapper, "value", "sp1=2,sp10=1,sp11=1,sp12=1,sp14=1,sp15=2,sp16=1,sp18=1,sp19=1,sp20=2,sp22=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp29=2,sp30=1,sp31=1,sp32=2,sp33=1,sp34=1,sp35=2,sp36=2,sp37=2,sp38=2,sp39=2,sp40=2,sp41=2,sp42=2");
		hiddenTraitStash.printLksMap();
		
		String lambdasToStatesString = "0,1,2";
		Double lambda1 = 0.08885618; // 0A
		Double lambda2 = 0.1509081; // 0B
		Double lambda3 = 3.597251; // 1B
		Double[] lambdas = { lambda1, lambda2, lambda3 };
		RealParameter lambda = new RealParameter(lambdas);
		
		String musToStatesString = "0,1,2";
		Double mu1 = 0.0224435; // 0A
		Double mu2 = 0.05379272; // 0B
		Double mu3 = 1.809666E-08; // 1B
		Double[] mus = { mu1, mu2, mu3 };		
		RealParameter mu = new RealParameter(mus);

		Double[] piEs = new Double[totalNumberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[totalNumberOfStates];
		Arrays.fill(piDs, (1.0/totalNumberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.33333, 0.33333, 0.33333
		RealParameter pi = new RealParameter(pis);
		
		LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
		lambdaMuAssigner.initByName("totalNumberOfStates", 3, "nDistinctLambdas", 3, "nDistinctMus", 3, "lambdasToStates", lambdasToStatesString, "lambda", lambda, "musToStates", musToStatesString, "mu", mu, "pi", pi);
		System.out.println("Lambdas: " + Arrays.toString(lambdaMuAssigner.getLambdas()));
		System.out.println("Mus: " + Arrays.toString(lambdaMuAssigner.getMus()));
		System.out.println("Pis: " + Arrays.toString(lambdaMuAssigner.getPis()));
		
		boolean disallowDoubleTransitions = true;
		int symmetrifyAcrossDiagonal = -1;
		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString = "0.02597671 0.04404249 2.061154E-09 4.59653214"; // test 1 
		hirm.initByName("numberOfStates", 2, "numberOfHiddenStates", 1, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper);
		
		System.out.println("Qs:");
		hirm.printMatrix();
		
		String treeStr = "(((sp11:7.520431298,sp12:7.520431298):14.58559324,(((sp15:7.127781869,((sp35:0.4791031718,sp36:0.4791031718):0.6918706338,sp29:1.170973806):5.956808064):0.9713601967,sp14:8.099142066):13.81179333,(((sp10:11.2968662,sp16:11.2968662):4.233766343,(((sp27:1.270878296,sp28:1.270878296):3.103802529,sp25:4.374680825):1.378878877,(sp41:0.2553128505,sp42:0.2553128505):5.498246851):9.777072843):4.554919652,((sp26:3.972596005,(sp39:0.294787963,sp40:0.294787963):3.677808042):15.28504624,sp1:19.25764224):0.8279099527):1.825383202):0.1950891369):0.2486838349,((sp18:6.128125769,sp19:6.128125769):6.076105996,(sp20:10.59925438,((((sp30:1.085875942,sp31:1.085875942):4.632366157,sp22:5.718242099):3.38586537,(sp24:4.67364838,(sp33:0.5364403626,sp34:0.5364403626):4.137208017):4.43045909):0.07736434735,((sp37:0.4032131745,sp38:0.4032131745):0.3340409803,sp32:0.7372541548):8.444217662):1.417782565):1.604977382):10.15047661):0.0;";
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
    	System.out.println(negLnl); // -100.67039644731487

	}

	@Test
	public void againstDiversitreeHiSSE() {
		Assert.assertEquals(-100.67, negLnl, EPSILON); // difference due to precision and rounding
	}
}
