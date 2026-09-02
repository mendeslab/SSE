package test;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.CladoTriplet;
import SSE.CladogeneticSpeciationRateStash;
import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;
import SSE.HiddenStateDependentSpeciationExtinctionProcess;
import SSE.HiddenTraitStash;
import SSE.LambdaMuAssigner;
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

public class HSDSEPClaSSETest {
	final static double EPSILON1 = 1e-4;
	final static double EPSILON2 = 1e-10;
	private double negLnl;

	@Before
	public void setUp() throws Exception {
		int numberOfStates = 4; // 1:Null states, 2:A, 3:B, 4:A+B
		int numberOfHiddenStates = 0;
		int totalNumberOfStates = numberOfStates + numberOfHiddenStates;
		String[] spNames = new String[] { "sp4", "sp6", "sp10", "sp11", "sp12", "sp15", "sp16", "sp17", "sp18", "sp20", "sp21", "sp23", "sp24", "sp25", "sp26", "sp27", "sp28", "sp30", "sp31", "sp32", "sp33", "sp34", "sp35", "sp36", "sp37", "sp38", "sp39", "sp40", "sp41", "sp42" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		
		HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper(); // not used
		
		HiddenTraitStash hiddenTraitStash = new HiddenTraitStash();
		hiddenTraitStash.initByName("numberOfStates", numberOfStates, "numberOfHiddenStates", numberOfHiddenStates, "taxa", taxonSet, "hiddenObsStateMapper", stateMapper, "value", "sp4=1,sp6=4,sp10=2,sp11=3,sp12=1,sp15=1,sp16=4,sp17=1,sp18=4,sp20=2,sp21=3,sp23=1,sp24=1,sp25=2,sp26=2,sp27=2,sp28=2,sp30=2,sp31=2,sp32=2,sp33=2,sp34=2,sp35=4,sp36=2,sp37=2,sp38=2,sp39=1,sp40=1,sp41=2,sp42=2");
		hiddenTraitStash.printLksMap();
		
		double sympProb = 1.0; // DEC-like
		double subsympProb = 1.0 / 6.0;
		double vicProb = 1.0 / 6.0;
		double jProb = 0.0; // no jump dispersal
		double birthRate = 0.2;
		Double[] sSpeciationRate = { sympProb * birthRate };
		Double[] ssSpeciationRate = { subsympProb * birthRate };
		Double[] vSpeciationRate = { vicProb * birthRate };
		Double[] jSpeciationRate = { jProb * birthRate }; // 0.0
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
			
		String musToStatesString = "0,0,0,0";
		Double mu1 = 0.1;
		Double[] mus = { mu1 };
		RealParameter mu = new RealParameter(mus);
		
		Double[] piEs = new Double[totalNumberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[totalNumberOfStates];
		Arrays.fill(piDs, (1.0/totalNumberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25
		RealParameter pi = new RealParameter(pis);
		
		LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
		lambdaMuAssigner.initByName("totalNumberOfStates", 4, "nDistinctMus", 1, "musToStates", musToStatesString, "mu", mu, "pi", pi, "cladogeneticStash", csrt);
	
		System.out.println("Lambdas: " + Arrays.toString(lambdaMuAssigner.getLambdas()));
		System.out.println("Mus: " + Arrays.toString(lambdaMuAssigner.getMus()));
		
		boolean disallowDoubleTransitions = false; // not used
		int symmetrifyAcrossDiagonal = -1; // not used
		HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
		String flatQMatrixString = "0.0 0.0 0.0 0.01 0.0 0.01 0.01 0.0 0.01 0.0 0.01 0.01";
		hirm.initByName("numberOfStates", 4, "numberOfHiddenStates", 0, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper); // MuSSE
		
		System.out.println("Qs:");
		hirm.printMatrix();
		
		String treeStr = "(((sp39:0.518912972,sp40:0.518912972):19.54206195,((((sp25:3.198513788,(sp32:2.293402763,(sp41:0.1728996412,sp42:0.1728996412):2.120503122):0.9051110254):7.525323533,((sp20:9.577599427,(sp26:2.751623892,(sp30:2.405609293,sp31:2.405609293):0.3460145989):6.825975535):0.05909341512,(sp17:7.221384607,sp18:7.221384607):2.415308235):1.087144479):0.5715875464,sp10:11.29542487):6.453137462,((sp33:2.252609903,sp34:2.252609903):4.989398146,sp16:7.24200805):10.50655428):2.312412597):0.3952852439,(((((sp23:5.605262355,sp24:5.605262355):3.179619681,((sp37:1.329072526,sp38:1.329072526):1.780228265,(sp27:2.543803164,sp28:2.543803164)nd38:0.5654976265):5.675581245):6.165501477,sp6:14.95038351):0.6290423683,((sp11:8.298747349,(sp15:8.099808068,sp12:8.099808068):0.1989392817):3.788483262,(sp21:6.890801228,(sp35:1.989124199,sp36:1.989124199):4.901677029):5.196429383):3.492195269):2.878773551,sp4:18.45819943):1.998060738)0.0;";
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
		
        boolean incorporateCladogenesis = true;
        
        HiddenStateDependentSpeciationExtinctionProcess hsdsep = new HiddenStateDependentSpeciationExtinctionProcess();
        hsdsep.initByName(
        		"tree", myTree,
        		"hiddenTraitStash", hiddenTraitStash,
        		"lambdaMuAssigner", lambdaMuAssigner,
        		"hiddenInstantaneousRateMatrix", hirm,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	negLnl = hsdsep.calculateLogP();
    	System.out.println(negLnl); // 
        
	}

	@Test
	public void againstDiversitreeClaSSE() {
		Assert.assertEquals(-129.9762, negLnl, EPSILON1); 
	}

//	@Test
//	public void againstMyClaSSE() {
//		Assert.assertEquals(-129.97617749574326, negLnl, EPSILON2);
//	}


}
