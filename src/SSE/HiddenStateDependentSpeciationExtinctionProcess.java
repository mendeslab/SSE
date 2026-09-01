package SSE;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.Executors;
import beast.base.core.ProgramStatus;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.*;

@Description("More general class capable of inferring under ClaSSE model, but also (Mu)HiSSE.")
@Citation(value="Goldberg EE, Igic B (2012) Tempo and mode in plant breeding system evolution. Evolution 16(12):3701-3709",
year=2012, firstAuthorSurname="Goldberg", DOI="10.1111/j.1558-5646.2012.01730.x")
public class HiddenStateDependentSpeciationExtinctionProcess extends StateDependentSpeciationExtinctionProcess {

	final public Input<HiddenTraitStash> hiddenTraitStashInput = new Input<>("hiddenTraitStash", "TraitStash object containing the observed character state for each species.", Validate.REQUIRED);
	final public Input<HiddenInstantaneousRateMatrix> hirmInput = new Input<>("hiddenInstantaneousRateMatrix", "HiddenInstantaneousRateMatrix object containing anagenenetic rates for both observed and hidden states.", Validate.REQUIRED);
	final public Input<LambdaMuAssigner> lambdaMuAssignerInput = new Input<>("lambdaMuAssigner", "LambdaMuAssigner object that assigns distinct parameters to each state.", Validate.REQUIRED);

	// input
	protected HiddenTraitStash hiddenTraitStash;
	protected HiddenInstantaneousRateMatrix q;
	private LambdaMuAssigner lambdaMuAssigner;
	protected int numHiddenStates;
	protected int numStates; // obs + hidden
    
    public HiddenStateDependentSpeciationExtinctionProcess() {
    	traitStashInput.setRule(Input.Validate.FORBIDDEN);
    	irmInput.setRule(Input.Validate.FORBIDDEN);
    	lambdaInput.setRule(Input.Validate.FORBIDDEN);
    	muInput.setRule(Input.Validate.FORBIDDEN);
    	piInput.setRule(Input.Validate.FORBIDDEN);
    	cladoStashInput.setRule(Input.Validate.OPTIONAL);
	}
    
	@Override
	public void initAndValidate() {		
		tree = treeInput.get();
		hiddenTraitStash = hiddenTraitStashInput.get();
		q = hirmInput.get();
		lambdaMuAssigner = lambdaMuAssignerInput.get();
		incorporateCladogenesis = cladoFlagInput.get();
		numObsStates = q.getNumObsStates();
		numHiddenStates = q.getNumHiddenStates();
		numStates = numObsStates + numHiddenStates;
		rate = 1.0;
		
		if (incorporateCladogenesis) {
			cladoStash = lambdaMuAssigner.getCladoStash();
		}
		else { 
			lambda = lambdaMuAssigner.getLambdas();
		}
		
		pi = lambdaMuAssigner.getPis();
		mu = lambdaMuAssigner.getMus();
		
		// likelihood-related
		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
		scalingConstants = new double[tree.getNodeCount()]; // equivalent to diversitree's lq (but not in log-scale), these are used as denominators during the likelihood computation
		
		Arrays.fill(scalingConstants, 1.0);		
		finalLk = 0.0;
		finalLogLk = 0.0;
		
		// cache-related
		storedNodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2];  
		storedScalingConstants = new double[tree.getNodeCount()]; 
		branchLengths = new double[tree.getNodeCount()];
		storedBranchLengths = new double[tree.getNodeCount()];
		hasDirt = Tree.IS_FILTHY;

        useThreads = useThreadsInput.get() && (ProgramStatus.m_nThreads > 1);
		nrOfThreads = useThreads ? ProgramStatus.m_nThreads : 1;
		if (useThreads && maxNrOfThreadsInput.get() > 0) {
			nrOfThreads = Math.min(maxNrOfThreadsInput.get(), ProgramStatus.m_nThreads);
		}
		if (useThreads) {
		     exec = Executors.newFixedThreadPool(nrOfThreads);
		}

		startStates = new int[tree.getNodeCount()];
		endStates = new int[tree.getNodeCount()];
		nodeConditionalScaledLks = new double[tree.getNodeCount()][numStates*2];

		sampleCharacterHistory = false;
		branchPartialLks = new ArrayList[tree.getNodeCount()];
//        dt = tree.getRoot().getHeight() / numTimeSlices * 50.0; // why we multiply by 50? following RevBayes code
		dt = tree.getRoot().getHeight() / numTimeSlices; // why we multiply by 50? following RevBayes code
		nodeTransitionStates = new ArrayList[tree.getNodeCount()];
		nodeTransitionTimes = new ArrayList[tree.getNodeCount()];
		nodeTimeInState = new double[tree.getNodeCount()][numStates];
		numNodeStateChanges = 0;
		// numBranchStateChanges = 0;
		averageSpeciationRates = new double[tree.getNodeCount()];
		averageExtinctionRates = new double[tree.getNodeCount()];
	}
	
	@Override
	public double calculateLogP() {
		mu = lambdaMuAssigner.getMus();
		pi = lambdaMuAssigner.getPis();

		if (!incorporateCladogenesis) {
			lambda = lambdaMuAssigner.getLambdas();
		}

		// when a biogeographical parameter changes, tree is filthy, we can use threads
		// otherwise, caching allows us to only recompute part of the tree likelihood, and thread overhead not worth it
		if (hasDirt == Tree.IS_FILTHY && useThreads) {
			computeNodeLkUsingThreads();
		} else {
			computeNodeLk(tree.getRoot(), true);
		}
		logP = finalLogLk;
		return logP;
	}
	
	// possibly hidden states, so total number of states = number of obs + hidden states
	@Override
	protected int getTotalNumberStates() {
		return numStates;
	}
	
	@Override
	protected void initializeLeafLks(Node aNode, double[] aNodePartial) {
		System.arraycopy(hiddenTraitStash.getSpLks(aNode.getID()), 0, aNodePartial, 0, aNodePartial.length);
	}
	
	@Override
	protected void numericallyIntegrateProcess(double[] likelihoods, double beginAge, double endAge, boolean backwardTime, boolean extinctionOnly) {
		SSEODE ode = new HiddenSSEODE(mu, q, rate, incorporateCladogenesis, backwardTime, extinctionOnly);
		solveODE(likelihoods, beginAge, endAge, ode);
	}
	
    @Override
	protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;
        
		if ((hirmInput.get().isDirtyCalculation()) ||
			(lambdaMuAssigner.getLambdasRealParameter() != null && lambdaMuAssigner.getLambdasRealParameter().somethingIsDirty()) || 
			(lambdaMuAssigner.getMusRealParameter().somethingIsDirty()) ||
			(lambdaMuAssigner.getPisRealParameter().somethingIsDirty()) ||
			(lambdaMuAssigner.getCladoStash() != null && lambdaMuAssigner.getCladoStash().isDirtyCalculation())
			) 
		{
			hasDirt = Tree.IS_FILTHY;
			return true;
		}
        return treeInput.get().somethingIsDirty();
	}
    
}
