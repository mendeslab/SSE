package SSE;

import java.util.Arrays;
import java.util.concurrent.Executors;

import beast.base.core.ProgramStatus;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Tree;

public class BSSVSStateDependentSpeciationExtinctionProcess extends HiddenStateDependentSpeciationExtinctionProcess {

	final public Input<MasqueradeBall> masqueradeBallInput = new Input<>("masqueradeBall", "MaskeradeBall object that allows resetting lambdas, mus and Q's depending on the model mask.", Validate.REQUIRED);
	
	private MasqueradeBall masqueradeBall;
	
	public BSSVSStateDependentSpeciationExtinctionProcess() {
		hirmInput.setRule(Input.Validate.FORBIDDEN);
		lambdaMuAssignerInput.setRule(Input.Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
		masqueradeBall = masqueradeBallInput.get();
		tree = treeInput.get();
		hiddenTraitStash = hiddenTraitStashInput.get();
		incorporateCladogenesis = cladoFlagInput.get();
		
		if (incorporateCladogenesis) {
			cladoStash = masqueradeBall.getCladoStash();
		}
		else { 
			lambda = masqueradeBall.getLambdas();
		}
		
		pi = masqueradeBall.getPis();
		mu = masqueradeBall.getMus();
		q = masqueradeBall.getHIRM();
		
		numObsStates = q.getNumObsStates();
		numHiddenStates = q.getNumHiddenStates();
		numStates = numObsStates + numHiddenStates;
		rate = 1.0;
		
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
	}
	
	@Override
	public double calculateLogP() {
		pi = masqueradeBall.getPis();
		mu = masqueradeBall.getMus();
		q = masqueradeBall.getHIRM();

		// no support for cladogenetic changes yet
		if (!incorporateCladogenesis) {
			lambda = masqueradeBall.getLambdas();
		}
		
		// making sure we updated total number of states and of hidden states
		numStates = masqueradeBall.getNumStates(); // obs + hidden
		numHiddenStates = masqueradeBall.getNumHiddenStatesInMask();
		
  		// need to resize partial likelihood vectors
		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
		storedNodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2];  
		
//		System.out.println("Pis in rjhdsep: " + Arrays.toString(pi));
//		System.out.println("Lambdas in rjhdsep: " + Arrays.toString(lambda));
//		System.out.println("Mus in rjhdsep: " + Arrays.toString(mu));
//		System.out.println("Q in rjhdsep");
//		q.printMatrix();
		
		// when a character state parameter changes, tree is filthy, we can use threads
		// otherwise, caching allows us to only recompute part of the tree likelihood, and thread overhead not worth it
		if (hasDirt == Tree.IS_FILTHY && useThreads) {
			computeNodeLkUsingThreads();
		} else {
			computeNodeLk(tree.getRoot(), true);
		}
		
		logP = finalLogLk;
		return logP;
	}
	
	// used for testing
	public void setMask(Integer[] aStatesMaskPart, Integer[] aCIDMaskPart) {
		masqueradeBall.setMask(aStatesMaskPart, aCIDMaskPart);
	}
	
	@Override
	protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;
        
        if ( (masqueradeBallInput.get().isDirtyCalculation()) ||
        	 (masqueradeBall.getQRealParameter() != null && masqueradeBall.getQRealParameter().somethingIsDirty()) ||
        	 (masqueradeBall.getLambdasRealParameter() != null && masqueradeBall.getLambdasRealParameter().somethingIsDirty()) ||
        	 (masqueradeBall.getMusRealParameter() != null && masqueradeBall.getMusRealParameter().somethingIsDirty()) ||
        	 (masqueradeBall.getPisRealParameter().somethingIsDirty())
           ) 
		{
			hasDirt = Tree.IS_FILTHY;
			return true;
		}
        return treeInput.get().somethingIsDirty();
	}
}
