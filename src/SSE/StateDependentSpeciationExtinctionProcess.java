package SSE;

import java.util.*; // TODO Change this to all of them listed, * can add random things that confuse
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;
import java.util.stream.Stream;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import beast.base.core.ProgramStatus;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.tree.*;

@Description("Cladogenetic State change Speciation and Extinction (ClaSSE) model")
@Citation(value="Goldberg EE, Igic B (2012) Tempo and mode in plant breeding system evolution. Evolution 16(12):3701-3709",
year=2012, firstAuthorSurname="Goldberg", DOI="10.1111/j.1558-5646.2012.01730.x")
public class StateDependentSpeciationExtinctionProcess extends Distribution {

	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<TraitStash> traitStashInput = new Input<>("traitStash", "TraitStash object containing the observed character state for each species.", Validate.REQUIRED);
	final public Input<InstantaneousRateMatrix> irmInput = new Input<>("instantaneousRateMatrix", "InstantaneousRateMatrix object containing anagenenetic rates.", Validate.REQUIRED);
	final public Input<CladogeneticSpeciationRateStash> cladoStashInput = new Input<>("cladogeneticStash", "CladogeneticSpeciationRateStash object that generates event map.");
	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each state (if cladogenetic events are not considered).", Validate.XOR, cladoStashInput);
	final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each state.", Validate.REQUIRED);
	final public Input<RealParameter> piInput = new Input<>("pi", "Equilibrium frequencies at root.", Validate.REQUIRED);
	final public Input<Boolean> cladoFlagInput = new Input<>("incorporateCladogenesis", "Whether or not to incorporate cladogenetic events.", Validate.REQUIRED);

	final public Input<Boolean> useThreadsInput = new Input<>("useThreads", "calculated the distributions in parallel using threads (default false)", false);
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

	// input
	protected Tree tree;
	private TraitStash traitStash;
	private InstantaneousRateMatrix q;
	protected CladogeneticSpeciationRateStash cladoStash;
	protected Double[] lambda;
	protected Double[] mu;
	protected Double[] pi; // root eq freqs
	protected int numObsStates;
	protected double rate;
	protected boolean incorporateCladogenesis;
	
	/* Original version: sliced branches into chunks, aimed at fixed step size ODE solvers */ 
	// private double dt; // time slice size (ctor populates)
	// private double rootAge;
	// private int numTimeSlices;
	
	// members used for lk computation
	protected double integratorMinStep = 1.0e-8; // for ODE
	protected double integratorTolerance = 1.0e-6; // for ODE
	protected double[][] nodePartialScaledLksPostOde;
	protected double[] scalingConstants;
	protected double finalLogLk;
	protected double finalLk;
	protected int rootIdx;
	
	// cache used for lk computation
	protected double[][] storedNodePartialScaledLksPostOde;
	protected double[] storedScalingConstants;
	
	protected boolean useThreads;
	protected int nrOfThreads;
	
	/**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] branchLengths;
    protected double[] storedBranchLengths;
	
    /**
     * Flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /* 
     * Stochastic character mapping variables below
     */
    // sampled ancestral states
	public int[] startStates;
	public int[] endStates;

	protected double[][] nodeConditionalScaledLks;
	// (exclusive to forward pass, is initialized to 0's and 1's depending on
	// sampling results, and then acts as the "pi" every time we sample 

	// the paremeters below are used for drawStochasticCharacterMapping and mostly have to do with the activity on branches
	protected boolean sampleCharacterHistory; // RB convention. This is used during drawStochastic, where we sample states along a branch (not just at internal nodes). 
	protected ArrayList<double[]>[] branchPartialLks; // likelihoods along the branch above a node. Arraylist is used because number of samples varies depending on branch length
	protected int numTimeSlices = 500;
	protected double dt; // dt is calculated as a function of the number of time slices
	public ArrayList[] nodeTransitionStates; // tracks transitions along the branch above the node
	public ArrayList[] nodeTransitionTimes; // tracks how long the state was occupied before transitioning
	public double[][] nodeTimeInState; // total time of the node in the state
	// public int numBranchStateChanges; // TODO remove
	public int numNodeStateChanges; // TODO remove
	protected double[] averageSpeciationRates; // over all states in the branch, the speciation rate
	protected double[] averageExtinctionRates; // over all states in the branch, the extinction rate

	@Override
	public void initAndValidate() {		
		super.initAndValidate();
		tree = treeInput.get();
		traitStash = traitStashInput.get();
		q = irmInput.get();
		mu = muInput.get().getValues();
		pi = piInput.get().getValues();
		numObsStates = q.getNumObsStates();
		rate = 1.0;
		incorporateCladogenesis = cladoFlagInput.get();
		int myTotalNumberOfStates = getTotalNumberStates();

		if (incorporateCladogenesis) {
			cladoStash = cladoStashInput.get();
		}
		else { lambda = lambdaInput.get().getValues(); }
		
		// likelihood-related
		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][myTotalNumberOfStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above)
		scalingConstants = new double[tree.getNodeCount()]; // equivalent to diversitree's lq (but not in log-scale), these are used as denominators during the likelihood computation		
		Arrays.fill(scalingConstants, 1.0);		
		
		finalLk = 0.0;
		finalLogLk = 0.0;
		
		// cache-related 
		storedNodePartialScaledLksPostOde = new double[tree.getNodeCount()][myTotalNumberOfStates*2];
		storedScalingConstants = new double[tree.getNodeCount()]; 
		branchLengths = new double[tree.getNodeCount()];
		storedBranchLengths = new double[tree.getNodeCount()];
		hasDirt = Tree.IS_FILTHY;

		// thread-related
        useThreads = useThreadsInput.get() && (ProgramStatus.m_nThreads > 1);
		nrOfThreads = useThreads ? ProgramStatus.m_nThreads : 1;
		if (useThreads && maxNrOfThreadsInput.get() > 0) {
			nrOfThreads = Math.min(maxNrOfThreadsInput.get(), ProgramStatus.m_nThreads);
		}
		if (useThreads) {
		     exec = Executors.newFixedThreadPool(nrOfThreads);
		}

		// stochastic character mapping-related
		startStates = new int[tree.getNodeCount()];
		endStates = new int[tree.getNodeCount()];
		nodeConditionalScaledLks = new double[tree.getNodeCount()][myTotalNumberOfStates*2];

		sampleCharacterHistory = false;
        branchPartialLks = new ArrayList[tree.getNodeCount()];
		dt = tree.getRoot().getHeight() / numTimeSlices; // why we multiply by 50? following RevBayes code
		nodeTransitionStates = new ArrayList[tree.getNodeCount()];
		nodeTransitionTimes = new ArrayList[tree.getNodeCount()];
		nodeTimeInState = new double[tree.getNodeCount()][myTotalNumberOfStates];
		numNodeStateChanges = 0;
		// numBranchStateChanges = 0;
		averageSpeciationRates = new double[tree.getNodeCount()];
		averageExtinctionRates = new double[tree.getNodeCount()];
	}
	
	@Override
	public double calculateLogP() {
		muInput.get().getValues(mu); // every time this method is called, we need to update mu; instead of creating a new vector and assigning that to mu, we can just copy the contents of muInput into it (note that getValues is called with an argument here)
		piInput.get().getValues(pi); // same as above

		if (!incorporateCladogenesis) { lambdaInput.get().getValues(lambda); }

		// when a macroevolutionary/biogeographical parameter changes, tree is filthy, we can use threads
		// otherwise, caching allows us to only recompute part of the tree likelihood, and thread overhead not worth it
		if (hasDirt == Tree.IS_FILTHY && useThreads) {
			computeNodeLkUsingThreads();
		} else {
			computeNodeLk(tree.getRoot(), true);
		}
		logP = finalLogLk;
		return logP;
	}

    class ThreadRunner implements Runnable {

        ThreadRunner() {
        }

        @Override
		public void run() {
        	Integer nodeIdx = null;
            try {            
            	//int k = 0;
            	while (true) {
            		nodeIdx = next();
            		//k++;
            		if (nodeIdx == null) {
                        countDown.countDown();
                        //System.err.println("Done " + k);
            			return;
            		}
            		Node node = tree.getNode(nodeIdx);
            		computeNodeLk(node, false);
            		done[nodeIdx] = true;            		
            	}
            } catch (Exception e) {
                Log.err.println("Something went wrong in a calculation of " +
                		nodeIdx == null ? "something resulting in a null" : 
                			tree.getNode(nodeIdx).getID());
                e.printStackTrace();
                System.exit(1);
            }
        }

    } // CoreRunnable

    CountDownLatch countDown;
    ExecutorService exec;
    PriorityQueue<Integer> queue;
    boolean [] done;

    /**
     * comparator used by priority queue*
     */
    class TupleComparator implements Comparator<Integer> {
        @Override
		public int compare(final Integer o1, final Integer o2) {
        	Node node1 = tree.getNode(o1);
        	Node node2 = tree.getNode(o2);
        	
        	// leafs go first
        	if (node1.isLeaf() && !node2.isLeaf()) {
        		return -1;
        	}
        	if (!node1.isLeaf() && node2.isLeaf()) {
        		return 1;
        	}
        	
        	// both nodes are leafs OR both nodes are internal
            if (node1.getHeight() < node2.getHeight()) {
                return -1;
            }
            if (node1.getHeight() == node2.getHeight()) {
                return 0;
            }
            return 1;
        }
    }
    
    Integer next() {
    	synchronized(this) {
    		Integer i = queue.peek();
    		if (i == null) {
    			return null;
    		}
    		Node node = tree.getNode(i);
    		if (node.isLeaf()) {
        		return queue.poll();    			
    		}
    		if (!done[node.getChild(0).getNr()] || !done[node.getChild(1).getNr()]) {
    			return null;
    		}
    		return queue.poll();
    	}
    }

    protected void computeNodeLkUsingThreads() {
        try {
        	// set up queue
        	done = new boolean[tree.getNodeCount()];
        	queue = new PriorityQueue<>(tree.getNodeCount(), new TupleComparator());
        	for (int i = 0; i < tree.getNodeCount(); i++) {
        		queue.add(i);
        	}
        	
            countDown = new CountDownLatch(nrOfThreads);
            // kick off the threads
            for (int i = 0; i < nrOfThreads; i++) {
            	ThreadRunner coreRunnable = new ThreadRunner();
            	exec.execute(coreRunnable);
            }
            countDown.await();
        } catch (RejectedExecutionException | InterruptedException e) {
            useThreads = false;
            Log.err.println("Stop using threads: " + e.getMessage());
            computeNodeLk(tree.getRoot(), true); // in case of error, call regular non-thread computeNodeLk
        }
    }

	protected int computeNodeLk(Node node, boolean recurse) {
		int myTotalNumberOfStates = getTotalNumberStates();
		int nodeIdx = node.getNr();
		
		// cache-related stuff
        int update = (node.isDirty() | hasDirt);
        if (branchLengths[node.getNr()] != node.getLength()) {
        	update |= Tree.IS_DIRTY;
        }
        
        // may not be necessary if there is no (relaxed) clock
        branchLengths[node.getNr()] = node.getLength();

        double [] nodePartial = this.nodePartialScaledLksPostOde[nodeIdx];
        
        // leaves
        if (node.isLeaf()) {
        	if (update != Tree.IS_CLEAN) {
        		initializeLeafLks(node, nodePartial);
        	}
		}
		
		// internal nodes
		else {
			Node left = node.getChild(0);
			Node right = node.getChild(1);

			int leftIdx = left.getNr();
			int rightIdx = right.getNr();

			// recursion over here
			if (recurse) {
				// note that computeNodeLk returns update (int)
				update |= computeNodeLk(left, recurse); // if either left child dirty
				update |= computeNodeLk(right, recurse); // or right child dirty, I (parent) am dirty
			}
			
			/*
			 * if either daughter was dirty or if stochastic character mapping is active,
			 * we do stuff, otherwise skip (caching being done)
			 */
            if (update != Tree.IS_CLEAN || sampleCharacterHistory) {			
				final double[] leftLks = nodePartialScaledLksPostOde[leftIdx]; // at this point, left_lks has not yet been updated (scaled) for its parent merging
				final double[] rightLks = nodePartialScaledLksPostOde[rightIdx];

				HashMap<int[], Double> eventMap = new HashMap<int[], Double>();
				Double[] speciationRates = new Double[myTotalNumberOfStates];
				if (incorporateCladogenesis) {
					eventMap = cladoStash.getEventMap();
				}
				
				else {
					speciationRates = lambda;
				}
				
				// merge descendant lks
				for (int i = 0; i < myTotalNumberOfStates; ++i) {
					// E's
					nodePartial[i] = (leftLks[i] + rightLks[i])/2; // filling out Es using avg of children
							
					// merging with cladogenetic component
					if (incorporateCladogenesis) {
						double likeSum = 0.0;
						
						for (HashMap.Entry<int[], Double> entry: eventMap.entrySet()) {
							int[] states = entry.getKey();
							int j = states[1] - 1;
							int k = states[2] - 1;
							double speciationRate = entry.getValue();
							
							// if parent state (correcting for offset) is the state of the i-th node
							if ((states[0]-1) == i) {
								// System.out.println("Currently looking at state = " + Integer.toString(i + 1));
								// System.out.println("States are: " + Arrays.toString(states));
								// System.out.println("Lambda_ijk is: " + Double.toString(speciation_rate));
								
								// D's pre-merging steps								
								// System.out.println("Left lk array pre-scaling: " + Arrays.toString(leftLks));
								// System.out.println("Right lk array pre-scaling: " + Arrays.toString(rightLks));
								// System.out.println("Left j: " + Double.toString(leftLks[numStates + j]) + "\nRight k: " + Double.toString(rightLks[numStates + k]));
								// System.out.println("Left k: " + Double.toString(leftLks[numStates + k]) + "\nRight j: " + Double.toString(rightLks[numStates + j]));
								// System.out.println("Left scale factor: " + Double.toString(scaling_constants[leftIdx]) + "\nRight scale factor: " + Double.toString(scaling_constants[rightIdx]));
								// System.out.println("Scaled left j: " + Double.toString(leftLks[numStates + j]/scalingConstants[leftIdx]) + "\nScaled right k: " + Double.toString(right_lks[numStates + k]/scaling_constants[rightIdx]));
								// System.out.println("Scaled left k: " + Double.toString(leftLks[numStates + k]/scalingConstants[leftIdx]) + "\nScaled right j: " + Double.toString(right_lks[numStates + j]/scaling_constants[rightIdx]));
	
								/*
								 * Equation A3 (original version V0 took care of scaling here -- not anymore)
								 */
								double dnjDmk = (leftLks[myTotalNumberOfStates + j] ) * (rightLks[myTotalNumberOfStates + k]);
								double dnkDmj = (leftLks[myTotalNumberOfStates + k] ) *	(rightLks[myTotalNumberOfStates + j]);
								double dnjDmkPlusDnkDmj = dnjDmk + dnkDmj;
								likeSum += 0.5 * dnjDmkPlusDnkDmj * speciationRate;
								
								// System.out.println("Root age (height): " + Double.toString(rootAge));
								// System.out.println("Current scaled likelihood sum: " + Double.toString(likeSum));
							}
						}
						
						// finalizing merging for state Di
						nodePartial[myTotalNumberOfStates + i] = likeSum;
					}
					
					// merging w/o cladogenetic component
					else {		
						// D's pre-merging steps (scaling left and right children)
						nodePartial[myTotalNumberOfStates + i] = leftLks[myTotalNumberOfStates + i] *
								rightLks[myTotalNumberOfStates + i];
						nodePartial[myTotalNumberOfStates + i] *= speciationRates[i];
					}
				}
							
				// System.out.println("Normalized lks left: " + Arrays.toString(nodePartialScaledLksPostOde[leftIdx]));
				// System.out.println("Normalized lks right: " + Arrays.toString(nodePartialScaledLksPostOde[rightIdx]));
				// System.out.println("(Merged) Normalized lks from internal node: " + Arrays.toString(nodePartialScaledLksPostOde[nodeIdx]));
            }
		}
		
		/*
		 * not root yet, so numerical integration is carried out for all branches
		 * starting at this node up to its parent
		 */
		if (!node.isRoot()) {
            if (update != Tree.IS_CLEAN || sampleCharacterHistory) {
				// we are going from present (begin) to past (end)
				double beginAge = node.getHeight();
				double endAge = node.getParent().getHeight();

				boolean backwardTime = true;
				boolean extinctionOnly = false; // "normal" backward pass
				if (!sampleCharacterHistory) {
					numericallyIntegrateProcess(nodePartial, beginAge, endAge, backwardTime, extinctionOnly);
				}
				
				// if stochastic character mapping is active
				else {
					ArrayList<double[]> branchLks = new ArrayList<>();
					int numSteps = 0;
					double curDtStart, curDtEnd;
					
					// split the branch into small chunks to store partial lks at those chunks for forward pass later
					while (numSteps * dt + beginAge < endAge) {
						curDtStart = numSteps * dt + beginAge;
						curDtEnd = (numSteps + 1) * dt + beginAge;
						if (curDtEnd > endAge) {
							curDtEnd = endAge;
						}
						
						// "normal" backward pass, and now storing partial lks of each chunk
						numericallyIntegrateProcess(nodePartial, curDtStart, curDtEnd, backwardTime, extinctionOnly);
						double[] nodePartialCopy = new double[nodePartial.length];
						System.arraycopy(nodePartial, 0, nodePartialCopy, 0, nodePartial.length);
						branchLks.add(nodePartialCopy);
						numSteps += 1;
					}
					// TODO Unlike node likelihoods, we do not need to scale branch likelihoods since the ratios are what matter
					branchPartialLks[nodeIdx] = branchLks;
				}
				
				// scaling is done not at time of merging as original V0 version, but prior to returning when recurring
				double dScalingConstant = sum(nodePartial, myTotalNumberOfStates, nodePartial.length, -1, false); // -1 means don't ignore any item
				scalingConstants[nodeIdx] = dScalingConstant;
				for (int i = 0; i < myTotalNumberOfStates; i++) {
					nodePartial[myTotalNumberOfStates + i] /= dScalingConstant;
				}
				
            }
		}
		
		// root (no more numerical integration, children have already been joined above, now multiply by log)
		else {
			rootIdx = nodeIdx;			
			double prob = 0.0;
			
			for (int i = 0; i < myTotalNumberOfStates; ++i) {
				prob += pi[myTotalNumberOfStates + i] * nodePartial[myTotalNumberOfStates + i];
			}
			
			boolean takeLog = true;
			finalLogLk = Math.log(prob) + sum(scalingConstants, 0, scalingConstants.length, rootIdx, takeLog);
		}
		
		return update; // this is the reason why computeNodeLk isn't void() as in the original V0 version (we need update to carry out caching)
	}

	public void setIntegratorMinStep(double minStep) {
    	integratorMinStep = minStep;
	}

	public void setIntegratorTolerance(double tolerance) {
    	integratorTolerance = tolerance;
	}

	// no hidden states, so total number of states = number of obs states
	protected int getTotalNumberStates() {
		return numObsStates;
	}
	
	// used by unit tests
    public double[][] getNodePartialScaledLksPostOde() {
    	return nodePartialScaledLksPostOde;
	}
    
	protected void initializeLeafLks(Node aNode, double[] aNodePartial) {
		System.arraycopy(traitStash.getSpLks(aNode.getID()), 0, aNodePartial, 0, aNodePartial.length);

	}

	protected void numericallyIntegrateProcess(double[] likelihoods, double beginAge, double endAge, boolean backwardTime, boolean extinctionOnly) {
		if (beginAge > endAge) {
			throw new IllegalArgumentException("Improper integration. beginAge is greater than endAge");
		}
		SSEODE ode = new SSEODE(mu, q, rate, incorporateCladogenesis, backwardTime, extinctionOnly);
		solveODE(likelihoods, beginAge, endAge, ode);
	}

	protected void solveODE(double[] likelihoods, double beginAge, double endAge, SSEODE anOde) {
		FirstOrderIntegrator dp853 = new
				DormandPrince853Integrator(integratorMinStep, 100.0, integratorTolerance, integratorTolerance);

		if (incorporateCladogenesis) {
			HashMap<int[], Double> eventMap = cladoStash.getEventMap();
			anOde.setEventMap(eventMap);
			dp853.integrate(anOde, beginAge, likelihoods, endAge, likelihoods);
		}

		else {
			anOde.setSpeciationRates(lambda);
			dp853.integrate(anOde, beginAge, likelihoods, endAge, likelihoods);
		}
	}

	public int [] drawStochasticAncestralStateAtRoot() {
		/*
	     * Scaling is done in the backward pass every time we cross an internal node (because computeNodeLikelihood 
	     * scales right before returning). This means that the partial lks that we store at internal nodes for the
	     * forward pass won't suffer from underflow. No need to scale. In addition, only the ratio of state partial
	     * likelihoods matter for the sampling that is carried out during stochastic character mapping.
	     */
		// numNodeStateChanges = 0;
		// numBranchStateChanges = 0;

		// set flag for backwards pass to store likelihoods for all branch segments
        boolean sampleCharacterHistoryStashed = sampleCharacterHistory;
		sampleCharacterHistory = true;

		// backward pass
        this.calculateLogP();

		// sample root
		Node node = tree.getRoot();
		Node left = node.getChild(0);
		Node right = node.getChild(1);
		int leftIdx = left.getNr();
		int rightIdx = right.getNr();
		double[] piPrimitive = Stream.of(pi).mapToDouble(Double::doubleValue).toArray(); // just making Double[] into double[]
		int[] sampledStates = sampleAncestralState(nodePartialScaledLksPostOde[leftIdx], nodePartialScaledLksPostOde[rightIdx], piPrimitive);

		// RevBayes only tracks endStates (and call it character_histories) but I do not see why we shouldn't
		// use the paradigm in drawJointAncestralStateAtRoot
		endStates[rootIdx] = sampledStates[0]; // important one (sampled parent state, i.e., at internal node)
		startStates[leftIdx] = sampledStates[1];
		startStates[rightIdx] = sampledStates[2];

        recursivelyDrawStochasticAncestralState(left);
        recursivelyDrawStochasticAncestralState(right);

		// reset backward pass storing flag back to original
		sampleCharacterHistory = sampleCharacterHistoryStashed;
		return endStates;
	}

	private void recursivelyDrawStochasticAncestralState(Node node) {
        int nodeIdx = node.getNr();

        // Zero out E and condition D
        int startState = startStates[node.getNr()];
		
        // zero E's and condition D, then run blind backward ODE to calculate E
        initializeED(nodeConditionalScaledLks[nodeIdx], startState, false);

        /*
		 * in RB, not sure why don't store the E's from the regular backward pass... maybe
		 * this procedure is cheap enough that we can afford to do it... investigate
		 * later
		 */
        double parentAge = node.getParent().getHeight();
        numericallyIntegrateProcess(nodeConditionalScaledLks[nodeIdx], 0, parentAge, true, true);

		// set up sampling variables
		int numSteps = 0;
		double curDtStart, curDtEnd;
		int newState;
		int curState = startStates[nodeIdx]; // last sample 
		double branchLength = node.getLength();
		double[] nodeConditionalLk = nodeConditionalScaledLks[nodeIdx]; // "pi" for next step
		int numChunksInBranch = branchPartialLks[nodeIdx].size();

		// TODO check this variables are done correctly
        // set up logging/transition variables
		ArrayList<Integer> transitionStates = new ArrayList<>();
		ArrayList<Double> transitionTimes = new ArrayList<>(); // times in state before transition
		nodeTransitionStates[nodeIdx] = transitionStates;
		nodeTransitionTimes[nodeIdx] = transitionTimes;
		double[] timeInState = nodeTimeInState[nodeIdx];
		transitionStates.add(curState);
		double totalSpeciationRate = 0;
		double totalExtinctionRate = 0;
		double averageSpeciationRate, averageExtinctionRate;

		/*
		RevBayes code integrates and samples over the chunks up to the last full chunk (the leftover bit that is smaller than
		the other chunks). Instead of integrating over this small time interval (interval that is smaller than dt), they
		get the state for that chunk of time from the speciation event or tip state. They do not integrate over it either.
		We noticed that if we skip this <dt chunk of time, our posteriors get undesirable behaviors.
		Specifically, if the dt is too big, then the posteriors do not match. However, if we reduce the size of dt,
		we get matching posteriors. We do not want to have code performance dependent on this argument, thus in our
		code, we integrate over the last dt as well and handle the speciation event or tip state without sampling and
		conditioning the last chunk.
		 */
		
		// going down a branch and sampling at every chunk
		while ((numSteps+1) * dt < branchLength) {
		    // Run a small forward pass
			curDtStart = numSteps * dt;
			curDtEnd = (numSteps + 1) * dt;
			numericallyIntegrateProcess(nodeConditionalLk, curDtStart, curDtEnd, false, false); // now forward pass up to the end of this chunk

			// sample with freshly-calculated conditional and partial from backward pass
            int idxOfChunkWhosePartialToUse = numChunksInBranch - numSteps - 1; // since we are going opposite direction now (-1 is offset)
			double[] branchPartialLk = branchPartialLks[nodeIdx].get(idxOfChunkWhosePartialToUse);
            double[] branchPosterior = mergeArrays(branchPartialLk, nodeConditionalLk); // partial likelihood * "pi"
            // branchPosterior is what we will sample from (along a branch so no speciation rates involved)
            newState = sampleLksArray(branchPosterior) + 1; // +1 to correct for offset

            // log any transitions (not messing around with these variables anyway)
            if (newState != curState) {
                transitionStates.add(newState);
                double pastTransitionTimesSum = sum(transitionTimes);
				double transitionTime = curDtEnd - pastTransitionTimesSum;
                transitionTimes.add(transitionTime);
				curState = newState;
				// numBranchStateChanges++;
			}

			// Log general information (not messing around with these variables anyway)
			timeInState[curState - 1] += dt;
            totalSpeciationRate += calculateEffectiveLambda(curState);
            totalExtinctionRate += mu[curState - 1];

            // Condition D on the sampled state (but do not touch E)
			initializeED(nodeConditionalLk, curState, true);
            numSteps += 1;
		}

		// We integrate over the last chunk (then sample). This way the stochasticity of the event at the node is more pronounced
		// This bit is not in RB and seems to make a difference
		numericallyIntegrateProcess(nodeConditionalLk, numSteps * dt, branchLength, false, false);

		// we have gone down (up) a branch, and now reached a tip
		if (node.isLeaf()) {
			// last chunk will be observation
			newState = traitStash.getNodeState(node);

			// no observation, use curState for last chunk
			if (newState - 1 == -1) {
				newState = curState;
			}

			// Track transition activity
			if (newState != curState) {
				transitionStates.add(newState);
				double pastTransitionTimesSum = sum(transitionTimes);
				double transitionTime = branchLength - pastTransitionTimesSum;
				transitionTimes.add(transitionTime);
				curState = newState;
			}

			// Track branch activity
			timeInState[curState - 1] = branchLength % dt;  // the last chunk's length is the remainder
			totalSpeciationRate += calculateEffectiveLambda(curState);
			totalExtinctionRate += mu[curState - 1];

			// Finish average calculations
			averageSpeciationRate = totalSpeciationRate / numSteps;
			averageExtinctionRate = totalExtinctionRate / numSteps;
			averageSpeciationRates[nodeIdx] = averageSpeciationRate;
			averageExtinctionRates[nodeIdx] = averageExtinctionRate;

			// Set the final state
			endStates[nodeIdx] = curState;
		}
		
		// we have gone down (up) a branch, and now reached an internal node
		else {
		    // Last chunk will be sampled
			Node left = node.getChild(0);
			Node right = node.getChild(1);
			int leftIdx = left.getNr();
			int rightIdx = right.getNr();

			// sample using node likelihoods (TODO SHOULD WE BE USING THE LAST BRANCH PIECE?)
            double[] leftPartialLks = nodePartialScaledLksPostOde[leftIdx];
            double[] rightPartialLks = nodePartialScaledLksPostOde[rightIdx];
			int[] sampledStates = sampleAncestralState(leftPartialLks, rightPartialLks, nodeConditionalLk);
			newState = sampledStates[0];

			// Track transition activity
			if (curState != newState) {  //report transition
				transitionStates.add(newState);
				double pastTransitionTimesSum = sum(transitionTimes);
				double transitionTime = branchLength - pastTransitionTimesSum;
				transitionTimes.add(transitionTime);
				curState = newState;
			}

			// Track branch activity
			timeInState[curState - 1] = branchLength % dt;  // the last chunk's length
			totalSpeciationRate += calculateEffectiveLambda(curState);
			totalExtinctionRate += mu[curState - 1];

			// Finish average calculations
			averageSpeciationRate = totalSpeciationRate / numSteps;
			averageExtinctionRate = totalExtinctionRate / numSteps;
			averageSpeciationRates[nodeIdx] = averageSpeciationRate;
			averageExtinctionRates[nodeIdx] = averageExtinctionRate;

			// Set states and recurse
			endStates[nodeIdx] = curState;
			startStates[leftIdx] = sampledStates[1];
			startStates[rightIdx] = sampledStates[2];
			recursivelyDrawStochasticAncestralState(left);
			recursivelyDrawStochasticAncestralState(right);
		}
	}

	public int[] drawJointAncestralStatesAtRoot() {
	    /*
	     * Scaling is done in the backward pass every time we cross an internal node (because computeNodeLikelihood 
	     * scales right before returning). This means that the partial lks that we store at internal nodes for the
	     * forward pass won't suffer from underflow. No need to scale. In addition, only the ratio of state partial
	     * likelihoods matter for the sampling that is carried out during stochastic character mapping.
	     */
		// numNodeStateChanges = 0;
		// numBranchStateChanges = 0;

		// backward pass
    	this.calculateLogP();

		Node node = tree.getRoot();
		Node left = node.getChild(0);
		Node right = node.getChild(1);
		int leftIdx = left.getNr();
		int rightIdx = right.getNr();

		// sample root
        double[] piPrimitive = Stream.of(pi).mapToDouble(Double::doubleValue).toArray(); // just making Double[] into double[]
		int[] sampledStates = sampleAncestralState(nodePartialScaledLksPostOde[leftIdx], nodePartialScaledLksPostOde[rightIdx], piPrimitive); // sample triplet
		endStates[rootIdx] = sampledStates[0]; // important one (sampled parent state, i.e., at internal node)
		startStates[leftIdx] = sampledStates[1];
		startStates[rightIdx] = sampledStates[2];

		recursivelyDrawJointAncestralStates(left);
		recursivelyDrawJointAncestralStates(right);

		return endStates; // returning internal node sampled states (with tips too)
	}

	private void recursivelyDrawJointAncestralStates(Node node) {
    	if (node.isLeaf()) {
    	    int state = traitStash.getNodeState(node);

    	    // known tip states/observation
            if (state - 1 != -1) {
                endStates[node.getNr()] = state;
			}
            
            // should we not know the tip state (unobserved), it is -1, and so we actually sample it
            else {
				int nodeIdx = node.getNr();
				double parentAge = node.getParent().getHeight();
				double nodeAge = node.getHeight();

				int startState = startStates[node.getNr()];

				// zero E's and condition D, then run blind backward ODE to calculate E
				initializeED(nodeConditionalScaledLks[nodeIdx], startState, false);
				
				/*
				 * in RB, not sure why don't store the E's from the regular backward pass... maybe
				 * this procedure is cheap enough that we can afford to do it... investigate
				 * later
				 */
				numericallyIntegrateProcess(nodeConditionalScaledLks[nodeIdx], 0, parentAge, true, true);

				// now run forward to the tip
				numericallyIntegrateProcess(nodeConditionalScaledLks[nodeIdx], nodeAge, parentAge, false, false);

				// sample from the likelihoods for tip
				state = sampleLksArray(nodeConditionalScaledLks[nodeIdx]) + 1;
				endStates[nodeIdx] = state;
			}
		}
    	
    	// if internal node
    	else {
			int nodeIdx = node.getNr();
			Node left = node.getChild(0);
			Node right = node.getChild(1);
			int leftIdx = left.getNr();
			int rightIdx = right.getNr();
			double parentAge = node.getParent().getHeight();
			double nodeAge = node.getHeight();

			// zero Es and condition D, then run blind backward ODE to calculate E
			int startState = startStates[node.getNr()];
			initializeED(nodeConditionalScaledLks[nodeIdx], startState, false);
			
			/*
			 * in RB, not sure why don't store the E's from the regular backward pass... maybe
			 * this procedure is cheap enough that we can afford to do it... investigate
			 * later
			 */
			numericallyIntegrateProcess(nodeConditionalScaledLks[nodeIdx], 0, parentAge, true, true);

			// now run forward to the more recent internal node descending directly from current node
			numericallyIntegrateProcess(nodeConditionalScaledLks[nodeIdx], nodeAge, parentAge, false, false);

			/*
			 *  sample with freshly-calculated conditional likelihood and partial likelihoods from backwards pass;
			 *  We sample proportionally to a speciation rate * partial lks of left * partial lks of right * pi (nodeConditionalScaledLks)
			 *  This proportionality is dealt with by sampleAncestralState and differently for BiSSE/ClaSSE
			 */
			int[] sampledStates = sampleAncestralState(nodePartialScaledLksPostOde[leftIdx], nodePartialScaledLksPostOde[rightIdx], nodeConditionalScaledLks[nodeIdx]);
			endStates[nodeIdx] = sampledStates[0]; // important one (sampled parent state, i.e., at internal node)
			startStates[leftIdx] = sampledStates[1];
			startStates[rightIdx] = sampledStates[2];

			// Maybe useful at some point...
			//			if (startState != endStates[nodeIdx]) {
			//			    numNodeStateChanges++;
			//			}
			
			recursivelyDrawJointAncestralStates(left);
			recursivelyDrawJointAncestralStates(right);
		}
	}

	// helper
	public static double sum(double[] arr, int fromIdx, int toIdx, int idxToIgnore, boolean takeLog) {
		double result = 0.0;
		for (int i = fromIdx; i < toIdx; ++i) {
			if (i != idxToIgnore) {
				if (takeLog) {
					result += Math.log(arr[i]);
				}
				else { result += arr[i]; }
			}
		}
		// System.out.println("Lk array being summed to get normalizing constant: " + Arrays.toString(arr));
		// System.out.println("Sum is: " + Double.toString(result));
		return result;
	}
	
	public static double prod(double[] arr, int fromIdx, int toIdx, int idxToIgnore) {
		double result = 1.0;
		for (int i = fromIdx; i < toIdx; ++i) {
			if (i != idxToIgnore) {
				result *= arr[i];
			}
		}
		return result;
	}

	public static double sum(ArrayList<Double> arrayList) {
    	double ret = 0;
    	for (Double o: arrayList) {
    		ret += o;
		}
		return ret;
	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}
	
    @Override
	protected boolean requiresRecalculation() {
//    	if (true) {
//			hasDirt = Tree.IS_FILTHY;
//			return true;
//    	}
        hasDirt = Tree.IS_CLEAN;
		if ((irmInput.get().isDirtyCalculation()) ||
			(lambdaInput.get() != null && lambdaInput.get().somethingIsDirty()) || 
			(muInput.get().somethingIsDirty()) || 
			(piInput.get().somethingIsDirty()) ||
			(cladoStash != null && cladoStash.isDirtyCalculation())) {
			hasDirt = Tree.IS_FILTHY;
			return true;
		}
        return treeInput.get().somethingIsDirty();
        // return true;
	}
    
    @Override
    public void store() {
    	System.arraycopy(branchLengths, 0, storedBranchLengths, 0, branchLengths.length);
    	// for (int i = 0; i < nodePartialScaledLksPreOde.length; i++) {
    	//  	System.arraycopy(nodePartialScaledLksPreOde[i], 0, storedNodePartialScaledLksPreOde[i], 0, nodePartialScaledLksPreOde[i].length);    	     		
      	// }
    	
    	// before operator changes something, we store partial likelihoods of
    	// E's and D's
    	for (int i = 0; i < nodePartialScaledLksPostOde.length; i++) {
    	 	System.arraycopy(nodePartialScaledLksPostOde[i], 0, storedNodePartialScaledLksPostOde[i], 0, nodePartialScaledLksPostOde[i].length);    	     		
    	}
    	// and also store scaling constants, so that if operator changes yield worse posterior, we can go back to pre-operation
    	System.arraycopy(scalingConstants, 0, storedScalingConstants, 0, scalingConstants.length);

    	super.store();
    }
    
    @Override
    public void restore() {
    	double [] tmp = storedBranchLengths;
    	storedBranchLengths = branchLengths;
    	branchLengths = tmp;
    	
    	tmp = storedScalingConstants;
    	storedScalingConstants = scalingConstants;
    	scalingConstants = tmp;
    	
    	// double [][] tmp2 = storedNodePartialScaledLksPreOde;
    	// storedNodePartialScaledLksPreOde = nodePartialScaledLksPreOde;
    	// nodePartialScaledLksPreOde = tmp2;
    	
    	double [][] tmp2 = storedNodePartialScaledLksPostOde;
    	storedNodePartialScaledLksPostOde = nodePartialScaledLksPostOde;
    	nodePartialScaledLksPostOde = tmp2;
    	
    	super.restore();
    }

    /*
    * Used by all ASR methods for sampling states!
    * drawJointAncestralStatesAtRoot()
    * recursivelyDrawJointAncestralStates()
    * drawStochasticAncestralStatesAtRoot()
    * recursivelyDrawStochasticAncestralState()
    
    * This sampling will return a triplet in which each element is a number in {1, 2, ..., numStates}
    * since our states are indexed by 1
    */
    private int[] sampleAncestralState(double[] leftLks, double[] rightLks, double[] D) {
        // pick cladogenetic or anagenetic events
		int myTotalNumberOfStates = getTotalNumberStates();
		HashMap<int[], Double> eventMap = new HashMap<int[], Double>();
		Double[] speciationRates = new Double[myTotalNumberOfStates];
		if (incorporateCladogenesis) {
			eventMap = cladoStash.getEventMap();
		}
		else {
			speciationRates = lambda;
		}

		/*
		 * Set up the event probabilities with event <-> probability mapping
		 * where probability mapping is = speciation rate * partial lks left * partial lks right * pi (nodeConditionalScaledLks)
		 */
		HashMap<int[], Double> eventProb = new HashMap<>();
		double totalProb = 0;
        
		if (incorporateCladogenesis) {
			for (HashMap.Entry<int[], Double> entry: eventMap.entrySet()) {
				int[] states = entry.getKey();
				
				// correcting for offset (-1)
				int i = states[0] - 1; // parent state
				int j = states[1] - 1;
				int k = states[2] - 1;
				double speciationRate = entry.getValue(); // speciation rate associated to this particular kind of event
				
				// partial lks left * partial lks right * "pi"
				double lks = leftLks[myTotalNumberOfStates + j] * rightLks[myTotalNumberOfStates + k] * D[myTotalNumberOfStates + i];
				
				// now * speciation rate
				double prob = lks * speciationRate;

				// TODO: Do I need to handle left[k] and right[j]

				eventProb.put(states, prob);
				totalProb += prob; // used in sampling
			}
		}
        
		// just anagenetic
        else {
            for (int i = 0; i < myTotalNumberOfStates; i++) {
            	// partial lks left * partial lks right * "pi"
            	double lks = leftLks[myTotalNumberOfStates + i] * rightLks[myTotalNumberOfStates + i] * D[myTotalNumberOfStates + i];
            	
            	// now * speciation rate
            	double prob = 2 * lks * speciationRates[i]; // * 2 because there are we can flip left<->right
            	// (doesn't make a different to double the prob here, because all events get doubled)
            	int[] states = new int[]{i + 1, i + 1, i + 1}; // parent and children all have the same state (+1 for offset)s

				eventProb.put(states, prob);
				totalProb += prob; // used in sampling
			}
		} // finished building the event <-> probability map; we now use it to sample below

		// Sample from the events
		int[] triplet = new int[] { -1, -1, -1 };
        
		// RB: possibly dealing with underflow since they don't seem to scale
		if (totalProb <= 1e-8) {
        	// If we reach this point of the code, it is possible we need to scale our erronous code
            System.out.println("SAMPLING UNIFORMLY SINCE PROBABILITIES ARE SO SMALL EVEN THO IT MAY NOT BE RIGHT");
			int numEvents = eventProb.size();
			double randNum = Math.random(); // [0,1]
			
			// iterate over possible events
			for (HashMap.Entry<int[], Double> entry: eventProb.entrySet()) {
				triplet = entry.getKey();
				randNum -= 1.0 / numEvents;
				
				if (randNum < 0) {
					break;
				}
			}
		}
        
		// hopefully always enter here
        else {
			double randNum = Math.random() * totalProb; // Math.random() E [0,1]
			
			for (HashMap.Entry<int[], Double> entry: eventProb.entrySet()) {
				triplet = entry.getKey();
				double prob = entry.getValue();
				randNum -= prob;
				
				if (randNum < 0) {
					break;
				}
			}
		}
		
		return triplet;
	}

	/*
	Requires a likelihood array of format [E D]
    NOTE: THIS SAMPLING IS INDEX BY 0 SINCE ITS SAMPLING AN ARRAY INDEX
	IF SAMPLING STATE, PROBABLY NEED TO ADD 1 TO INDEX BY 1
    Used for branch sampling
	 */
	private int sampleLksArray(double[] lks) {
    	double totalProb = 0;
		int myTotalNumberOfStates = getTotalNumberStates();
    	for (int i = 0; i < myTotalNumberOfStates; i++) {
			totalProb += lks[myTotalNumberOfStates + i];
		}

        int ret = 0;
		if (totalProb <= 1e-8) {
			// If we reach this point of the code, it is possible we need to scale our erronous code
			System.out.println("SAMPLING UNIFORMLY SINCE PROBABILITIES ARE SO SMALL EVEN THO IT MAY NOT BE RIGHT");
			double randNum = Math.random();
			for (int i = 0; i < myTotalNumberOfStates; i++) {
				ret = i;
				randNum -= 1.0 / myTotalNumberOfStates;
				if (randNum < 0) {
					break;
				}
			}
		} else {
			double randNum = Math.random() * totalProb;
			for (int i = 0; i < myTotalNumberOfStates; i++) {
				ret = i;
				randNum -= lks[myTotalNumberOfStates + i];
				if (randNum < 0) {
					break;
				}
			}
		}
    	return ret;
	}

	/**
	 * a mapping from the internal node number (STARTING WITH INDEX 0) to the node's ID/name/label
	 * Useful for matching internal nodes with different programs who share node ID/name/label
	 * idxLabelMapper[i] is the label of node i (internal nodes indexed from 0)
	 */
	public String[] getNodeIndexNameMapper() {
    	String[] indexNameMap = new String[tree.getInternalNodeCount()];
    	Node root = tree.getRoot();
    	indexNameMap[rootIdx - tree.getLeafNodeCount()] = root.getID();

    	recursvelyGetNodeIndexNameMapper(root.getChild(0), indexNameMap);
		recursvelyGetNodeIndexNameMapper(root.getChild(1), indexNameMap);
    	return indexNameMap;
	}

	private void recursvelyGetNodeIndexNameMapper(Node node, String[] indexNameMap) {
        if (node.isLeaf()) {
        	return;
		}
		indexNameMap[node.getNr() - tree.getLeafNodeCount()] = node.getID();

		recursvelyGetNodeIndexNameMapper(node.getChild(0), indexNameMap);
		recursvelyGetNodeIndexNameMapper(node.getChild(1), indexNameMap);
	}

	/*
    Initializes extinction and likelihoods for a node (extinction likelihoods to 0, survival likelihoods to 0 except the actual state to 1)
    Useful to condition the node after sampling
	If ignoreE is true, then we do not 0 out the extinction likelihoods
	 */
	private double[] initializeED(double[] lks, int state, boolean ignoreE) {
		int myTotalNumberOfStates = getTotalNumberStates();
		
		for (int i = 0; i < myTotalNumberOfStates; i++) {
			if (!ignoreE)  {
				lks[i] = 0.0; // everyone is 0 but the sampled D
			}
			
			if (i + 1 == state) {
				lks[myTotalNumberOfStates + i] = 1.0; // sampled D
			}
			
			else {
				lks[myTotalNumberOfStates + i] = 0.0; // D's that were not sampled
			}
		}
		
		return lks;
	}

	public void setSampleCharacterHistory(boolean sch) {
    	sampleCharacterHistory = sch;
	}

	private double[] mergeArrays(double[] arr1, double[] arr2) {
	    // element-wise product
        double[] ret = new double[arr1.length];
    	for (int i = 0; i < arr1.length; i++) {
    		ret[i] = arr2[i] * arr1[i];
		}
		return ret;
	}

	/*
	 * Used in validation and unit tests
	 * Method for calling ASR many (numTrials) times
	 */
	public int[][] sampleStatesForTree(int numTrials, boolean joint) {
		int numNodes = tree.getNodeCount();
		int[][] samples = new int[numTrials][numNodes];

		// aample tree numTrial times
		for (int i = 0; i < numTrials; i++) {
			int[] drawnAncestralEnd;
		
			if (joint) {
				drawnAncestralEnd = drawJointAncestralStatesAtRoot();
			}
			
			else {
				drawnAncestralEnd = drawStochasticAncestralStateAtRoot();
			}
			
			System.arraycopy(drawnAncestralEnd, 0, samples[i], 0, numNodes);
		}

		return samples;
	}

	/*
	 * Used in BiSSE validation and unit tests
	 * Method for summarizing an array of samples (numTrials samples)
	 * Each sample has one sample state per node
	 * This function computes the frequency at which state 1 was sampled for all nodes across numTrials samples
	 */
	public double[] summarizeBiSSE(int[][] samples, int numTrials) {
		int numNodes = tree.getNodeCount();
		double[] posterior = new double[numNodes];
		int stateOneCount;
		
		for (int nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
			stateOneCount = 0;
			
			for (int nTrial = 0; nTrial < numTrials; nTrial++) {
				if (samples[nTrial][nodeIdx] == 1) {
					stateOneCount++;
				}
			}
			
			posterior[nodeIdx] = 1.0 * stateOneCount / numTrials;
		}
		
		return posterior;
	}

	/*
	 * Used in ClaSSE validation and unit tests
	 * Method for summarizing an array of samples (numTrials samples)
	 * Each sample has one sample state per node
	 * This function computes the frequency at which all state were sampled for all nodes across numTrials samples
	 */
	public double[][] summarizeCLaSSE(int[][] samples, int numTrials) {
		int myTotalNumberOfStates = getTotalNumberStates();
		int numNodes = tree.getNodeCount();
		double[][] posterior = new double[numNodes][myTotalNumberOfStates];
		
		for (int nodeIdx = 0; nodeIdx < numNodes; nodeIdx++) {
			for (int nTrial = 0; nTrial < numTrials; nTrial++) {
				int state = samples[nTrial][nodeIdx];
				posterior[nodeIdx][state - 1] ++;
			}
		}
		
		for (int i = 0; i < numNodes; i++) {
			for (int j = 0; j < myTotalNumberOfStates; j++) {
				posterior[i][j] *= 1.0 / numTrials;
			}
		}
		
		return posterior;
	}

	/*
	 * Used in BiSSE validation and unit tests
	 * Run the sampling many times (either drawJoint or drawStoc)
	 * Returns posterior probability of state 1 (the frequency state 1 was sampled in numTrials samples)
	 * of the sampling for all tips and internal nodes
	 * Important: tips and internal nodes!
	 */
	public double[] sampleAndSummarizeBiSSE(int numTrials, boolean joint) {
	    int[][] samples = sampleStatesForTree(numTrials, joint);
	    double[] posterior = summarizeBiSSE(samples, numTrials);

		if (joint) {
			System.out.println("Joint: Posterior probability of state 0: " + Arrays.toString(posterior));
		}
		
		else {
			System.out.println("Stoc: Posterior probability of state 0: " + Arrays.toString(posterior));
		}

		return posterior;
	}

	/*
	 * Used in ClaSSE validation and unit tests
	 * Run the sampling many times (either drawJoint or drawStoc)
	 * Returns posterior probability of each state (the frequency each state was sampled in numTrials samples)
	 * Important: tips and internal nodes!
	 */
	public double[][] sampleAndSummarizeCLaSSE(int numTrials, boolean joint) {
		int[][] samples = sampleStatesForTree(numTrials, joint);
		double[][] posterior = summarizeCLaSSE(samples, numTrials);

		return posterior;
	}

	public void setNumTimeSlices(int numSlices) {
        numTimeSlices = numSlices;
		dt = tree.getRoot().getHeight() / numTimeSlices * 50.0;
	}

	public double calculateEffectiveLambda(int state) {
		HashMap<int[], Double> eventMap;
		if (incorporateCladogenesis) {
			eventMap = cladoStash.getEventMap();
		} else {
		    return lambda[state - 1];
		}

		double lambda = 0;
		for (HashMap.Entry<int[], Double> entry: eventMap.entrySet()) {
			int[] states = entry.getKey();
			int i = states[0] - 1;
			double speciationRate = entry.getValue();
			if (i == state - 1) {
				lambda += speciationRate;
			}
		}
		return lambda;
	}
}
