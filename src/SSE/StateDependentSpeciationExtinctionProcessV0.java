package SSE;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Random;

import org.apache.commons.math3.ode.FirstOrderIntegrator;
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator;

import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.*;

@Description("Cladogenetic State change Speciation and Extinction (ClaSSE) model")
@Citation(value="Goldberg EE, Igic B (2012) Tempo and mode in plant breeding system evolution. Evolution 16(12):3701-3709",
year=2012, firstAuthorSurname="Goldberg", DOI="10.1111/j.1558-5646.2012.01730.x")
public class StateDependentSpeciationExtinctionProcessV0 extends Distribution {

	final public Input<Tree> treeInput = new Input<>("tree", "Tree object containing tree.", Validate.REQUIRED);
	final public Input<TraitStash> traitStashInput = new Input<>("traitStash", "TraitStash object containing the observed character state for each species.", Validate.REQUIRED);
	final public Input<InstantaneousRateMatrix> irmInput = new Input<>("instantaneousRateMatrix", "InstantaneousRateMatrix object containing anagenenetic rates.", Validate.REQUIRED);
	final public Input<CladogeneticSpeciationRateStash> cladoStashInput = new Input<>("cladogeneticStash", "CladogeneticSpeciationRateStash object that generates event map.");
	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each state (if cladogenetic events are not considered).", Validate.XOR, cladoStashInput);
	final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each state.", Validate.REQUIRED);
	final public Input<RealParameter> piInput = new Input<>("pi", "Equilibrium frequencies at root.", Validate.REQUIRED);
	final public Input<Boolean> cladoFlagInput = new Input<>("incorporateCladogenesis", "Whether or not to incorporate cladogenetic events.", Validate.REQUIRED);
	
	// input
	private Tree tree;
	private TraitStash traitStash;
	private InstantaneousRateMatrix q;
	private CladogeneticSpeciationRateStash cladoStash;
	private Double[] lambda;
	private Double[] mu;
	private Double[] pi; // root eq freqs
	private int numStates;
	private double rate;
	private double dt; // time slice size (ctor populates)
	private double rootAge;
	private int numTimeSlices;
	private boolean incorporateCladogenesis;
	
	// members used for lk computation
	private double[][] nodePartialScaledLksPreOde;
	private double[][] nodePartialScaledLksPostOde;
	private double[] scalingConstants;
	double finalLogLk;
	double finalLk;
	int rootIdx;
	
	@Override
	public void initAndValidate() {		
		super.initAndValidate();
		tree = treeInput.get();
		traitStash = traitStashInput.get();
		q = irmInput.get();
		mu = muInput.get().getValues();
		pi = piInput.get().getValues();
		numStates = q.getNumObsStates();
		rate = 1.0;
		incorporateCladogenesis = cladoFlagInput.get();
		
		if (incorporateCladogenesis) {
			cladoStash = cladoStashInput.get();
		}
		else { lambda = lambdaInput.get().getValues(); }
		
		// ode-related
		numTimeSlices = 1;
		rootAge = tree.getRoot().getHeight();
		dt = rootAge / ((double) (numTimeSlices * 50));
		
		// likelihood-related
		nodePartialScaledLksPreOde = new double[tree.getNodeCount()][numStates*2]; // tips have initialization lks, internal nodes (and root) just after merge
		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
		scalingConstants = new double[tree.getNodeCount()]; // equivalent to diversitree's lq (but not in log-scale), these are used as denominators during the likelihood computation
		Arrays.fill(scalingConstants, 1.0);		
		finalLk = 0.0;
		finalLogLk = 0.0;
	}
	
//	public StateDependentSpeciationExtinctionProcess(TreeParser tree, double[] lambda, double[] mu, double[] pi, int numStates,
//			TraitStash traitStash, CladogeneticSpeciationRateStash clado_stash, InstantaneousRateMatrix q, double rate,
//			boolean incorporateCladogenesis) {
//		this.lambda = lambda;
//		this.mu = mu;
//		this.pi = pi;
//		this.numStates = numStates;
//		this.traitStash = traitStash;
//		this.cladoStash = clado_stash;
//		Q = q;
//		this.rate = rate;
//		this.incorporateCladogenesis = incorporateCladogenesis;
//		
//		numTimeSlices = 1;
//		double rootAge = tree.getRoot().getHeight();
//		dt = rootAge / ((double) (numTimeSlices * 50));
//		// System.out.println("Root age (height): " + Double.toString(rootAge));
//		// System.out.println("dt: " + Double.toString(dt));
//		
//		// initializing members for lk computation
//		nodePartialScaledLksPreOde = new double[tree.getNodeCount()][numStates*2]; // tips have initialization lks, internal nodes (and root) just after merge
//		nodePartialScaledLksPostOde = new double[tree.getNodeCount()][numStates*2]; // tips and internal nodes have lks after the ODE went down their ancestral branches (root is special case, where it's just after merge, so the same as above) 
//		
//		scalingConstants = new double[tree.getNodeCount()];
//		Arrays.fill(scalingConstants, 1.0);
//		
//		finalLk = 0.0;
//		finalLogLk = 0.0;
//	}
	
	@Override
	public double calculateLogP() {
		muInput.get().getValues(mu); // every time this method is called, we need to update mu; instead of creating a new vector and assigning that to mu, we can just copy the contents of muInput into it (note that getValues is called with an argument here)
		piInput.get().getValues(pi); // same as above

		if (!incorporateCladogenesis) { lambdaInput.get().getValues(lambda); }

		computeNodeLk(tree.getRoot(), tree.getRoot().getNr());
		logP = finalLogLk;
		return logP;
	}

	// The smallest time slice that we will attempt to numerically integrate.
	// if the lattice of dt's falls closer than this to a node, that sliver of time is ignored in the integration.
	private static final double VERY_SMALL_TIME_SLIVER = 1e-15;

	private void computeNodeLk(Node node, int nodeIdx) {
		if (node.isLeaf()) {
			nodePartialScaledLksPreOde[nodeIdx] = traitStash.getSpLks(node.getID());
			nodePartialScaledLksPostOde[nodeIdx] = traitStash.getSpLks(node.getID()).clone();
			// System.out.println("Leaf " + node.getID() + " has node idx: " + node.getNr());
			// System.out.println("Leaf " + node.getID() + " scaled pre-ODE lks are: " + Arrays.toString(node_partial_normalized_lks_pre_ode[nodeIdx]));
			// System.out.println("Leaf " + node.getID() + " scaled post-ODE lks are: " + Arrays.toString(node_partial_normalized_lks_post_ode[nodeIdx]));
		}
		
		// this block initializes node_partial_lks for internal nodes
		else {
			Node left = node.getChild(0);
			Node right = node.getChild(1);

			int leftIdx = left.getNr();
			int rightIdx = right.getNr();

			// recursion over here
			computeNodeLk(left, leftIdx);
			computeNodeLk(right, rightIdx);
			
			// System.out.println("Recurring back to internal node " + Integer.toString(node.getNr()));
			
			double[] leftLks = nodePartialScaledLksPostOde[leftIdx].clone(); // at this point, left_lks has not yet been updated (scaled) for its parent merging
			double[] rightLks = nodePartialScaledLksPostOde[rightIdx].clone();
			
			double dLeftScalingConstant = sum(leftLks, numStates, leftLks.length, -1, false); // -1 means don't ignore any item
			scalingConstants[leftIdx] = dLeftScalingConstant;
			
			double dRightScalingConstant = sum(rightLks, numStates, rightLks.length, -1, false);
			scalingConstants[rightIdx] = dRightScalingConstant;
			
			HashMap<int[], Double> eventMap = new HashMap<int[], Double>();
			Double[] speciationRates = new Double[numStates];
			if (incorporateCladogenesis) {
				eventMap = cladoStash.getEventMap();
//				System.out.println("Event map inside computeNodeLk");
//				System.out.println(new PrettyPrintHashMap<int[], Double>(eventMap));
			}
			
			else {
				speciationRates = lambda;
			}
			
			// merge descendant lks
			for (int i = 0; i < numStates; ++i) {
				// E's
				// node_partial_unscaled_lks[nodeIdx][i] = (leftLks[i] + rightLks[i])/2; // filling out Es using avg of children
				nodePartialScaledLksPostOde[nodeIdx][i] = (leftLks[i] + rightLks[i])/2; // filling out Es using avg of children
				nodePartialScaledLksPreOde[nodeIdx][i] = nodePartialScaledLksPostOde[nodeIdx][i]; // same for pre-ODE
						
				// merging with cladogenetic component
				if (incorporateCladogenesis) {
					double likeSum = 0.0;
					
					for (HashMap.Entry<int[], Double> entry: eventMap.entrySet()) {
						int[] states = entry.getKey();
						int j = states[1] - 1;
						int k = states[2] - 1;
						double speciationRate = entry.getValue();
						
						if ((states[0]-1) == i) {
							// System.out.println("Currently looking at state = " + Integer.toString(i + 1));
							// System.out.println("States are: " + Arrays.toString(states));
							// System.out.println("Lambda_ijk is: " + Double.toString(speciation_rate));
							
							// D's pre-merging steps
							
							/*
							 * RevBayes version of Equation A3 (with scaling)
							 */
							// double likelihoods = leftLks[numStates + j] / scalingConstants[leftIdx] * 
							// 		 rightLks[numStates + k] / scalingConstants[rightIdx];
							// likeSum += likelihoods * speciationRate;
							
							// System.out.println("Left lk array pre-scaling: " + Arrays.toString(leftLks));
							// System.out.println("Right lk array pre-scaling: " + Arrays.toString(rightLks));
							// System.out.println("Left j: " + Double.toString(leftLks[numStates + j]) + "\nRight k: " + Double.toString(rightLks[numStates + k]));
							// System.out.println("Left k: " + Double.toString(leftLks[numStates + k]) + "\nRight j: " + Double.toString(rightLks[numStates + j]));
							// System.out.println("Left scale factor: " + Double.toString(scaling_constants[leftIdx]) + "\nRight scale factor: " + Double.toString(scaling_constants[rightIdx]));
							// System.out.println("Scaled left j: " + Double.toString(leftLks[numStates + j]/scalingConstants[leftIdx]) + "\nScaled right k: " + Double.toString(right_lks[numStates + k]/scaling_constants[rightIdx]));
							// System.out.println("Scaled left k: " + Double.toString(leftLks[numStates + k]/scalingConstants[leftIdx]) + "\nScaled right j: " + Double.toString(right_lks[numStates + j]/scaling_constants[rightIdx]));

							/*
							 * Equation A3
							 */
							double dnjDmk = (leftLks[numStates + j] / scalingConstants[leftIdx]) *
									(rightLks[numStates + k] / scalingConstants[rightIdx]);
							double dnkDmj = (leftLks[numStates + k] / scalingConstants[leftIdx]) *
									(rightLks[numStates + j] / scalingConstants[rightIdx]);
							double dnjDmkPlusDnkDmj = dnjDmk + dnkDmj;
							likeSum += 0.5 * dnjDmkPlusDnkDmj * speciationRate;
							
							// System.out.println("Root age (height): " + Double.toString(rootAge));
							// System.out.println("Current scaled likelihood sum: " + Double.toString(likeSum));
						}
					}

					nodePartialScaledLksPostOde[leftIdx][numStates + i] = leftLks[numStates + i] / scalingConstants[leftIdx];
					nodePartialScaledLksPostOde[rightIdx][numStates + i] = rightLks[numStates + i] / scalingConstants[rightIdx];
					
					// finalizing merging for state Di
					nodePartialScaledLksPostOde[nodeIdx][numStates + i] = likeSum;
					nodePartialScaledLksPreOde[nodeIdx][numStates + i] = nodePartialScaledLksPostOde[nodeIdx][numStates + i]; // this is a double, so no deep copy necessary
					// so pre_ode for internal nodes contains post-merging, but before ODE of 
				}
				
				// merging w/o cladogenetic component
				else {		
					// D's pre-merging steps (scaling left and right children)
					
					// scaling
					nodePartialScaledLksPostOde[leftIdx][numStates + i] = leftLks[numStates + i] /
							dLeftScalingConstant; // now we update left child lk with scaled value
					nodePartialScaledLksPostOde[rightIdx][numStates + i] = rightLks[numStates + i] / 
							dRightScalingConstant; // now we update right child lk with scaled value
					nodePartialScaledLksPostOde[nodeIdx][numStates + i] = nodePartialScaledLksPostOde[leftIdx][numStates + i] *
							nodePartialScaledLksPostOde[rightIdx][numStates + i];
					nodePartialScaledLksPostOde[nodeIdx][numStates + i] *= speciationRates[i];
					
					// keeping track of likelihoods right before ODE, at all nodes (so at internal nodes, it's post scaling and merging)
					nodePartialScaledLksPreOde[nodeIdx][numStates + i] = nodePartialScaledLksPostOde[nodeIdx][numStates + i]; // this is a double, so no deep copy necessary
				}
			}
						
			// System.out.println("Normalized lks left: " + Arrays.toString(nodePartialScaledLksPostOde[leftIdx]));
			// System.out.println("Normalized lks right: " + Arrays.toString(nodePartialScaledLksPostOde[rightIdx]));
			// System.out.println("(Merged) Normalized lks from internal node: " + Arrays.toString(nodePartialScaledLksPostOde[nodeIdx]));
		}
		
		// numerical integration is carried out for all branches starting at this node, up to its parent
		// but if root, then no more numerical integration
		if (!node.isRoot()) {
			// we are going from present (begin) to past (end)
			double beginAge = node.getHeight();
			double endAge = node.getParent().getHeight();
			
			// System.out.println("Initial conditions: " + Arrays.toString(node_partial_normalized_lks_post_ode[node_idx]));
			int currentDt = 0; // counter used to multiply dt
			while ((currentDt * dt) + beginAge < endAge) {
				double currentDtStart = (currentDt * dt) + beginAge;
				double currentDtEnd = ((currentDt + 1) * dt) + beginAge;
				
				if (currentDtEnd > endAge) {
					currentDtEnd = endAge;
				}

				double timeslice = currentDtEnd - currentDtStart;
				if (timeslice >= VERY_SMALL_TIME_SLIVER) {
					numericallyIntegrateProcess(nodePartialScaledLksPostOde[nodeIdx], currentDtStart, currentDtEnd, false, false);
				} else {
					// DO NOTHING BECAUSE TOO LITTLE TIME HAS PAST AND nodePartialScaledLksPostOde[nodeIdx] will be unaffected
				}
				
	            currentDt++;
			}
		}
		
		// if we reach root, no more numerical integration, children have already been joined above,
		// now multiply by prior, populate final_prob
		else {
			int rootIdx = nodeIdx;
			
			// for (int i = 0; i < nodePartialScaledLksPreOde.length; ++i) {
			//     System.out.println("Pre-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(nodePartialScaledLksPreOde[i]));
			// }
			// for (int i = 0; i < nodePartialScaledLksPostOde.length; ++i) {
			//     System.out.println("Post-ODE lks for node = " + Integer.toString(i) + ": " + Arrays.toString(nodePartialScaledLksPostOde[i]));
			// }
			
			double prob = 0.0;
			for (int i = 0; i < numStates; ++i) {
				prob += pi[numStates + i] * nodePartialScaledLksPostOde[nodeIdx][numStates + i];
			}
			
			boolean takeLog = true;
			finalLogLk = Math.log(prob) + sum(scalingConstants, 0, scalingConstants.length, rootIdx, takeLog);
			
			// System.out.println("Root node state = " + Double.toString(nodeIdx));
			// System.out.println("(Sum over states) Pi * lk of state = " + Double.toString(prob));
			// System.out.println("Normalizing constants = " + Arrays.toString(scalingConstants));
			// System.out.println("Lk: " + Double.toString(finalLk));
			// System.out.println("LnLk: " + Double.toString(finalLogLk));
		}
	}

	private void numericallyIntegrateProcess(double[] likelihoods, double beginAge, double endAge, boolean backwardTime, boolean extinctionOnly) {
		if (beginAge > endAge) {
			throw new IllegalArgumentException("Improper integration. beginAge is greater than endAge");
		}
		FirstOrderIntegrator dp853 = new DormandPrince853Integrator(1.0e-8, 100.0, 1.0e-10, 1.0e-10);
		SSE.SSEODE ode = new SSE.SSEODE(mu, q, rate, incorporateCladogenesis, backwardTime, extinctionOnly);

		if (incorporateCladogenesis) {
			HashMap<int[], Double> eventMap = cladoStash.getEventMap();
//			System.out.println("Event map inside ODE");
//			System.out.println(new PrettyPrintHashMap<int[], Double>(eventMap));
			
			ode.setEventMap(eventMap);
			dp853.integrate(ode, beginAge, likelihoods, endAge, likelihoods);
			// System.out.println("Conditions at time " + end_age + ": " + Arrays.toString(likelihoods));
		}
		
		else {
			ode.setSpeciationRates(lambda);
			dp853.integrate(ode, beginAge, likelihoods, endAge, likelihoods);
			// System.out.println("Conditions at time " + end_age + ": " + Arrays.toString(likelihoods));
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
}
