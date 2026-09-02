package SSE;

import java.util.regex.Pattern;
import beast.base.inference.CalculationNode;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

public class LambdaMuAssigner extends CalculationNode {	
	
	final public Input<Integer> totalNstatesInput = new Input<>("totalNumberOfStates", "How many states (observed + hidden) or geographical ranges can affect speciation and extinction.");
	final public Input<RealParameter> piInput = new Input<>("pi", "Equilibrium frequencies of each state.", Validate.REQUIRED);
	final public Input<Integer> nDistinctMusInput = new Input<>("nDistinctMus", "How many distinct mu values.");
	final public Input<String> musToStatesAssignerStringInput = new Input<>("musToStates", "Comma-separated integer string, one mu per state.");
	final public Input<RealParameter> muInput = new Input<>("mu", "Death rates for each state.", Validate.REQUIRED);	
	final public Input<CladogeneticSpeciationRateStash> cladoStashInput = new Input<>("cladogeneticStash", "CladogeneticSpeciationRateStash object that generates event map.");
	final public Input<RealParameter> lambdaInput = new Input<>("lambda", "Speciation rates for each state (if cladogenetic events are not considered).", Validate.XOR, cladoStashInput);
	final public Input<String> lambdasToStatesAssignerStringInput = new Input<>("lambdasToStates", "Comma-separated integer string, one lambda per state.", Validate.XOR, cladoStashInput);
	final public Input<Integer> nDistinctLambdasInput = new Input<>("nDistinctLambdas", "How many distinct lambda values.", Validate.XOR, cladoStashInput);
	
	private int totalNumberOfStates;
	private int numberOfDistinctLambdas;
	private int numberOfDistinctMus;
	private String lambdaToStatesString;
	private String muToStatesString;
	private int[] lambdaAssignments;
	private int[] muAssignments;
	private Double[] pisContent;
	private Double[] lambdasContent;
	private Double[] musContent;
	private Double[] pi;
	private Double[] lambda;
	private Double[] mu;
	private boolean assignerDirty = true;
	Pattern comma;
		
	// mcmc
	private int storedTotalNumberOfStates;
	private int storedNumberOfDistinctLambdas;
	private int storedNumberOfDistinctMus;
	private Double[] storedLambdasContent;
	private Double[] storedMusContent;
	private int[] storedLambdaAssignments;
	private int[] storedMuAssignments;
	
	@Override
	public void initAndValidate() {
		comma = Pattern.compile(",");

		/* new stuff */
		if (cladoStashInput.get() == null) {
			lambdaToStatesString = lambdasToStatesAssignerStringInput.get();
			lambdaAssignments = comma.splitAsStream(lambdaToStatesString).mapToInt(Integer::parseInt).toArray(); // comma-separated string comes from xml (done only once)
		}
		
		muToStatesString = musToStatesAssignerStringInput.get();
		muAssignments = comma.splitAsStream(muToStatesString).mapToInt(Integer::parseInt).toArray();
		
		totalNumberOfStates = totalNstatesInput.get();
		
		// updating lambdas (no incorporate cladogenesis support for now)
		if (cladoStashInput.get() == null) {
			numberOfDistinctLambdas = nDistinctLambdasInput.get();
			lambdasContent = lambdaInput.get().getValues();
			storedLambdasContent = new Double[lambdasContent.length];
			storedLambdaAssignments = new int[lambdaAssignments.length];
		}
		
		numberOfDistinctMus = nDistinctMusInput.get();
		musContent = muInput.get().getValues();
		
		// mcmc
		storedMusContent = new Double[musContent.length];
		storedMuAssignments = new int[muAssignments.length];
		
		populateAssigner(totalNumberOfStates, numberOfDistinctLambdas, numberOfDistinctMus, lambdasContent, musContent, lambdaAssignments, muAssignments);
		/* end new stuff */
		
		// populateAssigner(); // original
	}
	
	/*
	 * called by initAndValidate and getters when assigner is dirty 
	 */
	public void populateAssigner() {
		if (cladoStashInput.get() == null) {
			lambdaToStatesString = lambdasToStatesAssignerStringInput.get();
			lambdaAssignments = comma.splitAsStream(lambdaToStatesString).mapToInt(Integer::parseInt).toArray(); // comma-separated string comes from xml (done only once)
		}
		
		muToStatesString = musToStatesAssignerStringInput.get();
		muAssignments = comma.splitAsStream(muToStatesString).mapToInt(Integer::parseInt).toArray();
		
		totalNumberOfStates = totalNstatesInput.get();
		
		// updating pis
		pi = new Double[2 * totalNumberOfStates];
		pisContent = piInput.get().getValues();
		updatePis(pisContent);
		
		// updating lambdas (no incorporate cladogenesis support for now)
		if (cladoStashInput.get() == null) {
			numberOfDistinctLambdas = nDistinctLambdasInput.get();
			lambdasContent = lambdaInput.get().getValues();
			lambda = new Double[totalNumberOfStates];
			updateLambdasNoClado(lambdasContent);
		}
				
		// updating mus
		mu = new Double[totalNumberOfStates];
		numberOfDistinctMus = nDistinctMusInput.get();
		musContent = muInput.get().getValues();
		updateMus(musContent);
		
		assignerDirty = false; // we got the new values, not dirty anymore
	}
	
	/*
	 * called by masquerade ball when applying mask (and updating assigner)
	 */
	public void populateAssigner(int totalNStates, int nDistinctLambdas, int nDistinctMus, Double[] aLambdaContent, Double[] aMuContent, int[] aLambdaAssignment, int[] aMuAssignment) {
		totalNumberOfStates = totalNStates;
		numberOfDistinctLambdas = nDistinctLambdas;
		numberOfDistinctMus = nDistinctMus;
		lambdaAssignments = aLambdaAssignment;
		muAssignments = aMuAssignment;
		
		/*
		 *  updating pis (for now, the eq freq distribution is always flat)
		 *  note that I'm doing this everytime I call applyMask (and in some cases, when
		 *  totalNumberOfStates isn't changed, I'm doing this without having to)
		 */
		pi = new Double[2 * totalNumberOfStates]; // E's and D's
		pisContent = new Double[2 * totalNumberOfStates]; // E's and D's
		for (int i = 0; i<pisContent.length; ++i) {
			if (i < totalNumberOfStates) { 
				pisContent[i] = 0.0; // E's
			}
			else {
				pisContent[i] = 1.0 / totalNumberOfStates;	
			}
		}
		updatePis(pisContent);
		
		// updating lambdas (no incorporate cladogenesis support for now)
		if (cladoStashInput.get() == null) {
			lambda = new Double[totalNumberOfStates];
			lambdasContent = aLambdaContent;
			updateLambdasNoClado(aLambdaContent);
		}
		
		// updating mus
		mu = new Double[totalNumberOfStates];
		musContent = aMuContent;
		updateMus(aMuContent);
		
		assignerDirty = false; // we got the new values, not dirty anymore
	}
	
	// for now, always flat eq freq distribution, not doing complex assignments
	public void updatePis(Double[] aPiContent) {
		for (int i = 0; i<pisContent.length; ++i) {
			pi[i] = aPiContent[i];
		}
	}
	
	public void updateLambdasNoClado(Double[] aLambdaContent) {	
		for (int i = 0; i<totalNumberOfStates; ++i) {
			int lambdaAssignmentIdx = lambdaAssignments[i];
			lambda[i] = aLambdaContent[lambdaAssignmentIdx];
		}
	}
	
	public void updateMus(Double[] aMuContent) {					
		for (int i = 0; i<totalNumberOfStates; ++i) {
			int muAssignmentIdx = muAssignments[i];
			mu[i] = aMuContent[muAssignmentIdx];
		}
	}
	
	// getters
	// for (H)SDSEP requiresRecalculation()
	public RealParameter getPisRealParameter() {
		return piInput.get();
	}
	
	public Double[] getPis() {
		// if pis were operated on, things are dirty, we need to update (this isn't supported yet anyway)
		if (assignerDirty) {
			populateAssigner();
		}
		
		return pi;
	}
	
	// for (H)SDSEP requiresRecalculation()
	public RealParameter getLambdasRealParameter() {
		return lambdaInput.get();
	}
	
	// returns only lambdas that matter (after applying mask, for example)
	public Double[] getLambdas() {
		// if lambdas were operated on, things are dirty, we need to update
		if (assignerDirty) {
			populateAssigner();
		}
		
		return lambda;
	}
	
	// returns what was obtained from Input (always the same number of lambdas)
	public Double[] getLambdaContent() {
//		return lambdasContent; // original, working
		return lambdaInput.get().getValues();
	}
	
	// for (H)SDSEP requiresRecalculation()
	public RealParameter getMusRealParameter() {
		return muInput.get();
	}
		
	// returns only mus that matter (after applying mask, for example)
	public Double[] getMus() {
		// if mus were operated on, things are dirty, we need to update
		if (assignerDirty) {
			populateAssigner();
		}
		
		return mu;
	}
	
	// returns what was obtained from Input (always the same number of mus)
	public Double[] getMuContent() {
		// return musContent; // original, working
		return muInput.get().getValues();
	}
	
	public CladogeneticSpeciationRateStash getCladoStash() {
		return cladoStashInput.get();
	}
	
	// getters
	public int getNDistinctLambdas() {
		return numberOfDistinctLambdas;
	}
	
	public int getNDistinctMus() {
		return numberOfDistinctMus;
	}
	
	// mcmc
	protected boolean requiresRecalculation() {
		assignerDirty = true;
		return super.requiresRecalculation();
	}

	@Override
	protected void store() {
		// Cladogenetic models keep speciation rates in their stash, so this assigner has no lambda arrays.
		if (cladoStashInput.get() == null) {
			System.arraycopy(lambdasContent, 0, storedLambdasContent, 0, lambdasContent.length);
			System.arraycopy(lambdaAssignments, 0, storedLambdaAssignments, 0, lambdaAssignments.length);
		}
		System.arraycopy(musContent, 0, storedMusContent, 0, musContent.length);
		System.arraycopy(muAssignments, 0, storedMuAssignments, 0, muAssignments.length);
		storedTotalNumberOfStates = totalNumberOfStates;
		storedNumberOfDistinctLambdas = numberOfDistinctLambdas;
		storedNumberOfDistinctMus = numberOfDistinctMus;
		super.store();
	}
	
	// original, when working (no store necessary when using this restore)
//	@Override
//	protected void restore() {
//		populateAssigner();
//		super.restore();
//	}
	
	@Override
	protected void restore() {    	
		populateAssigner(storedTotalNumberOfStates, storedNumberOfDistinctLambdas, storedNumberOfDistinctMus, storedLambdasContent, storedMusContent, storedLambdaAssignments, storedMuAssignments);
		super.restore();
	}
}
