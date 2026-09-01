package SSE;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import SSE.CladoTriplet.speciationType;
import beast.base.inference.CalculationNode;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

public class CladogeneticSpeciationRateStash extends CalculationNode {
	
	final public Input<List<CladoTriplet>> cladogeneticEventsInput = new Input<>("cladoTriplets", "List of CladoTriplet objects that contain states of a parent and its children.", new ArrayList<>());
	final public Input<RealParameter> sympatricRateInput = new Input<>("sympatricRate", "Speciation rate for sympatric event.");
	final public Input<RealParameter> subSympatricRateInput = new Input<>("subsympatricRate", "Speciation rate for subsympatric event.");
	final public Input<RealParameter> vicariantRateInput = new Input<>("vicariantRate", "Speciation rate for vicariant event.");
	final public Input<RealParameter> jumpRateInput = new Input<>("jumpRate", "Speciation rate for jump dispersal event.");
	
	private double sympatricRate;
	private double subSympatricRate;
	private double vicariantRate;
	private double jumpRate;
	private HashMap<int[], Double> eventMap = new HashMap<>(); // ctor populates this
	private boolean stashDirty = true;
	
	@Override
	public void initAndValidate() {
		populateStash();
	}

	public void populateStash() {
		sympatricRate = sympatricRateInput.get().getValue();
		subSympatricRate = subSympatricRateInput.get().getValue();
		vicariantRate = vicariantRateInput.get().getValue();
		
		if (jumpRateInput.get() != null) {
			jumpRate = jumpRateInput.get().getValue();
		}

		List<CladoTriplet> cladoTripletList = cladogeneticEventsInput.get();

		for (CladoTriplet cladoTriplet : cladoTripletList) {
			int[] cladogeneticEvent = cladoTriplet.getCladogeneticEvent();
			speciationType speciationEvent = cladoTriplet.getSpeciationEvent();

			switch (speciationEvent) {
				case SYMPATRY:
					eventMap.put(cladogeneticEvent, sympatricRate);
					break;
				case SUBSYMPATRY:
					eventMap.put(cladogeneticEvent, subSympatricRate);
					break;
				case VICARIANCE:
					eventMap.put(cladogeneticEvent, vicariantRate);
					break;
				case JUMPDISPERSAL:
					eventMap.put(cladogeneticEvent, jumpRate);
					break;
			}
		}
		
		stashDirty = false; // after re-population of stash, things are clean
	}

	protected boolean requiresRecalculation() {
		stashDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		populateStash();
		super.restore();
	}

	// setters and getters
	HashMap<int[], Double> getEventMap() {
		if (stashDirty) {
			populateStash(); // only re-populate stash if some speciation rate was operated on
		}
		return eventMap;
	}
	
	// helper
	public void printEventMap() {
		System.out.println(new PrettyPrintHashMap<int[], Double>(eventMap));
	}
}