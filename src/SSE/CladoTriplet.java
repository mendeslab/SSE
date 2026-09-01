package SSE;

import beast.base.core.BEASTObject;
import beast.base.core.Input;

public class CladoTriplet extends BEASTObject {
		
	final public Input<Integer> parentStateInput = new Input<>("parentState", "Character state of parent.");
	final public Input<Integer> leftChildInput = new Input<>("leftChildState", "Character state of left child.");
	final public Input<Integer> rightChildInput = new Input<>("rightChildState", "Character state of right child.");	
	final public Input<speciationType> speciationTypeInput = new Input<>("speciationType", "Speciation type (sympatric, subsympatric, vicariance).", speciationType.SYMPATRY, speciationType.values());
	
	public enum speciationType { SYMPATRY, SUBSYMPATRY, VICARIANCE, JUMPDISPERSAL };
	private int[] cladogeneticEvent;
	private speciationType speciationEvent;
	
	@Override
	public void initAndValidate() {
		cladogeneticEvent = new int[3];
		cladogeneticEvent[0] = parentStateInput.get();
		cladogeneticEvent[1] = leftChildInput.get();
		cladogeneticEvent[2] = rightChildInput.get();
		speciationEvent = speciationTypeInput.get();
	}
	
	public int[] getCladogeneticEvent() {
		return cladogeneticEvent;
	}
	
	public speciationType getSpeciationEvent() {
		return speciationEvent;
	}
}
