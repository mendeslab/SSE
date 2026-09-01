package SSE;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.*;
import java.util.regex.Pattern;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import beast.base.core.BEASTObject;
import beast.base.core.Input;

public class HiddenObservedStateMapper extends BEASTObject {

	final public Input<String> hiddenStatesInFocalOrderInput = new Input<>("hiddenStates", "Comma-separated integer strings corresponding to hidden states, one per focal states (in the order of focal states).");

	private Map<Integer, Integer> obs2HiddenMap; // makeMaps populates this
	private Multimap<Integer, Integer> hidden2ObsMap;
	private int[] hiddenStateAssignments;
	Pattern comma;
	
	@Override
	public void initAndValidate() {
		comma = Pattern.compile(",");
		String hiddenStateString = hiddenStatesInFocalOrderInput.get();
		hiddenStateAssignments = comma.splitAsStream(hiddenStateString).mapToInt(Integer::parseInt).toArray();
		
		makeMaps();
	}

	public void makeMaps() {
		obs2HiddenMap = new HashMap<>();
		hidden2ObsMap = HashMultimap.create();
		
		for (int i=0; i<hiddenStateAssignments.length; i++) {
			obs2HiddenMap.put(i, hiddenStateAssignments[i]);
		}
		
		// printObs2HiddenMap(); // for debugging
		
		for (Entry<Integer, Integer> entry: obs2HiddenMap.entrySet()) {
			hidden2ObsMap.put(entry.getValue(), entry.getKey());
		}
		
		// printHidden2ObsMap(); // for debugging

	}
	
	// getters
	public int getHiddenFromObs(int obsIdx) {
		return obs2HiddenMap.get(obsIdx);
	}
	
	public Collection<Integer> getObsCollFromHidden(int hiddenIdx) {
		return hidden2ObsMap.get(hiddenIdx);
	}
	
	// setters
	public void setHiddenStateStrings(int[] aHiddenStateArray) {
		hiddenStateAssignments = aHiddenStateArray;
		// System.out.println("Updated hiddenStateAssignments to " + Arrays.toString(hiddenStateAssignments));
		makeMaps();
	}
	
	// helper
	public void printObs2HiddenMap() {
		for (Integer obs: obs2HiddenMap.keySet()) {
			String obsString = Integer.toString(obs);
			String hiddenString = obs2HiddenMap.get(obs).toString();
			// System.out.println("Obs=" + obsString + " hidden=" + hiddenString);
		}
	}
	
	public void printHidden2ObsMap() {
		for (Entry<Integer, Collection<Integer>> entry: hidden2ObsMap.asMap().entrySet()) {
			System.out.println("hidden=" + entry.getKey() + " obs=" + entry.getValue());
		}
	}
}
