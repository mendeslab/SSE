package SSE;

import java.util.HashMap;
import java.util.Map.Entry;

import beast.base.core.Input;

public class HiddenTraitStash extends TraitStash implements Cloneable {

	final public Input<Integer> nHiddenStatesInput = new Input<>("numberOfHiddenStates", "How many hidden states or geographical ranges can affect speciation and extinction.");
	final public Input<HiddenObservedStateMapper> HiddenObservedStateMapperInput = new Input<>("hiddenObsStateMapper", "Maps hidden states onto observed states and vice-versa.");
	
	private int numberOfHiddenStates;
	private int numberOfStates;
	private HiddenObservedStateMapper hiddenObsStateMapper;
	
	public HiddenTraitStash() {
		traitNameInput.setRule(Input.Validate.FORBIDDEN);
		dateTimeFormatInput.setRule(Input.Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
		map = new HashMap<>();
		numberOfObsStates = nStatesInput.get();
		numberOfHiddenStates = nHiddenStatesInput.get();
		numberOfStates = numberOfObsStates + numberOfHiddenStates;
		
		if (numberOfHiddenStates > 0 && HiddenObservedStateMapperInput.get() == null) {
			throw new IllegalArgumentException("Number of hidden states > 0, but no mapping between observed and hidden states was found.");
		}
		
		if (numberOfHiddenStates > 0) { hiddenObsStateMapper = HiddenObservedStateMapperInput.get(); }
		
        check();
        
        populateSpLksMap(numberOfStates, numberOfObsStates, numberOfHiddenStates);
	}
	
	public void populateSpLksMap(int aNumberOfStates, int aNumberOfObsStates, int aNumberOfHiddenStates) {
        numberOfStates = aNumberOfStates;
        numberOfObsStates = aNumberOfObsStates;
        numberOfHiddenStates = aNumberOfHiddenStates;
        
		int hiddenIdx;
        for (Entry<String, Integer> entry : map.entrySet()) {
        	double[] lks = new double[numberOfStates*2];
			String spName = entry.getKey();
			int spIdx = entry.getValue();
			spNameLksMap.put(spName, lks);
			spNameLksMap.get(spName)[numberOfStates - 1 + Integer.parseInt(taxonValues[spIdx])] = 1.0;
			
			if (numberOfHiddenStates > 0) {
				hiddenIdx = hiddenObsStateMapper.getHiddenFromObs(Integer.parseInt(taxonValues[spIdx])-1);
				
				// only if hiddenIdx is matched with some observed state
				if (hiddenIdx != -1) {
					spNameLksMap.get(spName)[numberOfStates + numberOfObsStates + hiddenIdx] = 1.0;
				}
			} // hidden states match observed states according to hiddenObsStateMapper
        } 
	}
	
	public void setHiddenStateAssignment(int[] aHiddenObsStateAssignment) {
		hiddenObsStateMapper.setHiddenStateStrings(aHiddenObsStateAssignment);
	}
}
