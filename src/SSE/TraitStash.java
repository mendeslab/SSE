package SSE;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map.Entry;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.*;

public class TraitStash extends TraitSet {

	final public Input<Integer> nStatesInput = new Input<>("numberOfStates", "How many states or geographical ranges can affect speciation and extinction.");
	
	// if Human=1,Chimp=1,Gorilla=1
	// inheriting taxonValues, map, and values variables
	// protected String[] taxonValues; // state values as str, e.g., ["1", "1", "1"]
	// Map<String, Integer> map; // "spname":sp index, e.g. {"Human":0, "Chimp":1, "Gorilla":2}
	// double[] values; // state values as doubles, e.g., [1.0, 1.0, 1.0]
	protected HashMap<String, double[]> spNameLksMap = new HashMap<String, double[]>(); // initialized by ctor
	protected int numberOfObsStates;
	
	public TraitStash() {
		traitNameInput.setRule(Input.Validate.FORBIDDEN);
		dateTimeFormatInput.setRule(Input.Validate.FORBIDDEN);
	}
	
	public void initAndValidate() {
		map = new HashMap<>();
		numberOfObsStates = nStatesInput.get();
        
		// checking
        check();
        
        // populating spNameLksMap
        for (Entry<String, Integer> entry : map.entrySet()) {
        	double[] lks = new double[numberOfObsStates*2];
			String sp_name = entry.getKey();
			int sp_idx = entry.getValue();
			spNameLksMap.put(sp_name, lks);
			// System.out.println(taxonValues[sp_idx]);
			spNameLksMap.get(sp_name)[numberOfObsStates - 1 + Integer.parseInt(taxonValues[sp_idx])] = 1.0;
        }
	}
	
	protected void check() {
		List<String> labels = taxaInput.get().asStringList();
        String[] traits = traitsInput.get().split(","); // ["Human=1", "Chimp=1", "Gorilla=1"]
        taxonValues = new String[labels.size()];
        values = new double[labels.size()];
        
        for (String trait : traits) {
            trait = trait.replaceAll("\\s+", " ");
            String[] strs = trait.split("=");
            
            if (strs.length != 2) {
                throw new IllegalArgumentException("Could not parse trait: " + trait);
            }
            
            String taxonID = normalize(strs[0]);
            int taxonNr = labels.indexOf(taxonID);
            
            if (taxonNr < 0) {
                throw new IllegalArgumentException("Trait (" + taxonID + ") is not a known taxon. Spelling error perhaps?");
            }
            
            taxonValues[taxonNr] = normalize(strs[1]);
            
            // converting state to double
            try {
            	values[taxonNr] = convertValueToDouble(taxonValues[taxonNr]);
            } catch (IllegalArgumentException ex) {
                Log.err.println("Failed to parse trait '" + taxonValues[taxonNr] + "'.");
                System.exit(1);
            }
            
            map.put(taxonID, taxonNr);
        }
        
        // num_states = numOfValues(taxonValues);
        
        // find extremes
        minValue = values[0];
        maxValue = values[0];
        for (double value : values) {
            minValue = Math.min(minValue, value);
            maxValue = Math.max(maxValue, value);
        }
        
        // minValue should be 1
        if (minValue == 0) {
        	Log.warning.println("WARNING: Trait state values should start at 1, not 0. Exiting...");
        	System.exit(1);
        }
        
        // if all species share the same state, we assume binary state
        // developers: this is useful for testing likelihood computations
        // if (numberOfStates == 1) {
        //    Log.warning.println("WARNING: All species had same trait state. Assuming binary trait.");
        //    numberOfStates = 2;
        // }
        
        // checking all distinct states are consecutive order (1,2,4 is not allowed)
         if ((double) numberOfObsStates < maxValue) {
             Log.warning.println("WARNING: The larger state is greater than number of states (i.e., you specified 3 states, but then one species is =4). Exiting...");
             System.exit(1);
         }
        
        // checking all species were specified some trait value
        for (int i = 0; i < labels.size(); i++) {
            if (taxonValues[i] == null) {
            	Log.warning.println("WARNING: No trait specified for " + labels.get(i) + ". Exiting...");
            	System.exit(1);
            }
        }
        
	}
	
	// getters and setters
	public double[] getSpLks(String spname) {
		return spNameLksMap.get(spname);
	}
	
	// helper
	public int numOfValues(String[] arr) {
	    return (int) Arrays.stream(arr).distinct().count();
	}
	
	public void printLksMap() {
		for (HashMap.Entry<String, double[]> entry : spNameLksMap.entrySet()) {	
			String spname = entry.getKey();
			double[] lks = entry.getValue();
			// System.out.println(spname + " initial lks = " + Arrays.toString(lks));
		}
	}

	public int getNodeState(Node node) {
	   return (int) values[node.getNr()];
    }
}
