package SSE;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import com.google.common.base.Joiner;
import beast.base.inference.CalculationNode;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;

/*
 * Note: getCell is where dirty check is done
 */
public class HiddenInstantaneousRateMatrix_original extends CalculationNode {

	final public Input<Integer> NstatesInput = new Input<>("numberOfStates", "How many states or geographical ranges can affect speciation and extinction.");
	final public Input<Integer> NHiddenStatesInput = new Input<>("numberOfHiddenStates", "How many hidden states or geographical ranges can affect speciation and extinction.");
	final public Input<RealParameter> FlatQmatrixInput = new Input<>("flatQMatrix", "Array (matrix whose rows were pasted) containing the instantaneous transition rate between character states.");
	final public Input<Boolean> DisallowDoubleTransitionsInput = new Input<>("disallowDoubleTransitions", "Whether or not to set double transition parameters to zero.", Validate.REQUIRED);
	final public Input<Integer> SymmetrifyAcrossDiagonalInput = new Input<>("symmetrifyAcrossDiagonal", "Whether or not to set transition from observed to hidden the same as transitions from hidden to observed.", Validate.REQUIRED);
	final public Input<HiddenObservedStateMapper> HiddenObservedStateMapperInput = new Input<>("hiddenObsStateMapper", "Maps hidden states onto observed states and vice-versa.");
	
	private int numberOfStates; // number of observed states (never changes)
	private int numberOfHiddenStates;
	private int totalNumberOfStates;
	private boolean disallowDoubleTransitions = true;
	private int symmetrifyAcrossDiagonalStateIdx;
	private HiddenObservedStateMapper hiddenObsStateMapper;
	private Map<Integer, List<Integer>> doubleTransitionIRMCellsMap;
	private boolean ignoreDiagonal = true; // right now, always true (we never query Qij for i=j in SSEODE)
	private boolean irmDirty = true;
	private Double[] matrixContent; // these transition rates will be used to populate 'q' 2D-array
	private Double[][] q;
	private Map<Integer, int[]> realParameterToQCell;
	
	@Override
	public void initAndValidate() {
		numberOfStates = NstatesInput.get();
		numberOfHiddenStates = NHiddenStatesInput.get(); // when rjMCMC implemented, this can change at different steps
		totalNumberOfStates = numberOfHiddenStates + numberOfStates; // this will also change as a result
				
		if (numberOfHiddenStates > 0 && HiddenObservedStateMapperInput.get() == null) {
			throw new IllegalArgumentException("Number of hidden states > 0, but no mapping between observed and hidden states was found.");
		}
		
		disallowDoubleTransitions = DisallowDoubleTransitionsInput.get();
		if (numberOfHiddenStates > 0 && disallowDoubleTransitions) {
			hiddenObsStateMapper = HiddenObservedStateMapperInput.get();
		}

		symmetrifyAcrossDiagonalStateIdx = SymmetrifyAcrossDiagonalInput.get(); // if -1, no symmetrifying done
		System.out.println("My symmetrifyAcrossDiagonalStateIdx is " + symmetrifyAcrossDiagonalStateIdx);
		realParameterToQCell = new HashMap<>();
		
		Double[] someMatrixContent = FlatQmatrixInput.get().getValues();
      
		populateIRM(ignoreDiagonal, disallowDoubleTransitions, symmetrifyAcrossDiagonalStateIdx, numberOfStates, numberOfHiddenStates, someMatrixContent);
	}
	
	public void populateIRM(boolean ignoreDiagonal, boolean disallowDoubleTransitions, int symmetrifyAcrossDiagonal, int nObsStates, int nHiddenStates, Double[] someMatrixContent) {
		matrixContent = someMatrixContent; 
		numberOfHiddenStates = nHiddenStates;
		totalNumberOfStates = nObsStates + nHiddenStates; // this will also change as a result
		q = new Double[totalNumberOfStates][totalNumberOfStates]; // and so we need to vary the size of q	
		int numberOfInputElements = FlatQmatrixInput.get().getDimension();
		int diagEntry = 0;
		
		// Making sure the right number of transition parameters was provided
		if (ignoreDiagonal) {
			if (!disallowDoubleTransitions) {
				String errorMsg = "Tried to fill transition matrix Q (ignoring diagonal elements, but not double transitions), but was given ";
				if (((totalNumberOfStates) * (totalNumberOfStates) - (totalNumberOfStates)) > numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too few values.");
				} else if (((totalNumberOfStates) * (totalNumberOfStates) - (totalNumberOfStates)) < numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too many values.");
				}
			}
			
			else if ((numberOfHiddenStates > 0) && (numberOfHiddenStates == numberOfStates) && disallowDoubleTransitions){
				String errorMsg = "Tried to fill transition matrix Q with Nobs=Nhidden (ignoring diagonal elements and also double transitions), but was given ";
				if (((totalNumberOfStates) * (totalNumberOfStates) - (totalNumberOfStates) - (2 * (numberOfStates * numberOfHiddenStates - numberOfHiddenStates))) > numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too few values."); // subtract diagonal and then double transitions (twice b/c it's top-right and bottom-left)
				} else if (((totalNumberOfStates) * (totalNumberOfStates) - (totalNumberOfStates) - (2 * (numberOfStates * numberOfHiddenStates - numberOfHiddenStates))) < numberOfInputElements) {
					throw new RuntimeException(errorMsg + "too many values.");
				}
			}
		}
				
		// grabbing cells with double transitions (populating doubleTransitionIRMCellsMap), which will be ignored when populating IRM and remain 0's
		if (numberOfHiddenStates > 0 && disallowDoubleTransitions) {
			doubleTransitionIRMCellsMap = new HashMap<>();
			findDoubleTransitionIRMCells(numberOfHiddenStates, totalNumberOfStates);
			// printDoubleTransitionIRMCellsMap(); // for debugging
		}
		
		int realParameterIdx = 0;
		List<Integer> hiddenIdxList = new ArrayList<Integer>(); // will store list of hidden idx associated to observed idx (and be rewritten in a loop)	
		for (int i=0; i<totalNumberOfStates; ++i) {
			for (int j=0; j<totalNumberOfStates; ++j) {
				if (ignoreDiagonal && j == diagEntry) { 
					setCell(i, j, 0.0);
					continue;
				} // ignore diagonal element (make it zero)
				
				if (numberOfHiddenStates > 0 && disallowDoubleTransitions) {
					hiddenIdxList = doubleTransitionIRMCellsMap.get(i);
					
					if (hiddenIdxList != null) {
						if (hiddenIdxList.contains(j)) {
							setCell(i, j, 0.0);
							continue;
						}
					}
				} // ignore double transitions (make it zero)
				
				setCell(i, j, someMatrixContent[realParameterIdx]); // after setting zeros (above), we finally set the matrix cells
				
				/* 
				 * recording which cells of the Q matrix we put some non-zero value in (mapping cell to index in RealParameter list
				 * so we can use it later to resize matrixContent in rjMCMC)
				 */
				if (realParameterToQCell.get(realParameterIdx) == null) {
					int[] nonZeroIJ = {i, j};
					realParameterToQCell.put(realParameterIdx, nonZeroIJ);
					// System.out.println("HIRM: recording RealParameter transition index " + realParameterIdx + " and its cell " + Arrays.toString(nonZeroIJ));
				}
				
				realParameterIdx += 1;
			}
				 
			diagEntry += 1; // index of the diagonal element we should ignore
		}
		
		// make q_hidden2obs = q_obs2hidden for all states
		if (symmetrifyAcrossDiagonal != -1) {
			symmetrifyAcrossDiagonal(symmetrifyAcrossDiagonal);
		}
		
		irmDirty = false; // after re-population of IRM, things are clean
	}
	
	/* Later I should make sure this is made just once for each k when rjMCMC is implemented, or even when any transition rate is operated on */
	public void findDoubleTransitionIRMCells(int numberOfHiddenStates, int totalNumberOfStates) {
		
		for (int i=0; i<totalNumberOfStates; ++i) {
			// i-th row

			// recording double transition cells from top-right of IRM
			if (i < numberOfStates) {
				for (int j=numberOfStates; j<totalNumberOfStates; ++j) {
					// j-th col

					int h = hiddenObsStateMapper.getHiddenFromObs(i); // hidden state idx
					
					// this is checking the diagonal of the hidden states sub-matrix (off-diagonal are 0!)
					// if double transition; or if Nhidden < Nobs, and this obs never transitions into a hidden state (h == -1)
					if ((j - numberOfStates) != h || h == -1) {
					
						// j is hidden state idx
						if (doubleTransitionIRMCellsMap.get(i) == null) {
							List<Integer> v = new ArrayList<>();
							doubleTransitionIRMCellsMap.put(i, v); // recording
						}
						doubleTransitionIRMCellsMap.get(i).add(j); // recording
					}
					
				}
			}
			
			// recording double transition cells from bottom-left of IRM
			if (i >= numberOfStates) {
				int h = i-numberOfStates; // hidden state idx
				Collection<Integer> obsIdxColl = hiddenObsStateMapper.getObsCollFromHidden(h); // collection of observed state idxs
				
				// left cols (j) of the bottom half
				for (int j=0; j < numberOfStates; ++j) {
					/* the observed states that never transition to/from hidden states are taken care of automatically below,
					 * as their indexes will not be obsIdxColl
					 */
					if (!obsIdxColl.contains(j)) {		
						if (doubleTransitionIRMCellsMap.get(i) == null) {
							List<Integer> v = new ArrayList<>();
							doubleTransitionIRMCellsMap.put(i, v); // recording
						}
						
						doubleTransitionIRMCellsMap.get(i).add(j); // recording
					}

					// System.out.println("Bottom-left i=" + i + " j=" +j);
				}
			}
		}
					
	}
	
	/* necessary for the intermediate step in rjMCMC */
	public void symmetrifyAcrossDiagonal(int obsStateToSymmIdx) {

		System.out.println("Symmetrifying observed/hidden state " + obsStateToSymmIdx);
		
		for (int col = numberOfStates; col < totalNumberOfStates; ++col) {
			if (q[obsStateToSymmIdx][col] > 0.0) {
				setCell(col, obsStateToSymmIdx, q[obsStateToSymmIdx][col]);
			}
		}
	}
	
	public void resetQ(int aNumberOfStates, int aNumberOfHiddenStatesInMask) {
		totalNumberOfStates = aNumberOfStates + aNumberOfHiddenStatesInMask;
		q = new Double[totalNumberOfStates][totalNumberOfStates];
	}
	
	public void setCell(int from, int to, double prob) {
        this.q[from][to] = prob;
	}

	/* 
	 * called by masquerade ball during rjMCMC
	 */
	public void setHiddenObsStateAssignment(int[] aHiddenToObsStateAssignment) {
		hiddenObsStateMapper.setHiddenStateStrings(aHiddenToObsStateAssignment);
	}
	
	// getters
	public int getNumStates() {
		return numberOfStates;
	}
	
	public int getNumHiddenStates() {
		return numberOfHiddenStates;
	}
	
	public int getTotalNumStates () {
		return totalNumberOfStates;
	}
	
	public double getCell(int from, int to, double rate) {
		synchronized (this) {
			if (irmDirty) {
				Double[] someMatrixContent = FlatQmatrixInput.get().getValues();
				populateIRM(ignoreDiagonal, disallowDoubleTransitions, symmetrifyAcrossDiagonalStateIdx, numberOfStates, numberOfHiddenStates, someMatrixContent); // only re-populate IRM if some transition rate was operated on
			}			
		}
		return q[from][to] * rate;
	}

	public Double[] getMatrixContent() {
		return matrixContent;
	}
	
	public Map<Integer, int[]> getRealParameterToQCellMap() {
		return realParameterToQCell;
	}
	
	/*
	 * only used by masquerade ball (not by ODE),
	 * so does not check for dirtyness nor populateIRM
	 * (getCell is the method that does that because that's what the ODE uses)
	 */
	public Double[][] getQ() {
		return q;
	}
	
	// helper
	public void printMatrix() {
		for (int i = 0; i < totalNumberOfStates; ++i) {
			for (int j = 0; j < totalNumberOfStates; ++j) {
				// System.out.print(Double.toString(Q.getMatrixValue(i, j)) + " ");
				System.out.print(Double.toString(q[i][j]) + " ");
			}
			System.out.println();
		}
	}
	
	public void printDoubleTransitionIRMCellsMap() {
		System.out.println("IRM cells that contain double transitions and that should get 0s.");
		Iterator<Entry<Integer, List<Integer>>> iter = doubleTransitionIRMCellsMap.entrySet().iterator();
		while (iter.hasNext()) {
			Entry<Integer, List<Integer>> entry = iter.next();
			String k = Integer.toString(entry.getKey());
			String v = Joiner.on(",").join(entry.getValue());
			System.out.println("i=" + k + " j(s)=" + v);
		}
	}
	
	public void printHiddenObservedStateMapper() {
		hiddenObsStateMapper.printHidden2ObsMap();
		hiddenObsStateMapper.printObs2HiddenMap();
	}
	
	protected boolean requiresRecalculation() {
		irmDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		Double[] someMatrixContent = FlatQmatrixInput.get().getValues();
		populateIRM(ignoreDiagonal, disallowDoubleTransitions, symmetrifyAcrossDiagonalStateIdx, numberOfStates, numberOfHiddenStates, someMatrixContent);
		super.restore();
	}
}