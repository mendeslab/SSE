package SSE;

import beast.base.inference.CalculationNode;
import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;

public class InstantaneousRateMatrix extends CalculationNode {

	final public Input<Integer> NstatesInput = new Input<>("numberOfStates", "How many states or geographical ranges can affect speciation and extinction.");
	final public Input<RealParameter> FlatQmatrixInput = new Input<>("flatQMatrix", "Array (matrix whose rows were pasted) containing the instantaneous transition rate between character states.");
	
	protected int numberOfStates;
	protected boolean ignoreDiagonal = true; // right now, always true (we never query Qij for i=j in SSEODE)
	protected boolean irmDirty = true;
	protected Double[][] q;
	
	@Override
	public void initAndValidate() {
		numberOfStates = NstatesInput.get();
		int numberOfInputElements = FlatQmatrixInput.get().getDimension();
		q = new Double[numberOfStates][numberOfStates];
		
		// Making sure the right number of transition parameters was provided
		if (ignoreDiagonal) {
			String errorMsg = "Tried to fill transition matrix Q (ignoring diagonal elements), but was given ";
			
			if ((numberOfStates*numberOfStates - numberOfStates) > numberOfInputElements) {
				throw new RuntimeException(errorMsg + "too few values.");
			} else if ((numberOfStates*numberOfStates - numberOfStates) < numberOfInputElements) {
				throw new RuntimeException(errorMsg + "too many values.");
			}
		}

		populateIRM(ignoreDiagonal);
		// Q.initByName("minordimension", numberOfStates);
	}
		
	public void populateIRM(boolean ignoreDiagonal) {
		Double[] matrixContent = FlatQmatrixInput.get().getValues();

		int q = 0;
		int diagEntry = 0;
		for (int i=0; i<numberOfStates; ++i) {
			for (int j=0; j<numberOfStates; ++j) {
				if (ignoreDiagonal && j == diagEntry) { 
					setCell(i, j, 0.0);
					continue;
				} // ignore diagonal element (make it zero)
				
				setCell(i, j, matrixContent[q]);
				q += 1;
			}
				 
			diagEntry += 1; // index of element we should ignore
		}
		
		irmDirty = false; // after re-population of IRM, things are clean
	}
	
	public void setCell(int from, int to, double prob) {
		// Q.setMatrixValue(from, to, prob);
        this.q[from][to] = prob;
	}

	// getters
	public int getNumObsStates() {
		return numberOfStates;
	}
	
	public double getCell(int from, int to, double rate) {
		// return Q.getMatrixValue(from, to) * rate;
		if (irmDirty) {
			populateIRM(ignoreDiagonal); // only re-populate IRM if some transition rate was operated on
		}
		return q[from][to] * rate;
	}
	
	// helper
	public void printMatrix() {
		for (int i = 0; i < numberOfStates; ++i) {
			for (int j = 0; j < numberOfStates; ++j) {
				// System.out.print(Double.toString(Q.getMatrixValue(i, j)) + " ");
				System.out.print(Double.toString(q[i][j]) + " ");
			}
			System.out.println();
		}
	}
	
	protected boolean requiresRecalculation() {
		irmDirty = true;
		return super.requiresRecalculation();
	}

	protected void restore() {
		populateIRM(ignoreDiagonal);
		super.restore();
	}
}