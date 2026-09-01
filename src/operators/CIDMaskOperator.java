package operators;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;

public class CIDMaskOperator extends Operator {
    final public Input<IntegerParameter> cidMaskInput =
            new Input<>("mask", "integer (1=CID, 0=Not-CID) to operate on.", Validate.REQUIRED);
    
	@Override
	public void initAndValidate() {
	}

	@Override
	public double proposal() {
		final IntegerParameter mask = cidMaskInput.get();
		final double value = mask.getValue(0);
		
		// if at this position is 0 or 2, we can only (and will!) go to 1
		if (value == 0) {
			mask.setValue(0, 1);
		} else {
			mask.setValue(0, 0);
		}
		
		return 0.0;
	}
}
