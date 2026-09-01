package operators;

import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;

public class StateMaskOperator extends Operator {
    final public Input<IntegerParameter> stateMaskInput =
            new Input<>("mask", "the array of doubles to operate a random walk on.", Validate.REQUIRED);
    
	@Override
	public void initAndValidate() {
	}

	@Override
	public double proposal() {
		final IntegerParameter mask = stateMaskInput.get();
		final int pos = Randomizer.nextInt(mask.getDimension());
		final double value = mask.getValue(pos);
		
		// if at this position is 0 or 2, we can only (and will!) go to 1
		if (value == 0 || value == 2) {
			mask.setValue(pos, 1);
			return Math.log(0.5);
		}
		
		else {
			if (Math.random() < 0.5) {
				mask.setValue(pos, 0); // half chance of going 0.0
			}
			else { mask.setValue(pos, 2); } // half chance of going 2.0
				
			return Math.log(2.0); // if value == 1.0
		}
	}
}
