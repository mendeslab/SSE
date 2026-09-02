package drivers;

import java.util.Arrays;

import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import operators.StateMaskOperator;

public class MaskOperatorTestDriver {

	public static void main(String[] args) {
		Integer[] init = new Integer[3];
        Arrays.fill(init, 0);
        IntegerParameter mask = new IntegerParameter(init);
        
		State state = new State();
		state.initByName("stateNode", mask);
		state.initialise();
		
		StateMaskOperator operator = new StateMaskOperator();
        operator.initByName("mask", mask, "weight", 1.0);
        
        int[][] count = new int[2][3]; // two obs/hidden states, each 0.0/1.0/2.0
        
		for (int i = 0; i < 1000000; i++) {
			operator.proposal();
            Integer[] aMask = mask.getValues();
            
            for (int pos = 0; pos < aMask.length-1; pos++) {
                double j = (double) aMask[pos];
                int v = (int) j; // 0 or 1 or 2
                count[pos][v] += 1; 
            }
		}
		
		for (int pos = 0; pos < count.length; pos++) {
			for (int v=0; v < count[pos].length; v++) {
				System.out.println("pos = " + pos + " v = " + v + " count = " + count[pos][v]);
			}
		}
	}
}
