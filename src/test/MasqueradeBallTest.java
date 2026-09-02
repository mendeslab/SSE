package test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import SSE.HiddenInstantaneousRateMatrix;
import SSE.HiddenObservedStateMapper;
import SSE.HiddenTraitStash;
import SSE.LambdaMuAssigner;
import SSE.MasqueradeBall;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;

public class MasqueradeBallTest {
	
	private MasqueradeBall maskBall;
	private TaxonSet taxonSet;
	
	String hiddenStatesString;
	String flatQMatrixString;
	String lambdasToStatesString;
	String musToStatesString;
	RealParameter lambda;
	RealParameter mu;
	RealParameter pi;
	boolean disallowDoubleTransitions;
	int symmetrifyAcrossDiagonal;
	int numberOfStates;
	int numberOfHiddenStates;
	
	IntegerParameter stateMask1;
	IntegerParameter cidMask1;
	IntegerParameter stateMask2;
	IntegerParameter cidMask2;
	IntegerParameter stateMask3;
	IntegerParameter cidMask3;
	IntegerParameter stateMask4;
	IntegerParameter cidMask4;
	IntegerParameter stateMask5;
	IntegerParameter cidMask5;
	IntegerParameter stateMask6;
	IntegerParameter cidMask6;
	IntegerParameter stateMask7;
	IntegerParameter cidMask7;
	IntegerParameter stateMask8;
	IntegerParameter cidMask8;
	IntegerParameter stateMask9;
	IntegerParameter cidMask9;
	IntegerParameter stateMask10;
	IntegerParameter cidMask10;
	IntegerParameter stateMask11;
	IntegerParameter cidMask11;
	IntegerParameter stateMask12;
	IntegerParameter cidMask12;
	ArrayList<IntegerParameter> stateMasks;
	ArrayList<IntegerParameter> cidMasks;

	RealParameter mask1;
	RealParameter mask2;
	RealParameter mask3;
	RealParameter mask4;
	RealParameter mask5;
	RealParameter mask6;
	RealParameter mask7;
	RealParameter mask8;
	RealParameter mask9;
	RealParameter mask10;
	RealParameter mask11;
	RealParameter mask12;
	ArrayList<RealParameter> masks;
	
	public interface Instance {
		Double[] getSp1Lk();
		
		Double[] getSp2Lk();
		
		Double[] getSp3Lk();
		
		Double[] getPis();
		
		Double[] getLambdas();

        Double[] getMus();

        Double[][] getQ();
    }
	
	Instance test1 = new Instance() {
		@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0 };
		}
		
		@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3 }, 
				{ 1.0, 0.0, 1.2, 1.3 }, 
				{ 2.0, 2.1, 0.0, 2.3 }, 
				{ 3.0, 3.1, 3.2, 0.0 } 
			};
		}
    };
    
    Instance test2 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
		}
		
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0,
					1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1, .2 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03, .06 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.4 }, 
				{ 1.0, 0.0, 1.2, 1.3, 0.0 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0 },
				{ 0.4, 0.0, 0.0, 0.0, 0.0 }
			};
		}
    };
    
    Instance test3 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
		}
    	
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1, .2, .3 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03, .06, .09 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 0.0, 1.5 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0, 0.0 },
				{ 0.4, 0.0, 0.0, 0.0, 0.0, 4.5 },
				{ 0.0, 1.5, 0.0, 0.0, 5.4, 0.0 }
			};
		}
    };
    
    Instance test4 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
		}
		
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					1.0/7.0, 1.0/7.0, 1.0/7.0, 1.0/7.0, 1.0/7.0, 1.0/7.0, 1.0/7.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1, .2, .3, .4 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03, .06, .09, .12 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.0, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 0.0, 1.5, 0.0 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0, 0.0, 2.6 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0, 0.0, 0.0 },
				{ 0.4, 0.0, 0.0, 0.0, 0.0, 4.5, 4.6 },
				{ 0.0, 1.5, 0.0, 0.0, 5.4, 0.0, 5.6 },
				{ 0.0, 0.0, 2.6, 0.0, 6.4, 6.5, 0.0 },
			};
		}
    };
    
    Instance test5 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
		
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1, .2, .3, .4, .2 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03, .06, .09, .12, .06 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.0, 0.0, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 0.0, 1.5, 0.0, 0.0 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0, 0.0, 2.6, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0, 0.0, 0.0, 3.7 },
				{ 0.4, 0.0, 0.0, 0.0, 0.0, 4.5, 4.6, 4.7 },
				{ 0.0, 1.5, 0.0, 0.0, 5.4, 0.0, 5.6, 5.7 },
				{ 0.0, 0.0, 2.6, 0.0, 6.4, 6.5, 0.0, 6.7 },
				{ 0.0, 0.0, 0.0, 3.7, 7.4, 7.5, 7.6, 0.0 }
			};
		}
    };
    
    Instance test6 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
    	
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .1, .1, .1, .2, .2, .2, .2 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .03, .03, .03, .06, .06, .06, .06 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.0, 0.0, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 0.0, 1.5, 0.0, 0.0 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0, 0.0, 2.6, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0, 0.0, 0.0, 3.7 },
				{ 0.4, 0.0, 0.0, 0.0, 0.0, 4.5, 4.6, 4.7 },
				{ 0.0, 1.5, 0.0, 0.0, 5.4, 0.0, 5.6, 5.7 },
				{ 0.0, 0.0, 2.6, 0.0, 6.4, 6.5, 0.0, 6.7 },
				{ 0.0, 0.0, 0.0, 3.7, 7.4, 7.5, 7.6, 0.0 }
			};
		}
    };
    
    Instance test7 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0 };
		}
    	
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0,
					1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0, 1.0/5.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1, .3 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03, .09 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 1.5 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0 },
				{ 0.0, 5.1, 0.0, 0.0, 0.0 }
			};
		}
    };
    
    Instance test8 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
		}
		
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .15, .2, .1, .2, .3 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .045, .06, .03, .06, .09 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 0.0, 1.5 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0, 0.0 },
				{ 0.4, 0.0, 0.0, 0.0, 0.0, 4.5 },
				{ 0.0, 5.1, 0.0, 0.0, 5.4, 0.0 }
			};
		}
    };
    
    Instance test9 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0 };
		}
    	
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .1, .1, .1, .3, .3 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .03, .03, .03, .09, .09 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.0, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 1.5, 0.0 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0, 3.7 },
				{ 0.0, 5.1, 0.0, 0.0, 0.0, 5.7 },
				{ 0.0, 0.0, 0.0, 3.7, 7.5, 0.0 }
			};
		}
    };
    
    Instance test10 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0 };
		}
    	
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 
					1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0, 1.0/6.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .1, .1, .1, .3, .3 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .03, .03, .03, .09, .09 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.0, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 1.5, 0.0 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0, 3.7 },
				{ 0.0, 5.1, 0.0, 0.0, 0.0, 5.7 },
				{ 0.0, 0.0, 0.0, 7.3, 7.5, 0.0 }
			};
		}
    };
    
    Instance test11 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
    	
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 1.0/4.0, 1.0/4.0, 1.0/4.0, 1.0/4.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .1, .1, .1 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .03, .03, .03 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3 }, 
				{ 1.0, 0.0, 1.2, 1.3 }, 
				{ 2.0, 2.1, 0.0, 2.3 }, 
				{ 3.0, 3.1, 3.2, 0.0 }
			};
		}
    };
    
    Instance test12 = new Instance() {
    	@Override
		public Double[] getSp1Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp2Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
		}
		
		@Override
		public Double[] getSp3Lk() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0 };
		}
    	
    	@Override
		public Double[] getPis() {
			return new Double[] { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0, 1.0/8.0 };
		}
    	
    	@Override
		public Double[] getLambdas() {
			return new Double[] { .1, .1, .1, .1, .2, .2, .2, .2 };
		}

		@Override
		public Double[] getMus() {
			return new Double[] { .03, .03, .03, .03, .06, .06, .06, .06 };
		}

		@Override
		public Double[][] getQ() {
			return new Double[][] { 
				{ 0.0, 0.1, 0.2, 0.3, 0.4, 0.0, 0.0, 0.0 }, 
				{ 1.0, 0.0, 1.2, 1.3, 0.0, 1.5, 0.0, 0.0 }, 
				{ 2.0, 2.1, 0.0, 2.3, 0.0, 0.0, 2.6, 0.0 }, 
				{ 3.0, 3.1, 3.2, 0.0, 0.0, 0.0, 0.0, 3.7 },
				{ 4.0, 0.0, 0.0, 0.0, 0.0, 4.5, 4.6, 4.7 },
				{ 0.0, 5.1, 0.0, 0.0, 5.4, 0.0, 5.6, 5.7 },
				{ 0.0, 0.0, 6.2, 0.0, 6.4, 6.5, 0.0, 6.7 },
				{ 0.0, 0.0, 0.0, 7.3, 7.4, 7.5, 7.6, 0.0 }
			};
		}
    };
    
	@Before
	public void setUp() throws Exception {		
		String[] spNames = new String[] { "sp1", "sp2", "sp3" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		taxonSet = new TaxonSet(taxaList);
		
		hiddenStatesString = "0,1,2,3";
				
		lambdasToStatesString = "0,1,2,3,4,5,6,7";
		Double lambda1 = 0.1; // 0A
		Double lambda2 = 0.15; // 1A
		Double lambda3 = 0.2; // 2A
		Double lambda4 = 0.1; // 3A
		Double lambda5 = 0.2; // 0B
		Double lambda6 = 0.3; // 1B
		Double lambda7 = 0.4; // 2B
		Double lambda8 = 0.2; // 3B
		Double[] lambdas = { lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, lambda7, lambda8 };
		lambda = new RealParameter(lambdas);
		
		musToStatesString = "0,1,2,3,4,5,6,7";
		Double mu1 = 0.03;
		Double mu2 = 0.045;
		Double mu3 = 0.06;
		Double mu4 = 0.03;
		Double mu5 = 0.06;
		Double mu6 = 0.09;
		Double mu7 = 0.12;
		Double mu8 = 0.06;
		Double[] mus = { mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8 };
		mu = new RealParameter(mus);
		
		Double[] pis = { 0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.25, 0.25 };
		pi = new RealParameter(pis);
		
		disallowDoubleTransitions = true; // not used
		symmetrifyAcrossDiagonal = -1;
		flatQMatrixString = "0.1 0.2 0.3 0.4 1.0 1.2 1.3 1.5 2.0 2.1 2.3 2.6 3.0 3.1 3.2 3.7 4.0 4.5 4.6 4.7 5.1 5.4 5.6 5.7 6.2 6.4 6.5 6.7 7.3 7.4 7.5 7.6";
		
//		Double[] mask1Array = { 0.0, 0.0, 0.0, 0.0, 0.0 }; // first four states, then CID/Not-CID	
//		mask1 = new RealParameter(mask1Array);
		Integer[] statesMask1Array = { 0, 0, 0, 0 };
		Integer[] cidMask1Array = { 0 };
		stateMask1 = new IntegerParameter(statesMask1Array);
		cidMask1 = new IntegerParameter(cidMask1Array);
		
//		Double[] mask2Array = { 1.0, 0.0, 0.0, 0.0, 0.0 };
//		mask2 = new RealParameter(mask2Array);
		Integer[] statesMask2Array = { 1, 0, 0, 0 };
		Integer[] cidMask2Array = { 0 };
		stateMask2 = new IntegerParameter(statesMask2Array);
		cidMask2 = new IntegerParameter(cidMask2Array);
		
//		Double[] mask3Array = { 1.0, 1.0, 0.0, 0.0, 0.0 };
//		mask3 = new RealParameter(mask3Array);
		Integer[] statesMask3Array = { 1, 1, 0, 0 };
		Integer[] cidMask3Array = { 0 };
		stateMask3 = new IntegerParameter(statesMask3Array);
		cidMask3 = new IntegerParameter(cidMask3Array);
		
//		Double[] mask4Array = { 1.0, 1.0, 1.0, 0.0, 0.0 };
//		mask4 = new RealParameter(mask4Array);
		Integer[] statesMask4Array = { 1, 1, 1, 0 };
		Integer[] cidMask4Array = { 0 };
		stateMask4 = new IntegerParameter(statesMask4Array);
		cidMask4 = new IntegerParameter(cidMask4Array);
		
//		Double[] mask5Array = { 1.0, 1.0, 1.0, 1.0, 0.0 };
//		mask5 = new RealParameter(mask5Array);
		Integer[] statesMask5Array = { 1, 1, 1, 1 };
		Integer[] cidMask5Array = { 0 };
		stateMask5 = new IntegerParameter(statesMask5Array);
		cidMask5 = new IntegerParameter(cidMask5Array);
		
//		Double[] mask6Array = { 1.0, 1.0, 1.0, 1.0, 1.0 };
//		mask6 = new RealParameter(mask6Array);
		Integer[] statesMask6Array = { 1, 1, 1, 1 };
		Integer[] cidMask6Array = { 1 };
		stateMask6 = new IntegerParameter(statesMask6Array);
		cidMask6 = new IntegerParameter(cidMask6Array);
		
//		Double[] mask7Array = { 0.0, 2.0, 0.0, 0.0, 0.0 };
//		mask7 = new RealParameter(mask7Array);
		Integer[] statesMask7Array = { 0, 2, 0, 0 };
		Integer[] cidMask7Array = { 0 };
		stateMask7 = new IntegerParameter(statesMask7Array);
		cidMask7 = new IntegerParameter(cidMask7Array);
		
//		Double[] mask8Array = { 1.0, 2.0, 0.0, 0.0, 0.0 };
//		mask8 = new RealParameter(mask8Array);
		Integer[] statesMask8Array = { 1, 2, 0, 0 };
		Integer[] cidMask8Array = { 0 };
		stateMask8 = new IntegerParameter(statesMask8Array);
		cidMask8 = new IntegerParameter(cidMask8Array);
		
//		Double[] mask9Array = { 0.0, 2.0, 0.0, 1.0, 1.0 };
//		mask9 = new RealParameter(mask9Array);
		Integer[] statesMask9Array = { 0, 2, 0, 1 };
		Integer[] cidMask9Array = { 1 };
		stateMask9 = new IntegerParameter(statesMask9Array);
		cidMask9 = new IntegerParameter(cidMask9Array);
		
//		Double[] mask10Array = { 0.0, 2.0, 0.0, 2.0, 1.0 };
//		mask10 = new RealParameter(mask10Array);
		Integer[] statesMask10Array = { 0, 2, 0, 2 };
		Integer[] cidMask10Array = { 1 };
		stateMask10 = new IntegerParameter(statesMask10Array);
		cidMask10 = new IntegerParameter(cidMask10Array);
		
//		Double[] mask11Array = { 0.0, 0.0, 0.0, 0.0, 1.0 };
//		mask11 = new RealParameter(mask11Array);
		Integer[] statesMask11Array = { 0, 0, 0, 0 };
		Integer[] cidMask11Array = { 1 };
		stateMask11 = new IntegerParameter(statesMask11Array);
		cidMask11 = new IntegerParameter(cidMask11Array);
		
//		Double[] mask12Array = { 2.0, 2.0, 2.0, 2.0, 1.0 };
//		mask12 = new RealParameter(mask12Array);
		Integer[] statesMask12Array = { 2, 2, 2, 2 };
		Integer[] cidMask12Array = { 1 };
		stateMask12 = new IntegerParameter(statesMask12Array);
		cidMask12 = new IntegerParameter(cidMask12Array);
		
//		masks = new ArrayList<RealParameter>();
//		masks.add(mask1);
//		masks.add(mask2);
//		masks.add(mask3);
//		masks.add(mask4);
//		masks.add(mask5);
//		masks.add(mask6);
//		masks.add(mask7);
//		masks.add(mask8);
//		masks.add(mask9);
//		masks.add(mask10);
//		masks.add(mask11);
//		masks.add(mask12);
		
		stateMasks = new ArrayList<IntegerParameter>();
		stateMasks.add(stateMask1);
		stateMasks.add(stateMask2);
		stateMasks.add(stateMask3);
		stateMasks.add(stateMask4);
		stateMasks.add(stateMask5);
		stateMasks.add(stateMask6);
		stateMasks.add(stateMask7);
		stateMasks.add(stateMask8);
		stateMasks.add(stateMask9);
		stateMasks.add(stateMask10);
		stateMasks.add(stateMask11);
		stateMasks.add(stateMask12);
		
		cidMasks = new ArrayList<IntegerParameter>();
		cidMasks.add(cidMask1);
		cidMasks.add(cidMask2);
		cidMasks.add(cidMask3);
		cidMasks.add(cidMask4);
		cidMasks.add(cidMask5);
		cidMasks.add(cidMask6);
		cidMasks.add(cidMask7);
		cidMasks.add(cidMask8);
		cidMasks.add(cidMask9);
		cidMasks.add(cidMask10);
		cidMasks.add(cidMask11);
		cidMasks.add(cidMask12);
	}

	Instance[] all = { test1, test2, test3, test4, test5, test6, test7, test8, test9, test10, test11, test12 };
	
	@Test
	public void testMaskeradeBall() {
		
		int i = 0;
		for (Instance test : all) {
			
			System.out.println();
			System.out.println("DOING TEST " + (i+1));
			System.out.println();
			
			LambdaMuAssigner lambdaMuAssigner = new LambdaMuAssigner();
			lambdaMuAssigner.initByName("totalNumberOfStates", 8, "nDistinctLambdas", 8, "nDistinctMus", 8, "lambdasToStates", lambdasToStatesString, "lambda", lambda, "musToStates", musToStatesString, "mu", mu, "pi", pi);

			HiddenObservedStateMapper stateMapper = new HiddenObservedStateMapper();
			stateMapper.initByName("hiddenStates", hiddenStatesString);
			
			HiddenInstantaneousRateMatrix hirm = new HiddenInstantaneousRateMatrix();
			hirm.initByName("numberOfStates", 4, "numberOfHiddenStates", 4, "flatQMatrix", flatQMatrixString, "disallowDoubleTransitions", disallowDoubleTransitions, "symmetrifyAcrossDiagonal", symmetrifyAcrossDiagonal, "hiddenObsStateMapper", stateMapper);
			
			HiddenTraitStash hts = new HiddenTraitStash();
			hts.initByName("numberOfStates", 4, "numberOfHiddenStates", 4, "taxa", taxonSet, "hiddenObsStateMapper", stateMapper, "value", "sp1=2,sp2=1,sp3=2");

			maskBall = new MasqueradeBall();
			maskBall.initByName(
					"stateMask", stateMasks.get(i),
					"cidMask", cidMasks.get(i),
					"hiddenInstantaneousRateMatrix", hirm, 
					"lambdaMuAssigner", lambdaMuAssigner,
					"hiddenTraitStash", hts);
			
			hts = maskBall.getHTS();
			
			Assert.assertArrayEquals(test.getPis(), maskBall.getPis());
			Assert.assertArrayEquals(test.getLambdas(), maskBall.getLambdas());
			Assert.assertArrayEquals(test.getMus(), maskBall.getMus());
			Assert.assertArrayEquals(test.getQ(), maskBall.getQs());
			Assert.assertArrayEquals(test.getSp1Lk(), ArrayUtils.toObject(hts.getSpLks("sp1")));
			Assert.assertArrayEquals(test.getSp2Lk(), ArrayUtils.toObject(hts.getSpLks("sp2")));
			Assert.assertArrayEquals(test.getSp3Lk(), ArrayUtils.toObject(hts.getSpLks("sp3")));

			++i;
		}
	}
	
}
