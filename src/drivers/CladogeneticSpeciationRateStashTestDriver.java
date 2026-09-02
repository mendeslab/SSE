package drivers;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import SSE.CladoTriplet;
import SSE.CladogeneticSpeciationRateStash;
import SSE.CladoTriplet.speciationType;
import beast.base.inference.parameter.RealParameter;

public class CladogeneticSpeciationRateStashTestDriver {

	public static void main(String[] args) {
		CladoTriplet sTriplet = new CladoTriplet();
		sTriplet.initByName("parentState", 1,
				"leftChildState", 1,
				"rightChildState", 1,
				"speciationType", speciationType.SYMPATRY);
		
		CladoTriplet ssTriplet = new CladoTriplet();
		ssTriplet.initByName("parentState", 2,
				"leftChildState", 2,
				"rightChildState", 1,
				"speciationType", speciationType.SUBSYMPATRY);
		
		List<CladoTriplet> cladoTripletList = new ArrayList<CladoTriplet>();
		Collections.addAll(cladoTripletList, sTriplet, ssTriplet);
		
		Double[] sSpeciationRate = {0.1};
		Double[] ssSpeciationRate = {0.2};
		Double[] vSpeciationRate = {0.0};
		Double[] jSpeciationRate = {0.0};
		RealParameter sympatricSpeciationRate = new RealParameter(sSpeciationRate);
		RealParameter subSympatricSpeciationRate = new RealParameter(ssSpeciationRate);
		RealParameter vicariantSpeciationRate = new RealParameter(vSpeciationRate);
		RealParameter jumpSpeciationRate = new RealParameter(jSpeciationRate);

		CladogeneticSpeciationRateStash csrt = new CladogeneticSpeciationRateStash();
		csrt.initByName("cladoTriplets", cladoTripletList,
				"sympatricRate", sympatricSpeciationRate,
				"subsympatricRate", subSympatricSpeciationRate,
				"vicariantRate", vicariantSpeciationRate,
				"jumpRate", jumpSpeciationRate);
		csrt.printEventMap();	
	}
}

