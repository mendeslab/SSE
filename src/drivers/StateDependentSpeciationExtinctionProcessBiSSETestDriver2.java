package drivers;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.ArrayUtils;

import SSE.InstantaneousRateMatrix;
import SSE.StateDependentSpeciationExtinctionProcess;
import SSE.TraitStash;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.TreeParser;

public class StateDependentSpeciationExtinctionProcessBiSSETestDriver2 {

	public static void main(String[] args) {
		
		// initializing parameter values
		int numberOfStates = 2; // BiSSE
		String[] spNames = new String[] { "sp4", "sp7", "sp9", "sp11", "sp14", "sp16", "sp21", "sp22", "sp25", "sp27", "sp29", "sp30", "sp33", "sp35", "sp36", "sp37", "sp38", "sp39", "sp42", "sp45", "sp46", "sp48", "sp50", "sp51", "sp52", "sp53", "sp54", "sp55", "sp57", "sp58", "sp59", "sp61", "sp64", "sp67", "sp68", "sp69", "sp70", "sp71", "sp72", "sp73", "sp74", "sp75", "sp79", "sp80", "sp81", "sp82", "sp84", "sp85", "sp86", "sp87", "sp88", "sp89", "sp90", "sp91", "sp92", "sp93", "sp94", "sp95", "sp96", "sp97", "sp98", "sp99", "sp100", "sp101", "sp102", "sp103", "sp104", "sp105", "sp106", "sp107", "sp108", "sp109", "sp110", "sp111", "sp112", "sp113", "sp114", "sp115", "sp117", "sp118", "sp119", "sp121", "sp122", "sp123", "sp124", "sp125", "sp126", "sp127", "sp128", "sp129", "sp130", "sp131", "sp132", "sp133", "sp134", "sp135", "sp136", "sp137", "sp138", "sp139", "sp140", "sp141", "sp142", "sp143", "sp144", "sp145", "sp146", "sp148", "sp149", "sp150", "sp151", "sp152", "sp153", "sp154", "sp155", "sp156", "sp158", "sp159", "sp160", "sp161", "sp162", "sp163", "sp164", "sp165", "sp166", "sp167", "sp168", "sp169", "sp170", "sp171", "sp172" };
		List<Taxon> taxaList = Taxon.createTaxonList(Arrays.asList(spNames));
		TaxonSet taxonSet = new TaxonSet(taxaList);
		TraitStash traitStash = new TraitStash();
		traitStash.initByName("numberOfStates", numberOfStates, "taxa", taxonSet, "value", "sp4=2,sp7=2,sp9=2,sp11=1,sp14=1,sp16=2,sp21=2,sp22=2,sp25=1,sp27=2,sp29=1,sp30=1,sp33=2,sp35=1,sp36=2,sp37=1,sp38=1,sp39=1,sp42=1,sp45=1,sp46=1,sp48=1,sp50=1,sp51=2,sp52=2,sp53=1,sp54=1,sp55=1,sp57=1,sp58=1,sp59=1,sp61=1,sp64=1,sp67=1,sp68=1,sp69=1,sp70=1,sp71=2,sp72=1,sp73=1,sp74=1,sp75=1,sp79=1,sp80=1,sp81=1,sp82=1,sp84=1,sp85=1,sp86=1,sp87=1,sp88=1,sp89=1,sp90=1,sp91=1,sp92=1,sp93=1,sp94=1,sp95=1,sp96=1,sp97=1,sp98=1,sp99=1,sp100=1,sp101=2,sp102=2,sp103=1,sp104=1,sp105=1,sp106=1,sp107=1,sp108=1,sp109=1,sp110=1,sp111=2,sp112=1,sp113=1,sp114=2,sp115=2,sp117=1,sp118=1,sp119=1,sp121=1,sp122=1,sp123=1,sp124=1,sp125=1,sp126=1,sp127=1,sp128=1,sp129=1,sp130=1,sp131=1,sp132=1,sp133=1,sp134=1,sp135=1,sp136=1,sp137=1,sp138=1,sp139=1,sp140=1,sp141=1,sp142=1,sp143=1,sp144=1,sp145=1,sp146=1,sp148=1,sp149=1,sp150=1,sp151=1,sp152=1,sp153=1,sp154=1,sp155=1,sp156=1,sp158=1,sp159=1,sp160=1,sp161=1,sp162=1,sp163=1,sp164=1,sp165=1,sp166=1,sp167=1,sp168=1,sp169=2,sp170=2,sp171=1,sp172=1");
		traitStash.printLksMap();
		
		Double birthRate = 0.222222222;
		Double deathRate = 0.1;
		Double[] mus = { deathRate, deathRate };
		System.out.println("Mus: " + Arrays.toString(mus));
		RealParameter mu = new RealParameter(mus);
		mu.initByName("minordimension", 1);
		
		Double[] lambdas = new Double[numberOfStates];
		Arrays.fill(lambdas, birthRate);
		RealParameter lambda = new RealParameter(lambdas);
		
		InstantaneousRateMatrix irm = new InstantaneousRateMatrix();
		irm.initByName("numberOfStates", numberOfStates, "flatQMatrix", "0.00112606341068786 0.00532672696281225");
		irm.printMatrix();
		
		Double[] piEs = new Double[numberOfStates];
		Arrays.fill(piEs, 0.0);
		Double[] piDs = new Double[numberOfStates];
		Arrays.fill(piDs, (1.0/numberOfStates));
		Double[] pis = ArrayUtils.addAll(piEs, piDs); // 0.0, 0.0, 0.5, 0.5
		System.out.println("Pi is: " + Arrays.toString(pis));
		RealParameter pi = new RealParameter(pis);
		pi.initByName("minordimension", 1);
		
		String treeStr = "(((((((sp35:11.2660568,(sp111:3.029999957,(sp130:1.678957889,sp131:1.678957889)nd162:1.351042068)nd96:8.236056845)nd76:3.287545126,(sp106:3.125173623,sp107:3.125173623)nd139:11.4284283)nd48:3.031435786,((sp25:14.41837763,(sp158:4.247559594,(((sp88:3.628751564,sp89:3.628751564)nd153:0.07168634326,(sp142:1.194129785,(sp153:0.7749388406,sp154:0.7749388406)nd171:0.4191909443)nd154:2.506308122)nd148:0.2297187887,(sp95:3.567812585,(sp125:1.910737645,sp126:1.910737645)nd156:1.657074939)nd149:0.3623441109)nd146:0.3174028984)nd78:10.17081803)nd67:0.7786887545,((sp54:7.595075704,sp117:7.595075704)nd84:5.375713622,((sp108:3.081812333,sp109:3.081812333)nd94:8.454780771,(sp167:0.1415819879,sp168:0.1415819879)nd95:11.39501112)nd85:1.434196222)nd68:2.226277053)nd49:2.387971333)nd39:3.084754819,sp14:20.66979253)nd10:19.8075379,((((sp80:4.618930136,(sp161:0.4888022647,sp162:0.4888022647)nd141:4.130127871)nd74:9.94429016,((sp143:1.165243043,sp144:1.165243043)nd155:2.437754153,sp90:3.602997196)nd86:10.9602231)nd60:1.092121885,(sp39:10.95759413,sp55:10.95759413)nd64:4.69774805)nd36:7.547880034,sp11:23.20322221)nd11:17.27410822)nd5:4.588299098,(((((sp57:7.575160858,sp58:7.575160858)nd40:15.66986694,((((sp71:6.197423101,(sp104:3.201930505,sp105:3.201930505)nd134:2.995492596)nd119:1.822520215,(sp69:6.260205479,sp70:6.260205479)nd120:1.759737837)nd111:1.212104964,sp48:9.23204828)nd91:9.125027173,((sp169:0.1265415494,sp170:0.1265415494)nd46:17.72345537,(sp46:13.63429359,(((sp163:0.2461234975,sp164:0.2461234975)nd169:1.367739936,sp132:1.613863434)nd89:10.91075561,((sp114:2.791598369,sp115:2.791598369)nd129:4.415921015,sp148:7.207519384)nd90:5.317099657)nd83:1.109674549)nd47:4.215703327)nd43:0.5070785359)nd30:4.88795235)nd16:7.186427673,(sp7:28.11747419,(((sp50:10.24167628,(((((sp123:2.175626885,sp124:2.175626885)nd144:2.423144159,sp81:4.598771045)nd138:0.722572006,sp75:5.321343051)nd137:0.7050515173,sp72:6.026394568)nd118:2.418937933,sp51:8.445332501)nd103:1.796343776)nd99:0.380084675,sp42:10.62176095)nd88:1.999802207,sp33:12.62156316)nd22:15.49591103)nd17:2.313981284)nd14:5.556044373,(((((((sp96:3.503382882,sp121:3.503382882)nd147:0.5948884911,sp84:4.098271373)nd132:3.319621945,(sp99:3.447749971,(sp118:2.673762934,sp119:2.673762934)nd158:0.7739870375)nd127:3.970143347)nd107:2.160504188,(sp59:7.298751579,(sp91:3.592804082,sp92:3.592804082)nd128:3.705947497)nd108:2.279645927)nd100:1.427454773,sp38:11.00585228)nd87:7.174946028,((sp101:3.32129392,sp102:3.32129392)nd54:12.88926118,(((((sp127:1.906218245,(sp140:1.213169745,sp141:1.213169745)nd168:0.6930485006)nd152:1.819465491,sp87:3.725683736)nd121:4.173248486,((sp103:3.304189812,(sp112:2.979259865,(sp149:0.8678359409,sp150:0.8678359409)nd163:2.111423924)nd160:0.3249299472)nd159:0.03886028826,sp100:3.343050101)nd122:4.555882121)nd114:0.6116589098,((sp128:1.772027808,sp129:1.772027808)nd165:0.7852728188,(sp171:0.07870700064,sp172:0.07870700064)nd166:2.478593626)nd140:5.953290506)nd65:6.697195178,sp27:15.20778631)nd55:1.002768789)nd45:1.970243209)nd23:9.118785461,(((sp52:7.789274111,(sp151:0.8039157969,sp152:0.8039157969)nd123:6.985358314)nd71:8.768872893,(((((sp133:1.541132206,sp134:1.541132206)nd161:1.535137313,sp110:3.076269519)nd124:4.675831111,sp53:7.752100629)nd105:1.885857559,sp79:9.637958188)nd62:5.979717457,(sp29:12.9456946,sp30:12.9456946)nd63:2.671981042)nd53:0.9404713593)nd50:0.1381960076,((sp93:3.590827513,sp94:3.590827513)nd72:11.03058536,((((sp135:1.358573873,sp136:1.358573873)nd116:7.101003744,(sp67:6.302563771,sp68:6.302563771)nd117:2.157013846)nd112:0.4225232452,sp61:8.882100862)nd109:5.389559377,(sp36:11.07871682,sp37:11.07871682)nd80:3.192943422)nd73:0.3497526307)nd69:2.074930142)nd38:10.60324076)nd15:8.68791608)nd7:5.108194193,((sp16:18.54666543,(sp21:15.72262716,sp22:15.72262716)nd41:2.824038273)nd18:11.78954007,(sp9:21.76971128,(sp73:5.752079865,sp74:5.752079865)nd81:16.01763142)nd19:8.566494225)nd8:10.75948853)nd6:3.969935486)nd3:0.7802161345,((((sp145:0.9925019153,sp146:0.9925019153)nd92:10.95453752,(sp45:10.14505727,(sp64:6.624713652,((sp155:0.7112266349,sp156:0.7112266349)nd135:5.335832559,(((sp86:3.760078256,((sp159:0.5275647987,sp160:0.5275647987)nd170:0.7069288561,sp139:1.234493655)nd151:2.525584601)nd150:0.1619695445,sp85:3.9220478)nd143:0.6610233846,sp82:4.583071185)nd136:1.463988009)nd133:0.5776544577)nd104:3.520343618)nd93:1.801982161)nd27:11.8679967,(((sp137:1.307786203,sp138:1.307786203)nd167:0.9893437962,sp122:2.297129999)nd164:0.622433374,sp113:2.919563373)nd131:20.89547276)nd12:15.08079005,(sp4:26.14399682,((sp165:0.1764771701,sp166:0.1764771701)nd101:15.6684002,(sp97:3.478361503,sp98:3.478361503)nd58:12.36651586)nd25:10.29911945)nd13:12.75182936)nd9:6.950019481)nd2;"
;
        TreeParser myTree = new TreeParser(treeStr, false, false, true, 0); // true b/c species are labelled, offset=0
        
        boolean incorporateCladogenesis = false;
                
        StateDependentSpeciationExtinctionProcess sdsep = new StateDependentSpeciationExtinctionProcess();
        sdsep.initByName(
        		"tree", myTree,
        		"traitStash", traitStash,
        		"instantaneousRateMatrix", irm,
        		"lambda", lambda,
        		"mu", mu,
        		"pi", pi,
        		"incorporateCladogenesis", incorporateCladogenesis
        		);
    	
    	System.out.println(sdsep.calculateLogP());
	}
}
