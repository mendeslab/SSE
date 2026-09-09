package test;

import SSE.*;
import beast.base.inference.State;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.distribution.Normal;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import org.junit.Assert;
import org.junit.After;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;

import static SSE.SSEUtils.everyOtherToHeadInPlace;

/*
 * Tests for the main methods in QuaSSEDistribution.
 * These tests were built after those in QuaSSEDistributionTest,
 * but expected outputs are sometimes the same.
 * Here we initialize QuaSSEDistribution instances, and let those
 * instances figure all necessary variables upon initialization.
 *
 * Note that expectations for D's in tests with PropagateCh(...)
 * come from R running on an 2014 iMac with Mac OS X Catalina.
 *
 * When initial D values are very small (because getNormalDensity returns
 * very tiny values as a result of both dx and nx being large), FFT outputs
 * will differ between machines with different CPU architectures. I tried
 * to avoid this by decreasing dx to 0.0005 in the 1024 and 4096 unit tests.
 *
 * These differences can be very large. After convolution,
 * the same D bin might have a positive value under one architecture, and a
 * negative value under another. These differences add up!
 */
public class QuaSSEDistributionTest {

    private final List<QuaSSEDistribution> opened = new ArrayList<>();

    // Run the same reference and MCMC-lifecycle checks on either FFT, without duplicating their data.
    // The native target selects FFT/integration implementations; reference tolerances stay unchanged.
    private QuaSSEDistribution newDistribution() {
        QuaSSEDistribution distribution = new QuaSSEDistribution();
        distribution.fftBackendInput.setValue(System.getProperty("test.quasse.fft", "sst"), distribution);
        distribution.integrationBackendInput.setValue(System.getProperty("test.quasse.integration", "java"), distribution);
        opened.add(distribution);
        return distribution;
    }

    // Release native resources even when an assertion or model initialization fails.
    @After
    public void closeDistributions() {
        for (QuaSSEDistribution distribution : opened) distribution.close();
        opened.clear();
    }

    final static Double EPSILON = 1e-6;

    static QuaSSEDistribution q1024, q32Dt002, q32Dt001, q32Dt0005, q32BifTree0025HeightDt0005, q32ThreeSpTreeDt0005, q32FifteenSp;
    static Tree bifTreeHeight0025, bifTreeHeight002, bifTreeHeight001, threeSpTreeHeight002, fifteenSpTree;

    static List<Double> data2Sp, data3Sp, data15Sp;
    static RealParameter driftRp, xMid00Rp, xMidFifteenSpRp, flankWidthScaler10Rp, flankWidthScaler5Rp;
    static RealParameter dxBin001Rp, dxBin00005Rp, dxBinFifteenSpRp;
    static RealParameter quTrait2SpRp, quTrait3SpRp, quTrait15SpRp;
    static RealParameter dt01Rp, dt001Rp, dt002Rp, dt0005Rp;
    static RealParameter tc100Rp, tc0005Rp, tcFifteenSpRp;
    static RealParameter diffusionRp0001, diffusionRp001;
    static IntegerParameter hiLoRatioIp, nXbins1024Ip, nXbins128Ip, nXbins48Ip, nXbins32Ip;
    static BooleanParameter dynDtbpTrue, dynDtbpFalse;
    static String rootPriorType;

    static double dt01, dt001, dt002, dt0005;
    static Double[] x0, y1, y0, r;
    static LogisticFunction lfn;
    static ConstantLinkFn cfn;
    static NormalCenteredAtObservedLinkFn nfn2Sp, nfn3Sp, nfn15Sp;

    @BeforeClass
    public static void setupParameters() {
        // tree
        String bifTreeStr001 = "(sp1:0.01,sp2:0.01);";
        String bifTreeStr002 = "(sp1:0.02,sp2:0.02);";
        String bifTreeStr0025 = "(sp1:0.025,sp2:0.025);";
        String trifTreeStr002 = "((sp2:0.01,sp3:0.01):0.01,sp1:0.02);";
        String fifteenTreeStr = "(sp2:13.77320255,(sp1:12.76688384,((((sp12:1.170387028,sp13:1.170387028)nd16:0.9837720325,sp9:2.154159061)nd11:5.451401092,((sp5:4.311645343,(sp14:0.8910055279,sp15:0.8910055279)nd14:3.420639815)nd9:2.536663776,((sp16:0.3011866125,sp17:0.3011866125)nd12:4.264383667,(sp6:3.95083843,sp7:3.95083843)nd13:0.6147318498)nd10:2.282738839)nd8:0.7572510339)nd5:2.554739141,((sp10:2.059478202,sp11:2.059478202)nd15:0.4198789018,sp8:2.479357104)nd6:7.68094219)nd4:2.60658455)nd3:1.006318707)nd1;";
        bifTreeHeight001 = new TreeParser(bifTreeStr001, false, false, true, 0);
        bifTreeHeight002 = new TreeParser(bifTreeStr002, false, false, true, 0);
        bifTreeHeight0025 = new TreeParser(bifTreeStr0025, false, false, true, 0);
        threeSpTreeHeight002 = new TreeParser(trifTreeStr002, false, false, true, 0);
        fifteenSpTree = new TreeParser(fifteenTreeStr, false, false, true, 0);

        // qu trait data
        String spNames2Sp = "sp1 sp2";
        String spNames3Sp = "sp1 sp2 sp3";
        String spNames15Sp = "sp1 sp2 sp5 sp6 sp7 sp8 sp9 sp10 sp11 sp12 sp13 sp14 sp15 sp16 sp17";
        data2Sp = Arrays.asList(0.0, 0.1);
        data3Sp = Arrays.asList(0.0, 0.1, 0.2);
        data15Sp = Arrays.asList(-0.05384594, -0.37091896, 0.59169195, 0.14947513, 0.46156791,
                0.27345680, 1.73358959, 0.38883347, 0.42233625, 1.55011787,
                1.17169681, 0.72422971, 0.84092251, 0.19645523, 0.48495092);
        quTrait2SpRp = new RealParameter();
        quTrait2SpRp.initByName("value", data2Sp, "keys", spNames2Sp);
        quTrait3SpRp = new RealParameter();
        quTrait3SpRp.initByName("value", data3Sp, "keys", spNames3Sp);
        quTrait15SpRp = new RealParameter();
        quTrait15SpRp.initByName("value", data15Sp, "keys", spNames15Sp);

        // qu trait stuff
        Double[] drift = new Double[] { 0.0 };
        driftRp = new RealParameter(drift);

        Double[] diffusion0001 = new Double[] { 0.001 };
        Double[] diffusion001 = new Double[] { 0.01 };
        diffusionRp0001 = new RealParameter(diffusion0001);
        diffusionRp001 = new RealParameter(diffusion001);

        // dimension stuff
        dt0005 = 0.005;
        dt001 = 0.01;
        dt002 = 0.02;
        dt01 = 0.1;
        Double[] dtD00005 = new Double[] { 0.005 };
        Double[] dtD001 = new Double[] { 0.01 };
        Double[] dtD002 = new Double[] { 0.02 };
        Double[] dtD01 = new Double[] { 0.1 };
        dt01Rp = new RealParameter(dtD01);
        dt001Rp = new RealParameter(dtD001);
        dt002Rp = new RealParameter(dtD002);
        dt0005Rp = new RealParameter(dtD00005);

        // adjust dt dynamically to maximize accuracy and match diversitree
        Boolean[] dynDtTrue = new Boolean[] { true };
        Boolean[] dynDtFalse = new Boolean[] { false };
        dynDtbpTrue = new BooleanParameter(dynDtTrue);
        dynDtbpFalse = new BooleanParameter(dynDtFalse);

        Double[] tc100 = new Double[] { 100.0 };
        Double[] tc0005 = new Double[] { 0.005 };
        // Double[] tcFifteenSp = new Double[] { 1.0 };
        Double[] tcFifteenSp = new Double[] { 1.37732 };
        tc100Rp = new RealParameter(tc100);
        tc0005Rp = new RealParameter(tc0005);
        tcFifteenSpRp = new RealParameter(tcFifteenSp);

        // logistic realparameter's for lambda
        x0 = new Double[] { 0.0 };
        y1 = new Double[] { 0.2 };
        y0 = new Double[] { 0.1 };
        r = new Double[] { 2.5 }; // logistic function parameters
        RealParameter y0rp = new RealParameter(y0);
        RealParameter y1rp = new RealParameter(y1);
        RealParameter x0rp = new RealParameter(x0);
        RealParameter rrp = new RealParameter(r);

        lfn = new LogisticFunction();
        lfn.initByName( "curveYBaseValue", y0rp, "curveMaxY", y1rp, "sigmoidMidpoint", x0rp, "logisticGrowthRate", rrp);

        // constant realparameter's for mu
        Double[] yValue = new Double[] { 0.03 };
        RealParameter yValuerp = new RealParameter(yValue);

        cfn = new ConstantLinkFn();
        cfn.initByName("yV", yValuerp);

        // link function for D's

        Double[] sdNormaQuTraitValue002 = new Double[] { 0.02 };
        Double[] sdNormaQuTraitValue005 = new Double[] { 0.05 };
        Double[] sdNormaQuTraitValue0005 = new Double[] { 0.005 };
        RealParameter sdNormaQuTraitValue002Rp = new RealParameter(sdNormaQuTraitValue002);
        RealParameter sdNormaQuTraitValue005Rp = new RealParameter(sdNormaQuTraitValue005);
        RealParameter sdNormaQuTraitValue0005Rp = new RealParameter(sdNormaQuTraitValue0005);
        nfn2Sp = new NormalCenteredAtObservedLinkFn();
        nfn2Sp.initByName("quTraits", quTrait2SpRp, "sdNormalQuTrValue", sdNormaQuTraitValue005Rp);
        nfn3Sp = new NormalCenteredAtObservedLinkFn();
        nfn3Sp.initByName("quTraits", quTrait3SpRp, "sdNormalQuTrValue", sdNormaQuTraitValue005Rp);
        nfn15Sp = new NormalCenteredAtObservedLinkFn();
        nfn15Sp.initByName("quTraits", quTrait15SpRp, "sdNormalQuTrValue", sdNormaQuTraitValue002Rp);

        // nfn2 = new NormalCenteredAtObservedLinkFn();
        // nfn2.initByName("quTraits", quTraitRp, "sdNormalQuTrValue", sdNormaQuTraitValuerp);

        Double[] xMid00 = new Double[] { 0.0 };
        // Double[] xMidFifteenSp = new Double[] { 0.0 };
        Double[] xMidFifteenSp = new Double[] { 0.6813353 };
        xMid00Rp = new RealParameter(xMid00);
        xMidFifteenSpRp = new RealParameter(xMidFifteenSp);

        Double[] dxBin001 = new Double[] { 0.01 };
        dxBin001Rp = new RealParameter(dxBin001);

        Double[] dxBin00005 = new Double[] { 0.0005 };
        dxBin00005Rp = new RealParameter(dxBin00005);

//        Double[] dxBinFifteenSpD = new Double[] { 0.1 };
//        dxBinFifteenSpRp = new RealParameter(dxBinFifteenSpD);
        Double[] dxBinFifteenSpD = new Double[] { 0.01027592 };
        dxBinFifteenSpRp = new RealParameter(dxBinFifteenSpD);

        Double[] flankWidthScaler10 = new Double[] { 10.0 };
        Double[] flankWidthScaler5 = new Double[] { 5.0 };
        flankWidthScaler10Rp = new RealParameter(flankWidthScaler10);
        flankWidthScaler5Rp = new RealParameter(flankWidthScaler5);

        Integer[] hiLoRatio = new Integer[] { 4 };
        hiLoRatioIp = new IntegerParameter(hiLoRatio);

        Integer[] nXbins1024 = new Integer[] { 1024 };
        nXbins1024Ip = new IntegerParameter(nXbins1024);

        Integer[] nXbins128 = new Integer[] { 128 };
        nXbins128Ip = new IntegerParameter(nXbins128);

        Integer[] nXbins48 = new Integer[] { 48 };
        nXbins48Ip = new IntegerParameter(nXbins48);

        Integer[] nXbins32 = new Integer[] { 32 };
        nXbins32Ip = new IntegerParameter(nXbins32);

        rootPriorType = "Observed";
    }

    /*
     * Setup before each test (initializing from scratch every time,
     * otherwise we modify the state of classes in place and a test
     * interferes with the next one.
     */
    @Before
    public void setupQuaSSELiks() {
        q1024 = newDistribution();
        q1024.initByName("dtMax", dt001Rp, "dynDt", dynDtbpTrue,
                "tc", tc100Rp,
                "nX", nXbins1024Ip, "dX", dxBin00005Rp, "xMid", xMid00Rp, "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight002,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);

        q32Dt001 = newDistribution();
        q32Dt001.initByName("dtMax", dt001Rp, "dynDt", dynDtbpTrue,
                "tc", tc100Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp, "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight001,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);

        q32Dt002 = newDistribution();
        q32Dt002.initByName("dtMax", dt002Rp, "dynDt", dynDtbpTrue,
                "tc", tc100Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp, "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight002,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);

        q32Dt0005 = newDistribution();
        q32Dt0005.initByName("dtMax", dt0005Rp, "dynDt", dynDtbpTrue,
                "tc", tc0005Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp, "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight001,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);

        q32BifTree0025HeightDt0005 = newDistribution();
        q32BifTree0025HeightDt0005.initByName("dtMax", dt0005Rp, "dynDt", dynDtbpTrue,
                "tc", tc0005Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp, "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight0025,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);

        q32ThreeSpTreeDt0005 = newDistribution();
        q32ThreeSpTreeDt0005.initByName("dtMax", dt0005Rp, "dynDt", dynDtbpTrue,
                "tc", tc0005Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp, "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", threeSpTreeHeight002,
                "q2d", nfn3Sp,
                "priorProbAtRootType", rootPriorType);

        q32FifteenSp = newDistribution();
        q32FifteenSp.initByName("dtMax", dt0005Rp, "dynDt", dynDtbpTrue,
                "tc", tcFifteenSpRp,
                "nX", nXbins1024Ip, "dX", dxBinFifteenSpRp, "xMid", xMidFifteenSpRp, "flankWidthScaler", flankWidthScaler5Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", fifteenSpTree,
                "q2d", nfn15Sp,
                "priorProbAtRootType", rootPriorType);
    }

    /*
     * Checks that QuaSSE likelihood class it setting its own
     * dimensions correctly, and then
     *
     * (1) Applying a logistic function to convert the quantitative trait
     * into the birth rate (lambda)
     * (2) Applying a constant function to convert the quantitative trait
     * into the death rate (mu)
     */
    @Test
    public void testDimensions() {
        int nXbinsLo = q1024.getnXbins(true);
        int nXbinsHi = q1024.getnXbins(false);
        int nUsefulXbinsLo = q1024.getNUsefulXbins(true);
        int nUsefulXbinsHi = q1024.getNUsefulXbins(false);
        int nLeftFlanksLo = q1024.getNLeftFlanks(true);
        int nLeftFlanksHi = q1024.getNLeftFlanks(false);
        int nRightFlanksLo = q1024.getNRightFlanks(true);
        int nRightFlanksHi = q1024.getNRightFlanks(false);
        double xMinLo = q1024.getXMinLo();
        double xMinHi = q1024.getXMinHi();
        double[] xLo = q1024.getX(true);
        double[] xHi = q1024.getX(false);
        double[] lambdaLo = q1024.getLambda(true);
        double[] lambdaHi = q1024.getLambda(false);
        double[] muLo = q1024.getMu(true);
        double[] muHi = q1024.getMu(false);

        // expected x rulers
        double[] expectedXLoFirst10 = new double[] { -0.2235, -0.223, -0.2225, -0.222, -0.2215, -0.221, -0.2205, -0.22, -0.2195, -0.219 };
        double[] expectedXLoLast10 = new double[] { 0.219, 0.2195, 0.22, 0.2205, 0.221, 0.2215, 0.222, 0.2225, 0.223, 0.2235 };
        double[] expectedXHiFirst10 = new double[] { -0.223875, -0.22375, -0.223625, -0.2235, -0.223375, -0.22325, -0.223125, -0.223, -0.222875, -0.22275 };
        double[] expectedXHiLast10 = new double[] { 0.22275, 0.222875, 0.223, 0.223125, 0.22325, 0.223375, 0.2235, 0.223625, 0.22375, 0.223875 };

        // expected lambdas
        double[] expectedLambdaLoFirst10 = new double[] { 0.13638367349015, 0.136412610857295, 0.13644155805569, 0.136470515067719, 0.136499481875741, 0.136528458462084, 0.136557444809052, 0.13658644089892, 0.136615446713936, 0.136644462236322 };
        double[] expectedLambdaLoLast10 = new double[] { 0.163355537763678, 0.163384553286064, 0.16341355910108, 0.163442555190948, 0.163471541537916, 0.163500518124259, 0.163529484932281, 0.16355844194431, 0.163587389142705, 0.16361632650985 };
        double[] expectedLambdaHiFirst10 = new double[] { 0.136361976927131, 0.136369208498794, 0.136376440686558, 0.13638367349015, 0.136390906909294, 0.136398140943716, 0.136405375593141, 0.136412610857295, 0.136419846735901, 0.136427083228686 };
        double[] expectedLambdaHiLast10 = new double[] { 0.163572916771314, 0.163580153264099, 0.163587389142705, 0.163594624406859, 0.163601859056284, 0.163609093090706, 0.16361632650985, 0.163623559313442, 0.163630791501206, 0.163638023072869 };

        // expected mus for both lo and hi, first and last 10
        double[] expectedMuLoHiFirstLast10 = new double[] { 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03 };

        Assert.assertEquals(nXbinsLo, 1024, 0.0);
        Assert.assertEquals(nXbinsHi, 1024 * 4, 0.0);
        Assert.assertEquals(895, nUsefulXbinsLo, 0.0);
        Assert.assertEquals(3583, nUsefulXbinsHi, 0.0);
        Assert.assertEquals(64, nLeftFlanksLo, 0.0);
        Assert.assertEquals(256, nLeftFlanksHi, 0.0);
        Assert.assertEquals(64, nRightFlanksLo, 0.0);
        Assert.assertEquals(256, nRightFlanksHi, 0.0);
        Assert.assertEquals(-0.2235, xMinLo, EPSILON);
        Assert.assertEquals(-0.223875, xMinHi, EPSILON);
        Assert.assertEquals(895, xLo.length, 0.0);
        Assert.assertEquals(3583, xHi.length, 0.0);
        Assert.assertArrayEquals(expectedXLoFirst10, Arrays.copyOfRange(xLo, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedXLoLast10, Arrays.copyOfRange(xLo, 885, 895), EPSILON);
        Assert.assertArrayEquals(expectedXHiFirst10, Arrays.copyOfRange(xHi, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedXHiLast10, Arrays.copyOfRange(xHi, 3573, 3583), EPSILON);
        Assert.assertEquals(895, lambdaLo.length, 0.0);
        Assert.assertEquals(3583, lambdaHi.length, 0.0);
        Assert.assertArrayEquals(expectedLambdaLoFirst10, Arrays.copyOfRange(lambdaLo, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedLambdaLoLast10, Arrays.copyOfRange(lambdaLo, 885, 895), EPSILON);
        Assert.assertArrayEquals(expectedLambdaHiFirst10, Arrays.copyOfRange(lambdaHi, 0, 10), EPSILON);
        Assert.assertArrayEquals(expectedLambdaHiLast10, Arrays.copyOfRange(lambdaHi, 3573, 3583), EPSILON);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muLo, 0, 10), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muLo, 885, 895), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muHi, 0, 10), 0.0);
        Assert.assertArrayEquals(expectedMuLoHiFirstLast10, Arrays.copyOfRange(muHi, 3573, 3583), 0.0);
    }

    /*
     * Trait-dependent extinction must be evaluated on both x rulers, just like speciation. This
     * supplements the constant-μ checks and becomes redundant if link evaluation moves out of QuaSSE.
     */
    @Test
    public void testTraitDependentExtinctionRates() {
        QuaSSEDistribution distribution = newDistribution();
        distribution.initByName("dtMax", dt001Rp, "dynDt", dynDtbpTrue,
                "tc", tc100Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp,
                "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", lfn,
                "tree", bifTreeHeight001,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);

        Assert.assertArrayEquals(distribution.getLambda(true), distribution.getMu(true), 0.0);
        Assert.assertArrayEquals(distribution.getLambda(false), distribution.getMu(false), 0.0);
    }

    /*
     * Positive and negative drift require mirrored padding and backward-time kernel shifts. Zero-drift
     * tests cannot expose duplicated flank formulas; this remains useful for the constant-drift kernel.
     */
    @Test
    public void testNonzeroDriftPaddingAndKernelDirection() {
        RealParameter positiveDrift = new RealParameter(new Double[] { 1.0 });
        RealParameter negativeDrift = new RealParameter(new Double[] { -1.0 });
        QuaSSEDistribution positive = newDistribution();
        positive.initByName("dtMax", dt001Rp, "dynDt", dynDtbpTrue,
                "tc", tc100Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp,
                "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", positiveDrift, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight001,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);
        QuaSSEDistribution negative = newDistribution();
        negative.initByName("dtMax", dt001Rp, "dynDt", dynDtbpTrue,
                "tc", tc100Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp,
                "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", negativeDrift, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight001,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);

        Assert.assertEquals(5, positive.getNLeftFlanks(true));
        Assert.assertEquals(3, positive.getNRightFlanks(true));
        Assert.assertEquals(3, negative.getNLeftFlanks(true));
        Assert.assertEquals(5, negative.getNRightFlanks(true));
        Assert.assertEquals(20, positive.getNLeftFlanks(false));
        Assert.assertEquals(12, positive.getNRightFlanks(false));
        Assert.assertEquals(12, negative.getNLeftFlanks(false));
        Assert.assertEquals(20, negative.getNRightFlanks(false));

        positive.populatefY(dt001, true, false, true, false);
        negative.populatefY(dt001, true, false, true, false);
        double[] positiveKernel = positive.getfY(true);
        double[] negativeKernel = negative.getfY(true);
        for (int i = 0; i < nXbins32Ip.getValue(); i++) {
            int mirrored = (nXbins32Ip.getValue() - i) % nXbins32Ip.getValue();
            Assert.assertEquals(positiveKernel[2 * i], negativeKernel[2 * mirrored], 1e-15);
        }
        Assert.assertTrue(positiveKernel[2 * 31] > positiveKernel[0]);
        Assert.assertTrue(negativeKernel[2] > negativeKernel[0]);
    }

    /*
     * Checks that two species in a 2sp-tree get their D's correctly
     * initialized (using standard QuaSSE initialization, in high resolution).
     * "sp1" and "sp2" are assigned trait values of 0.0 and 0.01, respectively.
     */
    @Test
    public void testInitializationOfTips() {
        double[][][] esDsHi = q1024.getEsDs(false);

        double[] expectedSp1Ds = new double[] { 6.96077358436849, 6.95166528584696, 6.94252551477923, 6.93335442671583, 6.92415217758908, 6.91491892370871, 6.90565482175746, 6.89636002878667, 6.88703470221182, 6.87767899980818, 6.86829307970631, 6.85887710038768, 6.84943122068018, 6.83995559975374, 6.83045039711584, 6.82091577260705, 6.81135188639661, 6.80175889897795, 6.79213697116422, 6.78248626408384, 6.772806939176, 6.76309915818623, 6.75336308316186, 6.74359887644761, 6.73380670068105, 6.72398671878815, 6.71413909397875, 6.70426398974212, 6.69436156984244, 6.68443199831428, 6.67447543945815, 6.66449205783599, 6.65448201826663, 6.64444548582134, 6.63438262581928, 6.62429360382306, 6.61417858563415, 6.60403773728847, 6.59387122505179, 6.58367921541529, 6.57346187509105, 6.5632193710075, 6.55295187030495, 6.54265954033109, 6.53234254863645, 6.52200106296994, 6.51163525127429, 6.50124528168164 };
        double[] expectedSp2Ds = new double[] { 2.67859021074856, 2.68849414736158, 2.69841783805887, 2.70836123148143, 2.71832427571076, 2.72830691826807, 2.73830910611356, 2.74833078564564, 2.75837190270027, 2.76843240255022, 2.77851222990441, 2.78861132890721, 2.79872964313779, 2.80886711560949, 2.81902368876922, 2.82919930449678, 2.83939390410431, 2.84960742833573, 2.85983981736613, 2.87009101080125, 2.88036094767694, 2.89064956645866, 2.90095680504094, 2.91128260074695, 2.92162689032799, 2.93198960996305, 2.9423706952584, 2.95277008124712, 2.96318770238874, 2.97362349256884, 2.9840773850987, 2.9945493127149, 3.00503920757903, 3.01554700127736, 3.02607262482052, 3.03661600864323, 3.04717708260405, 3.05775577598507, 3.06835201749176, 3.07896573525268, 3.0895968568193, 3.10024530916586, 3.11091101868917, 3.12159391120842, 3.13229391196515, 3.14301094562307, 3.15374493626797, 3.16449580740766 };

        double[] sp1DsSST = new double[48];
        double[] sp2DsSST = new double[48];

        // TODO: later I need to update this and instead of making sp1DsSST look like the expected,
        // I should change the expected to alternate between values and zeros
        int j = 0;
        for (int i=0; i<48*2; i+=2) {
            sp1DsSST[j] = esDsHi[0][1][4000 + i];
            sp2DsSST[j] = esDsHi[1][1][4000 + i];
            j += 1;
        }

        Assert.assertArrayEquals(expectedSp1Ds, sp1DsSST, 1E-12);
        Assert.assertArrayEquals(expectedSp2Ds, sp2DsSST, 1E-12);

        // if jtransforms flag in QuaSSEDistribution initialization is set to true
        // Assert.assertArrayEquals(expectedSp1Ds, Arrays.copyOfRange(esDsHi[0][1], 2000, 2048), 1E-12);
        // Assert.assertArrayEquals(expectedSp2Ds, Arrays.copyOfRange(esDsHi[1][1], 2000, 2048), 1E-12);
    }

    /*
     * Checks that propagate methods in time for E's and D's
     * inside QuaSSE class are working.
     *
     * Differs from 'testPropagateTimeOneChQuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 0.01,
     * and nXbins = 32 (low res).
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0.
     *
     * This test, as is, uses SST's JavaFftService, which requires that
     * esDs arrays have elements in every other index; we can adjust this
     * test to not do that, which is what would be required to use JTransforms
     */
    @Test
    public void testIntegrateOneBranchLoRes32BinsOutsideClassJustT() {

        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are used in propagate in x
         */
        double[][] esDsHiAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D

        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];

        // just propagate in t, in place
        boolean jtransforms = false;
        q32Dt001.propagateTInPlace(esDsLoAtNode, scratchAtNode, dt001, true, jtransforms);

        // if jtransforms flag in QuaSSEDistribution initialization is set to true
        // boolean jtransforms = true;
        // q32Dt001.propagateTInPlace(esDsLoAtNode, scratchAtNode, dt001, true, jtransforms);

        esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);
        double[] esHiAtNode = esDsLoAtNode[0];
        double[] dsHiAtNode = esDsLoAtNode[1];

        double[] expectedInitialDs = new double[] { 0.709491856924629, 0, 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.9945493127149, 0, 3.88372109966426, 0, 4.83941449038287, 0, 5.79383105522966, 0, 6.66449205783599, 0, 7.36540280606647, 0, 7.82085387950912, 0, 7.97884560802865, 0, 7.82085387950912, 0, 7.36540280606647, 0, 6.66449205783599, 0, 5.79383105522966, 0, 4.83941449038287, 0, 3.88372109966426, 0, 2.9945493127149, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1EsAfterPropT = new double[] { 0.000299740440744498, 0, 0.000299739520470888, 0, 0.000299738597335782, 0, 0.00029973767161558, 0, 0.000299736743589792, 0, 0.000299735813540893, 0, 0.000299734881753716, 0, 0.000299733948515344, 0, 0.00029973301411462, 0, 0.00029973207884177, 0, 0.000299731142988294, 0, 0.000299730206846208, 0, 0.00029972927070804, 0, 0.000299728334866242, 0, 0.000299727399612858, 0, 0.000299726465239364, 0, 0.000299725532035927, 0, 0.000299724600291355, 0, 0.000299723670292635, 0, 0.000299722742324704, 0, 0.000299721816669763, 0, 0.000299720893607378, 0, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropT = new double[] { 0.708264611249434, 0, 1.07794488915543, 0, 1.57625248872262, 0, 2.21453845459856, 0, 2.98929572353432, 0, 3.87688349797518, 0, 4.83086427268124, 0, 5.78355856417974, 0, 6.65263439389469, 0, 7.35225216820117, 0, 7.80684129110362, 0, 7.96450018578724, 0, 7.80674374052118, 0, 7.35206845754541, 0, 6.65238511637526, 0, 5.78326972079126, 0, 4.83056283595408, 0, 3.87659337350287, 0, 2.98903491521347, 0, 2.21431781312691, 0, 1.57607596745573, 0, 1.07781089197137, 0, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esHiAtNode, 0, 64), 1E-9);
        Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(dsHiAtNode, 0, 64), 1E-14);

        // if jtransforms flag in QuaSSEDistribution initialization is set to true
        // double[] expectedInitialDs = new double[] { 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1EsAfterPropT = new double[] { 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1DsAfterPropT = new double[] { 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052118, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esHiAtNode, 0, 32), 1E-9);
        // Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(dsHiAtNode, 0, 32), 1E-14);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Uses JTransforms.
     *
     * Differs from 'testPropagateChOneCh1024QuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 0.01,
     * and nXbins = 32 (low res).
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes32BinsOutsideClassJustXJTransforms() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are used in propagate in x
         */
        double[][] esDsHiAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];

        /*
         * QuaSSEDistribution.java has jtransforms=false when initializing E's and D's,
         * which means E's and D's will be alternated with 0.0 values; in order for the
         * current unit test to work, we need to compress all values toward the head of
         * the array, removing the 0.0's
         */
        everyOtherToHeadInPlace(esDsLoAtNode[0], 32, 0, 0, 2, 1.0); // grabbing real part and scaling by 1/nXbins
        everyOtherToHeadInPlace(esDsLoAtNode[1], 32, 0, 0, 2, 1.0);

        // now we populate the arrays we will use in the assertions
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D

        /*
         * fY is computed and FFTed in initialization, by QuaSSEProcess;
         * here we are just grabbing it to verify its values in the asserts
         * below
         */
        q32Dt001.populatefY(0.01, true, true, true, true);
        double[] fftedfY = q32Dt001.getfY(true);

        // copying fY for assert (leaving original one inside class untouched)
        double[] fftedfY4Assert = new double[fftedfY.length]; // just for test, not used in propagate in X
        for (int i=0; i<fftedfY.length; i++) fftedfY4Assert[i] = fftedfY[i];
        everyOtherToHeadInPlace(fftedfY4Assert, q32Dt001.getnXbins(true),0, 0, 2, 1.0); // getting real part for assert below

        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];
        q32Dt001.propagateXInPlaceJTransforms(esDsLoAtNode, scratchAtNode, true);

        // looking at 'sp1'
        esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);
        double[] esLoAtNode = esDsLoAtNode[0];
        double[] dsLoAtNode = esDsLoAtNode[1];

        double[] expectedInitialDs = new double[] { 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedFFTedfY = new double[] { 1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216, 0.981614853950298, 0.979316029808302, 0.977301093995625, 0.975647479188405, 0.97441873269917, 0.973662074416229, 0.973406582192705, 0.973662074416229, 0.97441873269917, 0.975647479188405, 0.977301093995625, 0.979316029808302, 0.981614853950298, 0.984109224049216, 0.986703282961364, 0.989297342492751, 0.991791714355113, 0.994090541136289, 0.996105480062091, 0.997759097982437, 0.998987847110852, 0.999744507157237 };
        double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropX = new double[] { 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.99530083804265, 3.88416335936269, 4.83940600155165, 5.79327421787607, 6.66336349661662, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119625, 6.66336349661661, 5.79327421787607, 4.83940600155166, 3.88416335936269, 2.99530083804265, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 32), 1e-14);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 32), 1e-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 32), 1e-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 32), 1e-14);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Uses SST.
     *
     * Differs from 'testPropagateChOneCh1024QuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree over a single dt = 0.01,
     * and nXbins = 32 (low res).
     *
     * We look at just a single branch here, from 'sp1', whose trait value
     * is set to 0.0.
     */
    @Test
    public void testIntegrateOneBranchLoRes32BinsOutsideClassJustX() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);
        double[][] fftBufferEsDsAtNode = new double[2][esDsLoAtNode[0].length];

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are used in propagate in x
         */
        double[][] esDsHiAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D

        // SSEUtils.everyOtherToHeadInPlace(esDsHiAtNodeInitial[1], q32Dt001.getnXbins(true), 0, 0, 2, 1.0);

        /*
         * fY is computed and FFTed in initialization, by QuaSSEProcess;
         * here we are just grabbing it to verify its values in the asserts
         * below
         */
        q32Dt001.populatefY(0.01, true, true, true, false);
        double[] fftedfY = q32Dt001.getfftFY(true);

        // copying fY for assert (leaving original one inside class untouched)
        double[] fftedfY4Assert = new double[fftedfY.length]; // just for test, not used in propagate in X
        for (int i=0; i<fftedfY.length; i++) fftedfY4Assert[i] = fftedfY[i];

        // everyOtherToHeadInPlace(fftedfY4Assert, q32Dt001.getnXbins(true),0, 0, 2, 1.0); // getting real part for assert below

        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];
        q32Dt001.propagateXInPlace(esDsLoAtNode, fftBufferEsDsAtNode, scratchAtNode, true);

        // looking at 'sp1'
        esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);
        double[] esLoAtNode = esDsLoAtNode[0];
        double[] dsLoAtNode = esDsLoAtNode[1];

        double[] expectedInitialDs = new double[] { 0.709491856924629, 0, 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.9945493127149, 0, 3.88372109966426, 0, 4.83941449038287, 0, 5.79383105522966, 0, 6.66449205783599, 0, 7.36540280606647, 0, 7.82085387950912, 0, 7.97884560802865, 0, 7.82085387950912, 0, 7.36540280606647, 0, 6.66449205783599, 0, 5.79383105522966, 0, 4.83941449038287, 0, 3.88372109966426, 0, 2.9945493127149, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedFFTedfY = new double[] { 1, 0, 0.999744507157237, 0, 0.998987847110852, 0, 0.997759097982437, 0, 0.996105480062091, 0, 0.994090541136289, 0, 0.991791714355113, 0, 0.989297342492751, 0, 0.986703282961364, 0, 0.984109224049216, 0, 0.981614853950298, 0, 0.979316029808302, 0, 0.977301093995625, 0, 0.975647479188405, 0, 0.97441873269917, 0, 0.973662074416229, 0, 0.973406582192705, 0, 0.973662074416229, 0, 0.97441873269917, 0, 0.975647479188405, 0, 0.977301093995625, 0, 0.979316029808302, 0, 0.981614853950298, 0, 0.984109224049216, 0, 0.986703282961364, 0, 0.989297342492751, 0, 0.991791714355113, 0, 0.994090541136289, 0, 0.996105480062091, 0, 0.997759097982437, 0, 0.998987847110852, 0, 0.999744507157237, 0 };
        double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropX = new double[] { 0.709491856924629, 0, 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.99530083804265, 0, 3.88416335936269, 0, 4.83940600155165, 0, 5.79327421787607, 0, 6.66336349661662, 0, 7.36377090119626, 0, 7.81887626199731, 0, 7.97674483551017, 0, 7.81887626199731, 0, 7.36377090119625, 0, 6.66336349661661, 0, 5.79327421787607, 0, 4.83940600155166, 0, 3.88416335936269, 0, 2.99530083804265, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedInitialDs = new double[] { 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedFFTedfY = new double[] { 1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216, 0.981614853950298, 0.979316029808302, 0.977301093995625, 0.975647479188405, 0.97441873269917, 0.973662074416229, 0.973406582192705, 0.973662074416229, 0.97441873269917, 0.975647479188405, 0.977301093995625, 0.979316029808302, 0.981614853950298, 0.984109224049216, 0.986703282961364, 0.989297342492751, 0.991791714355113, 0.994090541136289, 0.996105480062091, 0.997759097982437, 0.998987847110852, 0.999744507157237 };
        // double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1DsAfterPropX = new double[] { 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.99530083804265, 3.88416335936269, 4.83940600155165, 5.79327421787607, 6.66336349661662, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119625, 6.66336349661661, 5.79327421787607, 4.83940600155166, 3.88416335936269, 2.99530083804265, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        // debugging
//        DecimalFormat df = new DecimalFormat("##.########");
//        df.setRoundingMode(RoundingMode.DOWN);
//        System.out.println("    D's        expected D's");
//        for (int i=0; i<dsLoAtNode.length; i++) {
//            if (dsLoAtNode[i] <= 1e-14) {
//                String expVal = "";
//                if (!(expectedSp1DsAfterPropX[i] > 0) || !(expectedSp1DsAfterPropX[i] < 0))  {
//                    expVal = "0.00000000";
//                }
//                System.out.println("i=" + (i+1) + " 0.00000000" + " " + expVal);
//            } else {
//                System.out.println("i=" + (i + 1) + " " + df.format(dsLoAtNode[i]) + " " + df.format(expectedSp1DsAfterPropX[i]));
//            }
//        }

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 64), 1e-14);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 64), 1e-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 64), 1e-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 64), 1e-14);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Differs from 'testPropagateChOneCh1024QuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree (but just looking at a single
     * branch here, "sp1"), over a single dt = 0.01, with dx = 0.0005,
     * and nXbins = 4096 (high res). The trait value of "sp1" is set to 0.0.
     *
     * The method for FFTing used here is JavaFftService, so arrays
     * have values spread over while ignoring even indices
     */
    @Test
    public void testIntegrateOneBranchLoRes1024BinsOutsideClassJustX() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        // we'll test the integration outside the class
        esDsLoAtNode = q1024.getEsDsAtNode(nodeIdx, true);
        double[][] fftBufferEsDsAtNode = new double[2][esDsLoAtNode[0].length];

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are used in propagate in x
         */
        double[][] esDsHiAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D

        /*
         * fY is computed and FFTed in initialization, by QuaSSEProcess;
         * here we are just grabbing it to verify its values in the asserts
         * below
         */
        double aDt = 0.01;
        q1024.populatefY(aDt, true, true, true, false);
        double[] fftedfY = q1024.getfftFY(true);

        // copying fY for assert (leaving original one inside class untouched)
        double[] fftedfY4Assert = new double[fftedfY.length]; // just for test, not used in propagate in X
        for (int i=0; i<fftedfY.length; i++) fftedfY4Assert[i] = fftedfY[i];

        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];
        q1024.propagateXInPlace(esDsLoAtNode, fftBufferEsDsAtNode, scratchAtNode, true);

        // if using JTransforms
        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        // double[][] scratchAtNode = q1024.getScratchAtNode(nodeIdx, true); // sp1
        // q1024.propagateXInPlaceJTransforms(esDsLoAtNode, scratchAtNode, true);

        esDsLoAtNode = q1024.getEsDsAtNode(nodeIdx, true);
        double[] esLoAtNode = esDsLoAtNode[0];
        double[] dsLoAtNode = esDsLoAtNode[1];

        // System.out.println("Final D's: " + Arrays.toString(dsLoAtNode));

        double[] expectedInitialDs = new double[] { 0.000365714984190948, 0, 0.000382414193856355, 0, 0.000399835934138456, 0, 0.00041800955800901, 0, 0.000436965523466329, 0, 0.000456735431302938, 0, 0.000477352064003592, 0, 0.000498849425801072, 0, 0.000521262783917567, 0,0.000544628711019852,0,0.000568985128916886,0,0.000594371353528846,0,0.000620828141157005,0,0.000648397736084276,0,0.000677123919536558,0,0.000707052060035464,0,0.000738229165173324,0,0.000770703934841743,0,0.000804526815945299,0,0.000839750058632349,0,0.000876427774075162,0,0.000914615993832026,0,0.000954372730824099,0,0.000995758041960245,0,0.0010388340924432,0,0.0010836652217908,0,0.00113031801160615,0,0.0011788613551308,0,0.00122936652861539,0,0.00128190726454212,0,0.00133655982673381,0,0.00139340308738429,0,0.00145251860604505,0,0.00151399071060322,0,0.00157790658028586,0,0.00164435633072572,0,0.00171343310112364,0,0.00178523314354266,0,0.00185985591436892,0,0.00193740416797439,0,0.00201798405261629,0,0.00210170520860801,0,0.00218868086879601,0,0.00227902796137729,0,0.00237286721509132,0,0.00247032326682047,0,0.00257152477163242,0,0.00267660451529771,0 };
        double[] expectedFFTedfY = new double[] { 1, 0, 0.999247292368191, 0, 0.996992567179915, 0, 0.993245992007481, 0, 0.988024427909734, 0, 0.981351303026827, 0, 0.973256437474405, 0, 0.963775821368549, 0, 0.952951348293497, 0, 0.940830506973933, 0, 0.92746603432656, 0, 0.912915533436717, 0, 0.897241060330147, 0, 0.880508683684162, 0, 0.862788021843104, 0, 0.844151761668242, 0, 0.824675163860569, 0, 0.804435559446082, 0, 0.783511842107391, 0, 0.761983960984176, 0, 0.739932418450195, 0, 0.717437777209001, 0, 0.694580180837756, 0, 0.671438891652661, 0, 0.64809184947515, 0, 0.624615254550175, 0, 0.601083177512084, 0, 0.577567198915383, 0, 0.554136080452835, 0, 0.530855469577788, 0, 0.507787638837061, 0, 0.484991260810779, 0, 0.462521219151724, 0, 0.440428455824044, 0, 0.418759854264325, 0, 0.39755815783139, 0, 0.376861922578424, 0, 0.356705503075537, 0, 0.3371190697353, 0, 0.318128655850257, 0, 0.299756232341542, 0, 0.282019808042489, 0, 0.264933553200873, 0, 0.248507943778199, 0, 0.232749924053504, 0, 0.217663085001486, 0, 0.203247855908884, 0, 0.189501706717032, 0 };
        double[] expectedSp1EsFirst48AfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsFirst48AfterPropX = new double[] { 0.000365714984190948, 0, 0.000382414193856355, 0, 0.000399835934138456, 0, 0.00041800955800901, 0, 0.000436965523466329, 0, 0.000456735431302938, 0, 0.000477352064003592, 0, 0.000498849425801072, 0, 0.000521262783917567, 0, 0.000544628711019852, 0, 0.000568985128916886, 0, 0.000594371353528846, 0, 0.000620828141157005, 0, 0.000648397736084276, 0, 0.000677123919536558, 0, 0.000707052060035464, 0, 0.000738229165173324, 0, 0.000770703934841743, 0, 0.000804526815945299, 0, 0.000839750058632349, 0, 0.000876427774075162, 0, 0.000914615993832026, 0, 0.000954372730824099, 0, 0.000995758041960245, 0, 0.0010388340924432, 0, 0.0010836652217908, 0, 0.00113031801160615, 0, 0.0011788613551308, 0, 0.00122936652861539, 0, 0.00128190726454212, 0, 0.00133655982673381, 0, 0.00139340308738429, 0, 0.00145251860604505, 0, 0.00151399071060322, 0, 0.00157790658028586, 0, 0.00164435633072572, 0, 0.00171343310112364, 0, 0.00178523314354266, 0, 0.00185985591436892, 0, 0.00193740416797439, 0, 0.00201798405261629, 0, 0.00210170520860801, 0, 0.00218868086879601, 0, 0.00227902796137729, 0, 0.00237286721509132, 0, 0.00247032326682047, 0, 0.00257152477163242, 0, 0.00267660451529771, 0 };
        double[] expectedSp1DsLater48AfterPropX = new double[] { 0.00267660451529771, 0, 0.00257152477163242, 0, 0.00247032326682047, 0, 0.00237286721509132, 0, 0.0022790279613773, 0, 0.00218868086879601, 0, 0.00210170520860801, 0, 0.0020179840526163, 0, 0.00193740416797439, 0, 0.00185985591436892, 0, 0.00178523314354266, 0, 0.00171343310112364, 0, 0.00164435633072573, 0, 0.00157790658028586, 0, 0.00151399071060322, 0, 0.00145251860604505, 0, 0.00139340308738429, 0, 0.00133655982673381, 0, 0.00128190726454212, 0, 0.00122936652861539, 0, 0.0011788613551308, 0, 0.00113031801160615, 0, 0.0010836652217908, 0, 0.0010388340924432, 0, 0.000995758041960245, 0, 0.000954372730824099, 0, 0.000914615993832026, 0, 0.000876427774075162, 0, 0.000839750058632349, 0, 0.000804526815945299, 0, 0.000770703934841743, 0, 0.000738229165173324, 0, 0.000707052060035464, 0, 0.000677123919536558, 0, 0.000648397736084276, 0, 0.000620828141157005, 0, 0.000594371353528846, 0, 0.000568985128916886, 0, 0.000544628711019852, 0, 0.000521262783917567, 0, 0.000498849425801072, 0, 0.000477352064003592, 0, 0.000456735431302938, 0, 0.000436965523466329, 0, 0.00041800955800901, 0, 0.000399835934138456, 0, 0.000382414193856355, 0, 0.000365714984190948, 0 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 96), 1E-17);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 96), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsFirst48AfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 96), 1E-16);
        Assert.assertArrayEquals(expectedSp1DsFirst48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 96), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsLater48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 1694, 1790), 1E-14);

        // if using JTranforms
        // double[] expectedInitialDs = new double[] { 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771 };
        // double[] expectedFFTedfY = new double[] { 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368548, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577788, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032 };
        // double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1DsFirst48AfterPropX = new double[] { 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771 };
        // double[] expectedSp1DsLast48AfterPropX = new double[] { 0.00267660451529771, 0.00257152477163242, 0.00247032326682047, 0.00237286721509132, 0.0022790279613773, 0.00218868086879601, 0.00210170520860801, 0.0020179840526163, 0.00193740416797439, 0.00185985591436892, 0.00178523314354266, 0.00171343310112364, 0.00164435633072573, 0.00157790658028586, 0.00151399071060322, 0.00145251860604505, 0.00139340308738429, 0.00133655982673381, 0.00128190726454212, 0.00122936652861539, 0.0011788613551308, 0.00113031801160615, 0.0010836652217908, 0.0010388340924432, 0.000995758041960245, 0.000954372730824099, 0.000914615993832026, 0.000876427774075162, 0.000839750058632349, 0.000804526815945299, 0.000770703934841743, 0.000738229165173324, 0.000707052060035464, 0.000677123919536558, 0.000648397736084276, 0.000620828141157005, 0.000594371353528846, 0.000568985128916886, 0.000544628711019852, 0.000521262783917567, 0.000498849425801072, 0.000477352064003592, 0.000456735431302938, 0.000436965523466329, 0.00041800955800901, 0.000399835934138456, 0.000382414193856355, 0.000365714984190948 };
        // Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 48), 1E-17);
        // Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 48), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 48), 1E-16);
        // Assert.assertArrayEquals(expectedSp1DsFirst48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 48), 1E-17);
        // Assert.assertArrayEquals(expectedSp1DsLast48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 847, 895), 1E-16);
    }

    /*
     * Checks that propagate methods in quantitative trait
     * value for E's and D's inside QuaSSE class are working.
     *
     * Differs from 'testPropagateChOneCh4096QuaSSETest' inside
     * 'PropagatesQuaSSETest' because it relies on the QuaSSE class
     * correctly initializing all its dimensions and E's and D's.
     *
     * Test is done on a bifurcating tree (but just looking at a single
     * branch here, "sp1"), over a single dt = 0.01, with dx = 0.0005,
     * and nXbins = 4096 (high res). The trait value of "sp1" is set to 0.0.
     *
     * The method for FFTing used here is JavaFftService, so arrays
     * have values spread over while ignoring even indices
     */
    @Test
    public void testIntegrateOneBranchHiRes4096BinsOutsideClassJustX() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsHiAtNode;

         // we'll test the integration outside the class
        esDsHiAtNode = q1024.getEsDsAtNode(nodeIdx, false);
        double[][] fftBufferEsDsAtNode = new double[2][esDsHiAtNode[0].length];

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are used in propagate in x
         */
        double[][] esDsHiAtNodeInitial = new double[esDsHiAtNode.length][esDsHiAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsHiAtNode[0], esDsHiAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsHiAtNode[1], esDsHiAtNode[1].length); // D

        /*
         * fY is computed and FFTed in initialization, by QuaSSEProcess;
         * here we are just grabbing it to verify its values in the asserts
         * below
         */
        double aDt = 0.01;
        q1024.populatefY(aDt, true, true, false, false);
        double[] fftedfY = q1024.getfftFY(false);

        // copying fY for assert (leaving original one inside class untouched)
        double[] fftedfY4Assert = new double[fftedfY.length]; // just for test, not used in propagate in X
        for (int i=0; i < fftedfY.length; i++) fftedfY4Assert[i] = fftedfY[i];

        // everyOtherToHeadInPlace(fftedfY4Assert, q1024.getnXbins(false),0, 0, 2, 1.0); // getting real part for assert below

        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        double[][] scratchAtNode = new double[2][esDsHiAtNode[0].length];
        q1024.propagateXInPlace(esDsHiAtNode, fftBufferEsDsAtNode, scratchAtNode, false);

        // if using JTransforms
        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        // double[][] scratchAtNode = q1024.getScratchAtNode(nodeIdx, false); // sp1
        // q1024.propagateXInPlaceJTransforms(esDsHiAtNode, scratchAtNode, false);

        esDsHiAtNode = q1024.getEsDsAtNode(nodeIdx, false);
        double[] esHiAtNode = esDsHiAtNode[0];
        double[] dsHiAtNode = esDsHiAtNode[1];

        // System.out.println("Final D's: " + Arrays.toString(dsHiAtNode));

        double[] expectedInitialDsFirst48 = new double[] { 0.000353647683540531,0,0.000357627448646487,0,0.000361649739618714,0,0.000365714984190948,0,0.000369823614096463,0,0.000373976065102163,0,0.000378172777042955,0,0.000382414193856355,0,0.000386700763617376,0,0.000391032938573662,0,0.000395411175180895,0,0.000399835934138456,0,0.000404307680425366,0,0.000408826883336479,0,0.000413394016518956,0,0.00041800955800901,0,0.000422673990268911,0,0.00042738780022428,0,0.000432151479301657,0,0.000436965523466329,0,0.000441830433260469,0,0.000446746713841521,0,0.000451714875020898,0,0.000456735431302938,0,0.000461808901924177,0,0.00046693581089287,0,0.000472116687028841,0,0.000477352064003592,0,0.000482642480380733,0,0.000487988479656683,0,0.000493390610301674,0,0.000498849425801072,0,0.000504365484696964,0,0.000509939350630071,0,0.000515571592381964,0,0.000521262783917567,0,0.00052701350442799,0,0.000532824338373647,0,0.000538695875527712,0,0.000544628711019852,0,0.000550623445380317,0,0.000556680684584297,0,0.000562801040096648,0,0.000568985128916886,0,0.000575233573624552,0,0.000581547002424852,0,0.000587926049194663,0,0.000594371353528846,0 };
        double[] expectedFFTedfY = new double[] { 1, 0, 0.999247292368191, 0, 0.996992567179915, 0, 0.993245992007481, 0, 0.988024427909734, 0, 0.981351303026827, 0, 0.973256437474405, 0, 0.963775821368549, 0, 0.952951348293497, 0, 0.940830506973933, 0, 0.92746603432656, 0, 0.912915533436717, 0, 0.897241060330147, 0, 0.880508683684162, 0, 0.862788021843104, 0, 0.844151761668242, 0, 0.824675163860569, 0, 0.804435559446082, 0, 0.783511842107391, 0, 0.761983960984176, 0, 0.739932418450195, 0, 0.717437777209001, 0, 0.694580180837756, 0, 0.671438891652661, 0, 0.64809184947515, 0, 0.624615254550175, 0, 0.601083177512084, 0, 0.577567198915383, 0, 0.554136080452835, 0, 0.530855469577789, 0, 0.507787638837061, 0, 0.484991260810779, 0, 0.462521219151724, 0, 0.440428455824044, 0, 0.418759854264325, 0, 0.39755815783139, 0, 0.376861922578424, 0, 0.356705503075537, 0, 0.3371190697353, 0, 0.318128655850257, 0, 0.299756232341542, 0, 0.282019808042489, 0, 0.264933553200873, 0, 0.248507943778199, 0, 0.232749924053504, 0, 0.217663085001486, 0, 0.203247855908884, 0, 0.189501706717032, 0, 0.17641935863023, 0, 0.163993000606569, 0, 0.152212509447405, 0, 0.14106567132108, 0, 0.130538402692698, 0, 0.120614968781996, 0, 0.111278197832449, 0, 0.102509689644048, 0, 0.0942900169966733, 0, 0.0865989187681513, 0, 0.0794154837283513, 0, 0.0727183241657203, 0, 0.0664857386734408, 0, 0.0606958635869983, 0, 0.0553268127217662, 0, 0.0503568052068668, 0, 0.045764281348922, 0, 0.0415280065854677, 0, 0.0376271637021404, 0, 0.0340414335898102, 0, 0.0307510649074428, 0, 0.0277369330935997, 0, 0.0249805892342995, 0, 0.02246429934781, 0, 0.0201710746882624, 0, 0.0180846937003869, 0, 0.0161897162778361, 0, 0.0144714909882638, 0, 0.0129161559303678, 0, 0.0115106338823706, 0, 0.0102426223887588, 0, 0.00910057941344994, 0, 0.00807370516375022, 0, 0.0071519206614005, 0, 0.0063258436054819, 0, 0.00558676203776515, 0, 0.00492660628496708, 0, 0.00433791961500649, 0, 0.00381382800634398, 0, 0.00334800939141312, 0, 0.00293466269748287, 0, 0.0025684769714751, 0, 0.00224460083965355, 0, 0.00195861251901208, 0, 0.00170649056486081, 0, 0.001484585508736, 0, 0.00128959251247564, 0, 0.00111852513821068, 0 };
        double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsFirst48AfterPropX = new double[] { 0.000353647683540531, 0, 0.000357627448646487, 0, 0.000361649739618714, 0, 0.000365714984190948, 0, 0.000369823614096463, 0, 0.000373976065102163, 0, 0.000378172777042955, 0, 0.000382414193856355, 0, 0.000386700763617376, 0, 0.000391032938573662, 0, 0.000395411175180895, 0, 0.000399835934138456, 0, 0.000404307680425366, 0, 0.000408826883336479, 0, 0.000413394016518956, 0, 0.00041800955800901, 0, 0.000422673990268911, 0, 0.00042738780022428, 0, 0.000432151479301657, 0, 0.000436965523466329, 0, 0.000441830433260469, 0, 0.000446746713841521, 0, 0.000451714875020898, 0, 0.000456735431302938, 0, 0.000461808901924177, 0, 0.00046693581089287, 0, 0.000472116687028841, 0, 0.000477352064003592, 0, 0.000482642480380733, 0, 0.000487988479656683, 0, 0.000493390610301674, 0, 0.000498849425801072, 0, 0.000504365484696964, 0, 0.000509939350630071, 0, 0.000515571592381964, 0, 0.000521262783917567, 0, 0.00052701350442799, 0, 0.000532824338373647, 0, 0.000538695875527712, 0, 0.000544628711019852, 0, 0.000550623445380317, 0, 0.000556680684584297, 0, 0.000562801040096648, 0, 0.000568985128916886, 0, 0.000575233573624552, 0, 0.000581547002424852, 0, 0.000587926049194663, 0, 0.000594371353528846, 0 };
        double[] expectedSp1DsLater48AfterPropX = new double[] { 1.13579683037393, 0, 1.14139979445304, 0, 1.14702325798403, 0, 1.15266725177838, 0, 1.15833180641492, 0, 1.16401695223792, 0, 1.16972271935522, 0, 1.17544913763623, 0, 1.18119623671012, 0, 1.18696404596383, 0, 1.19275259454019, 0, 1.19856191133597, 0, 1.20439202499999, 0, 1.21024296393119, 0, 1.21611475627669, 0, 1.22200742992987, 0, 1.22792101252845, 0, 1.23385553145257, 0, 1.23981101382282, 0, 1.24578748649837, 0, 1.25178497607497, 0, 1.25780350888308, 0, 1.26384311098589, 0, 1.26990380817741, 0, 1.27598562598051, 0, 1.28208858964501, 0, 1.28821272414573, 0, 1.29435805418053, 0, 1.30052460416843, 0, 1.30671239824758, 0, 1.31292146027342, 0, 1.31915181381665, 0, 1.32540348216137, 0, 1.33167648830306, 0, 1.33797085494669, 0, 1.34428660450477, 0, 1.35062375909539, 0, 1.35698234054031, 0, 1.36336237036299, 0, 1.36976386978664, 0, 1.37618685973233, 0, 1.38263136081698, 0, 1.38909739335148, 0, 1.39558497733871, 0, 1.40209413247163, 0, 1.40862487813131, 0, 1.41517723338502, 0, 1.42175121698428, 0 };

        Assert.assertArrayEquals(expectedInitialDsFirst48, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 96), 1E-16);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 192), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esHiAtNode, 0, 96), 1E-16);
        Assert.assertArrayEquals(expectedSp1DsFirst48AfterPropX, Arrays.copyOfRange(dsHiAtNode, 0, 96), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsLater48AfterPropX, Arrays.copyOfRange(dsHiAtNode, 2000, 2096), 1E-12); // here things start to get different between R and Java and different CPU architectures (the fft-ed fY's multiplying the D's are very small!)

        // if using JTransforms
        // double[] expectedInitialDsFirst48 = new double[] { 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846 };
        // double[] expectedFFTedfY = new double[] { 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368549, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577789, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032 };
        // double[] expectedSp1EsAfterPropX = new double[] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1DsFirst48AfterPropX = new double[] { 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846 };
        // double[] expectedSp1DsLater48AfterPropX = new double[] { 1.13579683037393, 1.14139979445304, 1.14702325798403, 1.15266725177838, 1.15833180641492, 1.16401695223792, 1.16972271935522, 1.17544913763623, 1.18119623671012, 1.18696404596383, 1.19275259454019, 1.19856191133597, 1.20439202499999, 1.21024296393119, 1.21611475627669, 1.22200742992987, 1.22792101252845, 1.23385553145257, 1.23981101382282, 1.24578748649837, 1.25178497607497, 1.25780350888308, 1.26384311098589, 1.26990380817741, 1.27598562598051, 1.28208858964501, 1.28821272414573, 1.29435805418053, 1.30052460416843, 1.30671239824758, 1.31292146027342, 1.31915181381665, 1.32540348216137, 1.33167648830306, 1.33797085494669, 1.34428660450477, 1.35062375909539, 1.35698234054031, 1.36336237036299, 1.36976386978664, 1.37618685973233, 1.38263136081698, 1.38909739335148, 1.39558497733871, 1.40209413247163, 1.40862487813131, 1.41517723338502, 1.42175121698428 };
        // Assert.assertArrayEquals(expectedInitialDsFirst48, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 48), 1E-16);
        // Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 48), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropX, Arrays.copyOfRange(esLoAtNode, 0, 48), 1E-16);
        // Assert.assertArrayEquals(expectedSp1DsFirst48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 0, 48), 1E-16);
        // Assert.assertArrayEquals(expectedSp1DsLater48AfterPropX, Arrays.copyOfRange(dsLoAtNode, 1000, 1048), 1E-12); // here things start to get different between R and Java and different CPU architectures (the fft-ed fY's multiplying the D's are very small!)
}

    /*
     * Checks that propagate methods in both time and quantitative trait value (X)
     * for E's and D's inside QuaSSE class are working.
     *
     * Test is done on a bifurcating tree (but just looking at a single
     * branch here, "sp1"), over a single dt = 0.01, with dx = 0.0005,
     * and nXbins = 48 (low res). The trait value of "sp1" is set to 0.0.
     *
     * The method for FFTing used here is JavaFftService, so arrays
     * have values spread over while ignoring even indices
     */
    @Test
    public void testBothXandTPropagateMethodsInsideClassOneBranchLoRes48Bins() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);
        double[][] fftBufferEsDsAtNode = new double[2][esDsLoAtNode[0].length];

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are the starting point of everything
         */
        double[][] esDsLoAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        esDsLoAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsLoAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D

        // propagating in t
        // preparing scratch for propagate in t
        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];

        // just propagate in t, in place
        boolean jtransforms = false;
        q32Dt001.propagateTInPlace(esDsLoAtNode, scratchAtNode, dt001, true, jtransforms);

        // grabbing intermediate to see if it's all good
        double[][] esDsLoAtNodeAfterPropT = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        // making deep copy
        for (int ithDim=0; ithDim < 2; ithDim++) {
            for (int i = 0; i < esDsLoAtNode[ithDim].length; i++) {
                esDsLoAtNodeAfterPropT[ithDim][i] = esDsLoAtNode[ithDim][i];
            }
        }

        // propagating in x
        /*
         * fY is computed and FFTed in initialization, by QuaSSEProcess;
         * here we are just grabbing it to verify its values in the asserts
         * below
         */
        double aDt = 0.01;
        boolean forceRecalcKernel = true;
        q32Dt001.populatefY(aDt, forceRecalcKernel, true, true, false);
        double[] fftedfY = q32Dt001.getfftFY(true);

        // copying fY for assert (leaving original one inside class untouched)
        double[] fftedfY4Assert = new double[fftedfY.length]; // just for test, not used in propagate in X
        for (int i=0; i < fftedfY.length; i++) fftedfY4Assert[i] = fftedfY[i];

        // everyOtherToHeadInPlace(fftedfY4Assert, q32Dt001.getnXbins(true),0, 0, 2, 1.0); // getting real part for assert below

        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        q32Dt001.propagateXInPlace(esDsLoAtNode, fftBufferEsDsAtNode, scratchAtNode, true);
        esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);

        // if using JTransforms
        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        // q32Dt001.propagateXInPlaceJTransforms(esDsLoAtNode, scratchAtNode, true);
        // esDsLoAtNode = q32Dt001.getEsDsAtNode(nodeIdx, true);

        double[] expectedInitialDs = new double[] { 0.709491856924629, 0, 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.9945493127149, 0, 3.88372109966426, 0, 4.83941449038287, 0, 5.79383105522966, 0, 6.66449205783599, 0, 7.36540280606647, 0, 7.82085387950912, 0, 7.97884560802865, 0, 7.82085387950912, 0, 7.36540280606647, 0, 6.66449205783599, 0, 5.79383105522966, 0, 4.83941449038287, 0, 3.88372109966426, 0, 2.9945493127149, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedFFTedfY = new double[] { 1, 0, 0.999744507157237, 0, 0.998987847110852, 0, 0.997759097982437, 0, 0.996105480062091, 0, 0.994090541136289, 0, 0.991791714355113, 0, 0.989297342492751, 0, 0.986703282961364, 0, 0.984109224049216, 0, 0.981614853950298, 0, 0.979316029808302, 0, 0.977301093995625, 0, 0.975647479188405, 0, 0.97441873269917, 0, 0.973662074416229, 0, 0.973406582192705, 0, 0.973662074416229, 0, 0.97441873269917, 0, 0.975647479188405, 0, 0.977301093995625, 0, 0.979316029808302, 0, 0.981614853950298, 0, 0.984109224049216, 0, 0.986703282961364, 0, 0.989297342492751, 0, 0.991791714355113, 0, 0.994090541136289, 0, 0.996105480062091, 0, 0.997759097982437, 0, 0.998987847110852, 0, 0.999744507157237, 0 };
        double[] expectedSp1EsAfterPropT = new double[] { 0.000299740440744498, 0, 0.000299739520470888, 0, 0.000299738597335782, 0, 0.00029973767161558, 0, 0.000299736743589792, 0, 0.000299735813540893, 0, 0.000299734881753716, 0, 0.000299733948515344, 0, 0.00029973301411462, 0, 0.00029973207884177, 0, 0.000299731142988294, 0, 0.000299730206846208, 0, 0.00029972927070804, 0, 0.000299728334866242, 0, 0.000299727399612858, 0, 0.000299726465239364, 0, 0.000299725532035927, 0, 0.000299724600291355, 0, 0.000299723670292635, 0, 0.000299722742324704, 0, 0.000299721816669763, 0, 0.000299720893607378, 0, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1EsAfterPropTandX = new double[] { 0.000299740440744498, 0, 0.000299739520470888, 0, 0.000299738597335782, 0, 0.00029973767161558, 0, 0.000299736743576341, 0, 0.000299735813529336, 0, 0.000299734881744068, 0, 0.000299733948507616, 0, 0.000299733014108822, 0, 0.000299732078837909, 0, 0.000299731142986375, 0, 0.000299730206846234, 0, 0.000299729270710011, 0, 0.000299728334870154, 0, 0.000299727399618708, 0, 0.000299726465247143, 0, 0.000299725532045626, 0, 0.000299724600302962, 0, 0.000299723670306136, 0, 0.000299722742324704, 0, 0.000299721816669763, 0, 0.000299720893607378, 0, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropT = new double[] { 0.708264611249434, 0, 1.07794488915543, 0, 1.57625248872262, 0, 2.21453845459856, 0, 2.98929572353432, 0, 3.87688349797518, 0, 4.83086427268124, 0, 5.78355856417974, 0, 6.65263439389469, 0, 7.35225216820117, 0, 7.80684129110362, 0, 7.96450018578724, 0, 7.80674374052118, 0, 7.35206845754541, 0, 6.65238511637526, 0, 5.78326972079126, 0, 4.83056283595408, 0, 3.87659337350287, 0, 2.98903491521347, 0, 2.21431781312691, 0, 1.57607596745573, 0, 1.07781089197137, 0, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropTandX = new double[] { 0.708264611249434, 0, 1.07794488915543, 0, 1.57625248872262, 0, 2.21453845459856, 0, 2.99004586159941, 0, 3.87732490267096, 0, 4.83085571964462, 0, 5.78300263831972, 0, 6.65150777531997, 0, 7.35062312893061, 0, 7.80486719135139, 0, 7.96240319031702, 0, 7.80476971649718, 0, 7.35043955518597, 0, 6.65125867066456, 0, 5.78271397425439, 0, 4.83055444185051, 0, 3.87703489789509, 0, 2.98978512541793, 0, 2.21431781312691, 0, 1.57607596745573, 0, 1.07781089197137, 0, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsLoAtNodeInitial[1], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[0], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[0], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[1], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[1], 0, 64), 1E-14);

        // // if using JTransforms
        // double[] expectedInitialDs = new double[] { 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedFFTedfY = new double[] { 1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216, 0.981614853950298, 0.979316029808302, 0.977301093995625, 0.975647479188405, 0.97441873269917, 0.973662074416229, 0.973406582192705, 0.973662074416229, 0.97441873269917, 0.975647479188405, 0.977301093995625, 0.979316029808302, 0.981614853950298, 0.984109224049216, 0.986703282961364, 0.989297342492751, 0.991791714355113, 0.994090541136289, 0.996105480062091, 0.997759097982437, 0.998987847110852, 0.999744507157237 };
        // // note how first nkl (=4) and last nkr (=4) are the same
        // double[] expectedSp1EsAfterPropT = new double[] { 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1EsAfterPropTandX = new double[] { 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743576341, 0.000299735813529336, 0.000299734881744068, 0.000299733948507616, 0.000299733014108822, 0.000299732078837909, 0.000299731142986375, 0.000299730206846234, 0.000299729270710011, 0.000299728334870154, 0.000299727399618708, 0.000299726465247143, 0.000299725532045626, 0.000299724600302962, 0.000299723670306136, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // note how first nkl (=4) and last nkr (=4) are the same
        // double[] expectedSp1DsAfterPropT = new double[] { 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052118, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1DsAfterPropTandX = new double[] { 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.99004586159941, 3.87732490267096, 4.83085571964462, 5.78300263831972, 6.65150777531997, 7.35062312893061, 7.80486719135139, 7.96240319031702, 7.80476971649718, 7.35043955518597, 6.65125867066456, 5.78271397425439, 4.83055444185051, 3.87703489789509, 2.98978512541793, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsLoAtNodeInitial[1], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[0], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[1], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[0], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1DsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[1], 0, 32), 1E-14);
    }

    /*
     * Checks that propagate methods in both time and quantitative trait value (X)
     * for E's and D's inside QuaSSE class are working.
     *
     * Test is done on a bifurcating tree (but just looking at a single
     * branch here, "sp1"), over one larger dt = 0.02, with dx = 0.0005,
     * and nXbins = 48 (low res). The trait value of "sp1" is set to 0.0.
     *
     * The method for FFTing used here is JavaFftService, so arrays
     * have values spread over while ignoring even indices.
     */
    @Test
    public void testBothXandTPropagateMethodsInsideClassOneBranchLoRes32BinsDt002() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        double[][] esDsLoAtNode;

        /*
         * we'll test the integration outside the class
         */
        esDsLoAtNode = q32Dt002.getEsDsAtNode(nodeIdx, true);
        double[][] fftBufferEsDsAtNode = new double[2][esDsLoAtNode[0].length];

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are the starting point of everything
         */
        double[][] esDsLoAtNodeInitial = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        esDsLoAtNodeInitial[0] = Arrays.copyOf(esDsLoAtNode[0], esDsLoAtNode[0].length); // E
        esDsLoAtNodeInitial[1] = Arrays.copyOf(esDsLoAtNode[1], esDsLoAtNode[1].length); // D

        // propagating in t
        // preparing scratch for propagate in t
        double[][] scratchAtNode = new double[2][esDsLoAtNode[0].length];

        // just propagate in t, in place
        boolean jtransforms = false;
        q32Dt002.propagateTInPlace(esDsLoAtNode, scratchAtNode, dt002,true, jtransforms);

        // grabbing intermediate to see if it's all good
        double[][] esDsLoAtNodeAfterPropT = new double[esDsLoAtNode.length][esDsLoAtNode[0].length];
        // making deep copy
        for (int ithDim=0; ithDim < 2; ithDim++) {
            for (int i = 0; i < esDsLoAtNode[ithDim].length; i++) {
                esDsLoAtNodeAfterPropT[ithDim][i] = esDsLoAtNode[ithDim][i];
            }
        }

        // propagating in x
        /*
         * fY is computed and FFTed in initialization, by QuaSSEProcess;
         * here we are just grabbing it to verify its values in the asserts
         * below
         */
        double aDt = 0.02;
        boolean forceRecalcKernel = true;
        q32Dt002.populatefY(aDt, forceRecalcKernel, true, true, false);
        double[] fftedfY = q32Dt002.getfftFY(true);

        // copying fY for assert (leaving original one inside class untouched)
        double[] fftedfY4Assert = new double[fftedfY.length]; // just for test, not used in propagate in X
        for (int i=0; i < fftedfY.length; i++) fftedfY4Assert[i] = fftedfY[i];

        // everyOtherToHeadInPlace(fftedfY4Assert, q32Dt002.getnXbins(true),0, 0, 2, 1.0); // getting real part for assert below

        q32Dt002.propagateXInPlace(esDsLoAtNode, fftBufferEsDsAtNode, scratchAtNode, true);
        esDsLoAtNode = q32Dt002.getEsDsAtNode(nodeIdx, true);

        // if using JTransforms
        // just propagate in x, in place
        // calling the actual method we want to test after making sure the FFTed fY and the initial D's are correct
        // q32Dt002.propagateXInPlaceJTransforms(esDsLoAtNode, scratchAtNode, true);
        // esDsLoAtNode = q32Dt002.getEsDsAtNode(nodeIdx, true);

        double[] expectedInitialDs = new double[] { 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.9945493127149, 0, 3.88372109966426, 0, 4.83941449038287, 0, 5.79383105522965, 0, 6.66449205783599, 0, 7.36540280606647, 0, 7.82085387950912, 0, 7.97884560802865, 0, 7.82085387950912, 0, 7.36540280606647, 0, 6.66449205783599, 0, 5.79383105522965, 0, 4.83941449038287, 0, 3.88372109966426, 0, 2.9945493127149, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedFFTedfY = new double[] { 1, 0, 0.997284635663215, 0, 0.989243568247244, 0, 0.976187735593861, 0, 0.958621745686731, 0, 0.937224046438699, 0, 0.91282033618305, 0, 0.886351315420652, 0, 0.858836097834951, 0, 0.831332753409919, 0, 0.80489754454845, 0, 0.780544437378578, 0, 0.759206428540079, 0, 0.741700129042481, 0, 0.728694899474875, 0, 0.720687643959832, 0, 0.717984152783717, 0, 0.720687643959832, 0, 0.728694899474875, 0, 0.741700129042481, 0, 0.759206428540079, 0, 0.780544437378578, 0, 0.80489754454845, 0, 0.831332753409919, 0, 0.858836097834951, 0, 0.886351315420652, 0, 0.91282033618305, 0, 0.937224046438699, 0, 0.958621745686731, 0, 0.976187735593861, 0, 0.989243568247244, 0, 0.997284635663215, 0 };

        // note how first nkl (=5) and last nkr (=5) are the same
        double[] expectedSp1EsAfterPropT = new double[] { 0.000598958856744621, 0, 0.000598955169219663, 0, 0.000598951471383546, 0, 0.000598947764353066, 0, 0.000598944049256322, 0, 0.000598940327231526, 0, 0.000598936599425321, 0, 0.000598932866991505, 0, 0.000598929131089762, 0, 0.000598925392884107, 0, 0.000598921653541243, 0, 0.000598917914229581, 0, 0.000598914176117241, 0, 0.000598910440370952, 0, 0.000598906708154486, 0, 0.000598902980627182, 0, 0.000598899258942585, 0, 0.000598895544247008, 0, 0.000598891837677968, 0, 0.000598888140363113, 0, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1EsAfterPropTandX = new double[] { 0.000598958856744621, 0, 0.000598955169219663, 0, 0.000598951471383546, 0, 0.000598947764353066, 0, 0.000598944049256322, 0, 0.000598940326823012, 0, 0.000598936599098334, 0, 0.000598932866746461, 0, 0.000598929130926967, 0, 0.000598925392803752, 0, 0.000598921653543447, 0, 0.000598917914314325, 0, 0.000598914176284426, 0, 0.000598910440620369, 0, 0.000598906708485821, 0, 0.000598902981040026, 0, 0.000598899258942585, 0, 0.000598895544247008, 0, 0.000598891837677968, 0, 0.000598888140363113, 0, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        // note how first nkl (=5) and last nkr (=5) are the same
        double[] expectedSp1DsAfterPropT = new double[] { 1.07607462857111, 0, 1.57350796410007, 0, 2.2106689156898, 0, 2.98405395410774, 0, 3.87006132459304, 0, 4.82233340368848, 0, 5.77330938669843, 0, 6.64080371901192, 0, 7.33913154877464, 0, 7.79286078082761, 0, 7.95018769794391, 0, 7.79266608851656, 0, 7.33876489743091, 0, 6.64030620874454, 0, 5.77273291052126, 0, 4.82173179361503, 0, 3.86948229161728, 0, 2.98353343056498, 0, 2.21022855749815, 0, 1.5731556614014, 0, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedSp1DsAfterPropTandX = new double[] { 1.07607462857111, 0, 1.57350796410007, 0, 2.2106689156898, 0, 2.98405395410774, 0, 3.87006132459304, 0, 4.82224125131593, 0, 5.76741044168543, 0, 6.62885082364537, 0, 7.32184914756387, 0, 7.77191831137883, 0, 7.92794195899393, 0, 7.77172522522963, 0, 7.3214854000561, 0, 6.62835697979005, 0, 5.76683776896087, 0, 4.8216430122779, 0, 3.86948229161728, 0, 2.98353343056498, 0, 2.21022855749815, 0, 1.5731556614014, 0, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsLoAtNodeInitial[1], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[0], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[0], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[1], 0, 64), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[1], 0, 64), 1E-14);

        // // if using JTransforms
        // double[] expectedInitialDs = new double[] { 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedFFTedfY = new double[] { 1, 0.997284635663215, 0.989243568247244, 0.976187735593861, 0.958621745686731, 0.937224046438699, 0.91282033618305, 0.886351315420652, 0.858836097834951, 0.831332753409919, 0.80489754454845, 0.780544437378578, 0.759206428540079, 0.741700129042481, 0.728694899474875, 0.720687643959832, 0.717984152783717, 0.720687643959832, 0.728694899474875, 0.741700129042481, 0.759206428540079, 0.780544437378578, 0.80489754454845, 0.831332753409919, 0.858836097834951, 0.886351315420652, 0.91282033618305, 0.937224046438699, 0.958621745686731, 0.976187735593861, 0.989243568247244, 0.997284635663215 };
        // // note how first nkl (=5) and last nkr (=5) are the same
        // double[] expectedSp1EsAfterPropT = new double[] { 0.000598958856744621, 0.000598955169219663, 0.000598951471383546, 0.000598947764353066, 0.000598944049256322, 0.000598940327231526, 0.000598936599425321, 0.000598932866991505, 0.000598929131089762, 0.000598925392884107, 0.000598921653541243, 0.000598917914229581, 0.000598914176117241, 0.000598910440370952, 0.000598906708154486, 0.000598902980627182, 0.000598899258942585, 0.000598895544247008, 0.000598891837677968, 0.000598888140363113, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1EsAfterPropTandX = new double[] { 0.000598958856744621, 0.000598955169219663, 0.000598951471383546, 0.000598947764353066, 0.000598944049256322, 0.000598940326823012, 0.000598936599098334, 0.000598932866746461, 0.000598929130926967, 0.000598925392803752, 0.000598921653543447, 0.000598917914314325, 0.000598914176284426, 0.000598910440620369, 0.000598906708485821, 0.000598902981040026, 0.000598899258942585, 0.000598895544247008, 0.000598891837677968, 0.000598888140363113, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // // note how first nkl (=5) and last nkr (=5) are the same
        // double[] expectedSp1DsAfterPropT = new double[] { 1.07607462857111, 1.57350796410007, 2.2106689156898, 2.98405395410774, 3.87006132459304, 4.82233340368848, 5.77330938669843, 6.64080371901192, 7.33913154877464, 7.79286078082761, 7.95018769794391, 7.79266608851656, 7.33876489743091, 6.64030620874454, 5.77273291052126, 4.82173179361503, 3.86948229161728, 2.98353343056498, 2.21022855749815, 1.5731556614014, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] expectedSp1DsAfterPropTandX = new double[] { 1.07607462857111, 1.57350796410007, 2.2106689156898, 2.98405395410774, 3.87006132459304, 4.82224125131593, 5.76741044168543, 6.62885082364537, 7.32184914756387, 7.77191831137883, 7.92794195899393, 7.77172522522963, 7.3214854000561, 6.62835697979005, 5.76683776896087, 4.8216430122779, 3.86948229161728, 2.98353343056498, 2.21022855749815, 1.5731556614014, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // Assert.assertArrayEquals(expectedInitialDs, Arrays.copyOfRange(esDsLoAtNodeInitial[1], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedFFTedfY, Arrays.copyOfRange(fftedfY4Assert, 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[0], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[0], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1DsAfterPropT, Arrays.copyOfRange(esDsLoAtNodeAfterPropT[1], 0, 32), 1E-14);
        // Assert.assertArrayEquals(expectedSp1DsAfterPropTandX, Arrays.copyOfRange(esDsLoAtNode[1], 0, 32), 1E-14);
    }

    /*
     * Checks that integrateBranch method inside QuaSSE class is working
     * (it propagates in both time and quantitative trait value (X).
     *
     * Test is done on a bifurcating tree (but just looking at a single
     * branch here, "sp1"), over two dt's of 0.01 = 2 * 0.01 = 0.02,
     * with dx = 0.0005, and nXbins = 192 (high res).
     * The trait value of "sp1" is set to 0.0.
     *
     * The method for FFTing used here is JavaFftService, so arrays
     * have values spread over while ignoring even indices
     */
    @Test
    public void testIntegrateOneBranchHiRes32BinsInsideClassBothXandTDt002() {
        // we're going to look at sp1
        int nodeIdx = 0; // sp1
        Node sp1Node = bifTreeHeight002.getNode(nodeIdx);

        /*
         * let us grab the E's and D's before integration
         */
        double[][] esDsHiAtNode;
        esDsHiAtNode = q32Dt002.getEsDsAtNode(nodeIdx, false);

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are the starting point of everything
         */
        double[][] esDsHiAtNodeInitial = new double[esDsHiAtNode.length][esDsHiAtNode[0].length];
        esDsHiAtNodeInitial[0] = Arrays.copyOf(esDsHiAtNode[0], esDsHiAtNode[0].length); // E
        esDsHiAtNodeInitial[1] = Arrays.copyOf(esDsHiAtNode[1], esDsHiAtNode[1].length); // D

        /*
         * now we integrate over 2 dt's = 2 * 0.01 = 0.02, inside class!
         * note that inside processBranch, we normalize Es and Ds
         */
        boolean jtransforms = false;
        boolean forceRecalcKernel = false;
        q32Dt002.processBranch(sp1Node, forceRecalcKernel, jtransforms);
        esDsHiAtNode = q32Dt002.getEsDsAtNode(nodeIdx, false);

        double[] expectedInitialDsFirst10 = new double[] { 0.791000831787404, 0, 0.879671919608544, 0, 0.975840371583655, 0, 1.07981933026376, 0, 1.19189412137632, 0, 1.31231629549353, 0, 1.44129748672436, 0, 1.57900316601788, 0, 1.72554637653023, 0, 1.88098154753774, 0 };
        double[] expectedInitialDsLater10 = new double[] { 1.88098154753774, 0, 1.72554637653023, 0, 1.57900316601788, 0, 1.44129748672436, 0, 1.31231629549353, 0, 1.19189412137632, 0, 1.07981933026376, 0, 0.975840371583657, 0, 0.879671919608545, 0, 0.791000831787405 };
        double[] expectedSp1EsAfterPropTandXFirst10 = new double[] { 0.000598961614956886, 0, 0.000598960696294835, 0, 0.000598959776885031, 0, 0.000598958856744621, 0, 0.00059895793589067, 0, 0.000598957014340337, 0, 0.000598956092110961, 0, 0.000598955169219663, 0, 0.000598954245683891, 0, 0.000598953321520975, 0 };
        double[] expectedSp1EsAfterPropTandXLater10 = new double[] { 0.000598889987793744, 0, 0.000598889063763026, 0, 0.000598888140363113, 0, 0.000598887217611451, 0, 0.000598886295525308, 0, 0.000598885374121974, 0, 0.000598884453418597, 0, 0.0005988835334325, 0, 0.000598882614180599, 0, 0.00059888169568006 };
        double[] expectedSp1DsAfterPropTandXFirst10 = new double[] { 0.816826443231728, 0, 0.908389791093209, 0, 1.00769467468308, 0, 1.11506438544145, 0, 1.23079348607031, 0, 1.35514165820029, 0, 1.48832736122564, 0, 1.63052138243079, 0, 1.78184036871647, 0, 1.9423404395558, 0 };
        double[] expectedSp1DsAfterPropTandXLater10 = new double[] { 1.94192952805763, 0, 1.78145241269027, 0, 1.63015631463554, 0, 1.4879849556747, 0, 1.35482154601338, 0, 1.23049517164847, 0, 1.11478726270473, 0, 1.00743804313279, 0, 0.908152871485605, 0, 0.816608392676392 };

        Assert.assertArrayEquals(expectedInitialDsFirst10, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 20), 1E-14);
        Assert.assertArrayEquals(expectedInitialDsLater10, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 154, 173), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropTandXFirst10, Arrays.copyOfRange(esDsHiAtNode[0], 0, 20), 1E-14);
        Assert.assertArrayEquals(expectedSp1EsAfterPropTandXLater10, Arrays.copyOfRange(esDsHiAtNode[0], 154, 173), 1E-14); // OK til here
        Assert.assertArrayEquals(expectedSp1DsAfterPropTandXFirst10, Arrays.copyOfRange(esDsHiAtNode[1], 0, 20), 1E-14);
        Assert.assertArrayEquals(expectedSp1DsAfterPropTandXLater10, Arrays.copyOfRange(esDsHiAtNode[1], 154, 173), 1E-14);

        // if using JTransforms
        // double[] expectedInitialDsFirst10 = new double[] { 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774 };
        // double[] expectedInitialDsLater10 = new double[] { 1.88098154753774, 1.72554637653023, 1.57900316601788, 1.44129748672436, 1.31231629549353, 1.19189412137632, 1.07981933026376, 0.975840371583657, 0.879671919608545, 0.791000831787405 };
        // double[] expectedSp1EsAfterPropTandXFirst10 = new double[] { 0.000598961614956886, 0.000598960696294835, 0.000598959776885031, 0.000598958856744621, 0.00059895793589067, 0.000598957014340337, 0.000598956092110961, 0.000598955169219663, 0.000598954245683891, 0.000598953321520975 };
        // double[] expectedSp1EsAfterPropTandXLater10 = new double[] { 0.000598889987793744, 0.000598889063763026, 0.000598888140363113, 0.000598887217611451, 0.000598886295525308, 0.000598885374121974, 0.000598884453418597, 0.0005988835334325, 0.000598882614180599, 0.00059888169568006 };
        // double[] expectedSp1DsAfterPropTandXFirst10 = new double[] { 0.816826443231728, 0.908389791093209, 1.00769467468308, 1.11506438544145, 1.23079348607031, 1.35514165820029, 1.48832736122564, 1.63052138243079, 1.78184036871647, 1.9423404395558 };
        // double[] expectedSp1DsAfterPropTandXLater10 = new double[] { 1.94192952805763, 1.78145241269027, 1.63015631463554, 1.4879849556747, 1.35482154601338, 1.23049517164847, 1.11478726270473, 1.00743804313279, 0.908152871485605, 0.816608392676392 };
        // Assert.assertArrayEquals(expectedInitialDsFirst10, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 0, 10), 1E-14);
        // Assert.assertArrayEquals(expectedInitialDsLater10, Arrays.copyOfRange(esDsHiAtNodeInitial[1], 77, 87), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropTandXFirst10, Arrays.copyOfRange(esDsHiAtNode[0], 0, 10), 1E-14);
        // Assert.assertArrayEquals(expectedSp1EsAfterPropTandXLater10, Arrays.copyOfRange(esDsHiAtNode[0], 77, 87), 1E-14); // OK til here
        // Assert.assertArrayEquals(expectedSp1DsAfterPropTandXFirst10, Arrays.copyOfRange(esDsHiAtNode[1], 0, 10), 1E-14);
        // Assert.assertArrayEquals(expectedSp1DsAfterPropTandXLater10, Arrays.copyOfRange(esDsHiAtNode[1], 77, 87), 1E-14);
    }

    /*
     * Checks that method for populating the root prior probabilities (on D's)
     * is doing the right thing when the chosen method is "Observed".
     * This is the default in diversitree, and assigns as the prior probability
     * for each bin of D the corresponding weight that D has with respect to the
     * sum of all D's at the root.
     */
    @Test
    public void testPriorProbAtRootObserved() {
        // if using JTransforms
        // double[] rootDs = new double[] {
        //         4.66345926030482e-05, 0.000175320546519251, 0.000608432288010798, 0.00197994686277693, 0.00584255914263498, 0.0159175108393565, 0.0400377135639036, 0.0929790087002835, 0.199352433863763, 0.394619574541096, 0.72119984930653, 1.21689198263463, 1.89568999042485, 2.72647162610174, 3.62036642668937, 4.43834919455334, 5.02351058794779, 5.24939946304328, 5.06439897859485, 4.51088490297466, 3.70945996798424, 2.81627005522746, 1.97008358430641, 1.27377969405492, 0.760235257721598, 0, 0, 0, 0, 0, 0, 0,
        //         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        // };

        double[] rootDs = new double[] { 4.66345926030482e-05, 0, 0.000175320546519251, 0, 0.000608432288010798, 0, 0.00197994686277693, 0, 0.00584255914263498, 0, 0.0159175108393565, 0, 0.0400377135639036, 0, 0.0929790087002835, 0, 0.199352433863763, 0, 0.394619574541096, 0, 0.72119984930653, 0, 1.21689198263463, 0, 1.89568999042485, 0, 2.72647162610174, 0, 3.62036642668937, 0, 4.43834919455334, 0, 5.02351058794779, 0, 5.24939946304328, 0, 5.06439897859485, 0, 4.51088490297466, 0, 3.70945996798424, 0, 2.81627005522746, 0, 1.97008358430641, 0, 1.27377969405492, 0, 0.760235257721598, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // the 14 extra 0' s at the end make the Ds array be 64-elements long, as it should be the case

        int nUsefulXBinsLowRes = q32Dt0005.getNUsefulTraitBins(true);
        int nXBinsRightRes = q32Dt0005.getnXbins(true); // not used in this test, but required by method

        q32Dt0005.initializePriorProbAtRoot(nUsefulXBinsLowRes); // initialize here so I don't have to do recursion
        q32Dt0005.setGot2LowRes();

        boolean jtransforms = false;
        q32Dt0005.populatePriorProbAtRoot(rootDs, q32Dt0005.getdXbin(), nXBinsRightRes, nUsefulXBinsLowRes, "Observed", jtransforms); // calculate root prior

        double[] priorProbAtRoot = q32Dt0005.getPriorProbsAtRoot("Observed");

        // the prior array inside the likelihood has nUsefulXBinsAtRightRes elements (where 'RightRes' is either high or low resolution)
        double[] expectedObservedPriorProbAtRoot = { 0.000101936764975178, 0.000383226449472287, 0.0013299487715952, 0.00432788980772701, 0.0127710256471166, 0.0347934756336926, 0.0875169004358447, 0.203239244270493, 0.435756829076961, 0.862583772672889, 1.57644305300714, 2.65995745025329, 4.14371594632748, 5.95968961768647, 7.9136199323705, 9.70161705011622, 10.9806988668854, 11.474460683722, 11.0700752296872, 9.86017006942937, 8.10836608266039, 6.15597655521472, 4.30633004611629, 2.7843061138813, 1.66176905311156 };

        Assert.assertArrayEquals(expectedObservedPriorProbAtRoot, priorProbAtRoot, 1E-13);
    }

    /*
     * Checks that method for populating the root prior probabilities (on D's)
     * is doing the right thing when the chosen method is "Flat".
     */
    @Test
    public void testPriorProbAtRootFlat() {
        // if using JTransforms
        // double[] rootDs = new double[] {
        //         4.66345926030482e-05, 0.000175320546519251, 0.000608432288010798, 0.00197994686277693, 0.00584255914263498, 0.0159175108393565, 0.0400377135639036, 0.0929790087002835, 0.199352433863763, 0.394619574541096, 0.72119984930653, 1.21689198263463, 1.89568999042485, 2.72647162610174, 3.62036642668937, 4.43834919455334, 5.02351058794779, 5.24939946304328, 5.06439897859485, 4.51088490297466, 3.70945996798424, 2.81627005522746, 1.97008358430641, 1.27377969405492, 0.760235257721598, 0, 0, 0, 0, 0, 0, 0,
        //         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
        // };

        double[] rootDs = new double[] { 0.000101936764975178, 0, 0.000383226449472287, 0, 0.0013299487715952, 0, 0.00432788980772701, 0, 0.0127710256471166, 0, 0.0347934756336926, 0, 0.0875169004358447, 0, 0.203239244270493, 0, 0.435756829076961, 0, 0.862583772672889, 0, 1.57644305300714, 0, 2.65995745025329, 0, 4.14371594632748, 0, 5.95968961768647, 0, 7.9136199323705, 0, 9.70161705011622, 0, 10.9806988668854, 0, 11.474460683722, 0, 11.0700752296872, 0, 9.86017006942937, 0, 8.10836608266039, 0, 6.15597655521472, 0, 4.30633004611629, 0, 2.7843061138813, 0, 1.66176905311156, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // the 14 extra 0' s at the end make the Ds array be 64-elements long, as it should be the case

        int nUsefulXBinsLowRes = q32Dt0005.getNUsefulTraitBins(true);
        int nXBinsLowRes = q32Dt0005.getnXbins(true);

        q32Dt0005.initializePriorProbAtRoot(nUsefulXBinsLowRes); // initialize here so I don't have to do recursion

        boolean jtransforms = false;
        q32Dt0005.populatePriorProbAtRoot(rootDs, q32Dt0005.getdXbin(), nXBinsLowRes, nUsefulXBinsLowRes, "Flat", jtransforms); // calculate root prior

        double[] priorProbAtRoot = q32Dt0005.getPriorProbsAtRoot("Flat");

        // the prior array inside the likelihood has nUsefulXBinsAtRightRes elements (where 'RightRes' is either high or low resolution)
        double[] expectedFlatPriorProbAtRoot = { 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129 };

        Assert.assertArrayEquals(expectedFlatPriorProbAtRoot, priorProbAtRoot, 1E-14);
    }

    /*
     * Checks correctness of calculations at the root after pruning is concluded.
     * Values of D's and E's come from R example and are hardcoded.
     * Prior probabilities at the root are computed by the likelihood class.
     */
    @Test
    public void testRootCalcProcedure() {
        double[] lambda = new double[] { 0.142555748318834, 0.143168001652175, 0.14378234991142, 0.144398610945538, 0.145016600268752, 0.145636131276292, 0.146257015465625, 0.146879062662624, 0.147502081252106, 0.148125878412146, 0.148750260351579, 0.149375032550049, 0.15, 0.150624967449951, 0.151249739648421, 0.151874121587854, 0.152497918747894, 0.153120937337376, 0.153742984534375, 0.154363868723708, 0.154983399731248, 0.155601389054462, 0.15621765008858, 0.156831998347825, 0.157444251681166 };

        // if using JTransforms
        // double[] rootDs = new double[] { 4.66345926030482e-05, 0.000175320546519251, 0.000608432288010798, 0.00197994686277693, 0.00584255914263498, 0.0159175108393565, 0.0400377135639036, 0.0929790087002835, 0.199352433863763, 0.394619574541096, 0.72119984930653, 1.21689198263463, 1.89568999042485, 2.72647162610174, 3.62036642668937, 4.43834919455334, 5.02351058794779, 5.24939946304328, 5.06439897859485, 4.51088490297466, 3.70945996798424, 2.81627005522746, 1.97008358430641, 1.27377969405492, 0.760235257721598, 0, 0, 0, 0, 0, 0, 0,
        //         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        // double[] rootEs = new double[] { 0.000299741357883797, 0.000299740440744553, 0.000299739520470931, 0.00029973859731951, 0.000299737671601055, 0.000299736743577048, 0.000299735813529964, 0.000299734881744545, 0.000299733948508, 0.000299733014109106, 0.000299732078838168, 0.000299731142986417, 0.000299730206846241, 0.00029972927070989, 0.000299728334869967, 0.000299727399618395, 0.000299726465246688, 0.000299725532045103, 0.000299724600302324, 0.000299723670305445, 0.000299722742339267, 0.000299721816686102, 0.00029972089360741, 0.000299719973413787, 0.000299719056361797, 0, 0, 0, 0, 0, 0, 0,
        //         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        double[] rootDs = new double[] { 4.66345926030482e-05, 0, 0.000175320546519251, 0, 0.000608432288010798, 0, 0.00197994686277693, 0, 0.00584255914263498, 0, 0.0159175108393565, 0, 0.0400377135639036, 0, 0.0929790087002835, 0, 0.199352433863763, 0, 0.394619574541096, 0, 0.72119984930653, 0, 1.21689198263463, 0, 1.89568999042485, 0, 2.72647162610174, 0, 3.62036642668937, 0, 4.43834919455334, 0, 5.02351058794779, 0, 5.24939946304328, 0, 5.06439897859485, 0, 4.51088490297466, 0, 3.70945996798424, 0, 2.81627005522746, 0, 1.97008358430641, 0, 1.27377969405492, 0, 0.760235257721598, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] rootEs = new double[] { 0.000299741357883797, 0, 0.000299740440744553, 0, 0.000299739520470931, 0, 0.00029973859731951, 0, 0.000299737671601055, 0, 0.000299736743577048, 0, 0.000299735813529964, 0, 0.000299734881744545, 0, 0.000299733948508, 0, 0.000299733014109106, 0, 0.000299732078838168, 0, 0.000299731142986417, 0, 0.000299730206846241, 0, 0.00029972927070989, 0, 0.000299728334869967, 0, 0.000299727399618395, 0, 0.000299726465246688, 0, 0.000299725532045103, 0, 0.000299724600302324, 0, 0.000299723670305445, 0, 0.000299722742339267, 0, 0.000299721816686102, 0, 0.00029972089360741, 0, 0.000299719973413787, 0, 0.000299719056361797, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        double[][] esDs = new double[2][rootDs.length];
        esDs[0] = rootEs;
        esDs[1] = rootDs;

        int nUsefulXBinsLowRes = q32Dt0005.getNUsefulTraitBins(true);
        int nXBinsLowRes = q32Dt0005.getnXbins(true); // not used in this test, but required by method

        q32Dt0005.initializePriorProbAtRoot(nUsefulXBinsLowRes); // initialize here so I don't have to do recursion

        boolean jtransforms = false;
        q32Dt0005.populatePriorProbAtRoot(rootDs, q32Dt0005.getdXbin(), nXBinsLowRes, nUsefulXBinsLowRes, "Observed", jtransforms); // calculate root prior

        double sumLogNormalizationFactors = -0.384716957244803;
        double logLik = q32Dt0005.getLogPFromRelevantObjects(esDs, sumLogNormalizationFactors, lambda, q32Dt001.getdXbin(), jtransforms);

        Assert.assertEquals(2.820698371976182, logLik, 1e-6); // Corrected R reference; original precision retained.
    }

    // Analytical root integrals detect spacing errors shared with diversitree references.
    // Both layouts need this check while root conditioning is implemented as grid quadrature.
    @Test
    public void testRootConditioningQuadrature() {
        for (boolean jtransforms : new boolean[] { false, true }) {
            for (int n : new int[] { 4, 8, 16 }) {
                double dx = 1.0 / n;
                for (boolean varying : new boolean[] { false, true }) {
                    // Fresh D arrays are required because root conditioning modifies them in place.
                    double[][] esDs = new double[2][2*n];
                    double[] lambda = new double[n];
                    for (int i = 0; i < n; ++i) {
                        int index = jtransforms ? i : 2*i;
                        esDs[1][index] = 2.0;
                        esDs[0][index] = varying && i % 2 == 1 ? .5 : 0.0;
                        lambda[i] = varying ? (i % 2 == 0 ? .2 : .6) : .15;
                    }
                    q32Dt0005.initializePriorProbAtRoot(n);
                    q32Dt0005.populatePriorProbAtRoot(esDs[1], dx, n, n, "Observed", jtransforms);
                    // D=2 implies p=1 over a unit interval. For varying rates, the conditioning
                    // integral is (.2*1^2 + .6*.5^2)/2=.175; the restored numerator is 2*3=6.
                    double expected = Math.log(6.0 / (varying ? .175 : .15));
                    Assert.assertEquals(expected, q32Dt0005.getLogPFromRelevantObjects(
                            esDs, Math.log(3), lambda, dx, jtransforms), 1e-12);
                }
            }
        }
    }

    /*
     * Checks log-likelihood and other internal quantities
     * for bifurcating tree (height = 0.01), with 32 quantitative trait
     * bins, and dtMax = 0.005.
     *
     * Half of the tree height is in high resolution, half in low resolution.
     */
    @Test
    public void testPruneBifTree32Bins() {

        // indices for high to low res transfer
        int[] hiLoIdxs4Transfer = q32Dt0005.getHiLoIdxs4Transfer();

        // let us grab the E's and D's before integration
        double[][] esDsHiAtNode0, esDsHiAtNode1;
        esDsHiAtNode0 = q32Dt0005.getEsDsAtNode(0, false); // sp1
        esDsHiAtNode1 = q32Dt0005.getEsDsAtNode(1, false); // sp2

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are the starting point of everything
         */
        double[][] esDsHiAtNodeInitial0 = new double[esDsHiAtNode0.length][esDsHiAtNode0[0].length];
        esDsHiAtNodeInitial0[0] = Arrays.copyOf(esDsHiAtNode0[0], esDsHiAtNode0[0].length); // E
        esDsHiAtNodeInitial0[1] = Arrays.copyOf(esDsHiAtNode0[1], esDsHiAtNode0[1].length); // D

        double[][] esDsHiAtNodeInitial1 = new double[esDsHiAtNode1.length][esDsHiAtNode1[0].length];
        esDsHiAtNodeInitial1[0] = Arrays.copyOf(esDsHiAtNode1[0], esDsHiAtNode1[0].length); // E
        esDsHiAtNodeInitial1[1] = Arrays.copyOf(esDsHiAtNode1[1], esDsHiAtNode1[1].length); // D

        int[] expectedHiLoIdxs4Transfer = new int[] { 6, 14, 22, 30, 38, 46, 54, 62, 70, 78, 86, 94, 102, 110, 118, 126, 134, 142, 150, 158, 166, 174, 182, 190, 198 };
        double[] expectedDsHiAtNodeInitialSp1 = { 0.308986942687903, 0, 0.350566009871371, 0, 0.396747087835907, 0, 0.447890605896858, 0, 0.504364398303888, 0, 0.566540754832024, 0, 0.634793036713348, 0, 0.709491856924629, 0, 0.791000831787404, 0, 0.879671919608544, 0, 0.975840371583655, 0, 1.07981933026376, 0, 1.19189412137632, 0, 1.31231629549353, 0, 1.44129748672436, 0, 1.57900316601788, 0, 1.72554637653023, 0, 1.88098154753774, 0, 2.04529849127956, 0, 2.21841669358911, 0, 2.40018001393971, 0, 2.59035191331783, 0, 2.7886113289072, 0, 2.9945493127149, 0, 3.20766654683839, 0, 3.42737184095615, 0, 3.65298170778044, 0, 3.88372109966426, 0, 4.11872537439949, 0, 4.35704354065101, 0, 4.59764281368466, 0, 4.83941449038287, 0, 5.08118112938378, 0, 5.3217049979751, 0, 5.55969772261993, 0, 5.79383105522966, 0, 6.02274864309609, 0, 6.24507866733522, 0, 6.45944719335828, 0, 6.66449205783599, 0, 6.85887710038768, 0, 7.04130653528599, 0, 7.21053924923296, 0, 7.36540280606647, 0, 7.50480693833876, 0, 7.62775630921048, 0, 7.73336233605698, 0, 7.82085387950912, 0, 7.88958661815778, 0, 7.93905094954024, 0, 7.96887828189528, 0, 7.97884560802865, 0, 7.96887828189528, 0, 7.93905094954024, 0, 7.88958661815778, 0, 7.82085387950912, 0, 7.73336233605698, 0, 7.62775630921048, 0, 7.50480693833876, 0, 7.36540280606647, 0, 7.21053924923296, 0, 7.04130653528599, 0, 6.85887710038768, 0, 6.66449205783599, 0, 6.45944719335828, 0, 6.24507866733522, 0, 6.02274864309609, 0, 5.79383105522965, 0, 5.55969772261993, 0, 5.32170499797509, 0, 5.08118112938378, 0, 4.83941449038287, 0, 4.59764281368466, 0, 4.35704354065101, 0, 4.11872537439949, 0, 3.88372109966426, 0, 3.65298170778044, 0, 3.42737184095615, 0, 3.20766654683839, 0, 2.9945493127149, 0, 2.7886113289072, 0, 2.59035191331783, 0, 2.40018001393971, 0, 2.21841669358911, 0, 2.04529849127956, 0, 1.88098154753774, 0, 1.72554637653023, 0, 1.57900316601788, 0, 1.44129748672436, 0, 1.31231629549353, 0, 1.19189412137632, 0, 1.07981933026376, 0, 0.975840371583655, 0, 0.879671919608544, 0, 0.791000831787404, 0, 0.709491856924628, 0, 0.634793036713349, 0, 0.566540754832024, 0, 0.504364398303888, 0, 0.447890605896858, 0, 0.396747087835907, 0, 0.350566009871371, 0, 0.308986942687903, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedDsHiAtNodeInitialSp2 = { 0.000254946647636669, 0, 0.000319674822138109, 0, 0.000399835934138456, 0, 0.000498849425801072, 0, 0.000620828141157003, 0, 0.000770703934841743, 0, 0.000954372730824099, 0, 0.0011788613551308, 0, 0.00145251860604505, 0, 0.00178523314354266, 0, 0.00218868086879601, 0, 0.00267660451529771, 0, 0.00326512817532484, 0, 0.00397310942785545, 0, 0.00482253160451986, 0, 0.0058389385158292, 0, 0.00705191364734891, 0, 0.00849560541101504, 0, 0.0102092994868837, 0, 0.0122380386022755, 0, 0.0146332892566062, 0, 0.0174536539009152, 0, 0.0207656259132282, 0, 0.0246443833694604, 0, 0.0291746160933349, 0, 0.0344513787810736, 0, 0.0405809611459954, 0, 0.0476817640292969, 0, 0.0558851682975889, 0, 0.0653363811239983, 0, 0.0761952419644361, 0, 0.08863696823876, 0, 0.102852818461079, 0, 0.119050648395517, 0, 0.137455333812279, 0, 0.158309031659599, 0, 0.181871250031821, 0, 0.208418696288452, 0, 0.238244872152104, 0, 0.271659384673712, 0, 0.308986942687903, 0, 0.350566009871371, 0, 0.396747087835906, 0, 0.447890605896858, 0, 0.504364398303888, 0, 0.566540754832024, 0, 0.634793036713348, 0, 0.709491856924629, 0, 0.791000831787404, 0, 0.879671919608544, 0, 0.975840371583655, 0, 1.07981933026376, 0, 1.19189412137632, 0, 1.31231629549353, 0, 1.44129748672436, 0, 1.57900316601788, 0, 1.72554637653023, 0, 1.88098154753774, 0, 2.04529849127956, 0, 2.21841669358911, 0, 2.40018001393971, 0, 2.59035191331783, 0, 2.7886113289072, 0, 2.9945493127149, 0, 3.20766654683839, 0, 3.42737184095615, 0, 3.65298170778044, 0, 3.88372109966426, 0, 4.11872537439949, 0, 4.35704354065101, 0, 4.59764281368466, 0, 4.83941449038287, 0, 5.08118112938378, 0, 5.32170499797509, 0, 5.55969772261993, 0, 5.79383105522965, 0, 6.02274864309609, 0, 6.24507866733522, 0, 6.45944719335828, 0, 6.66449205783599, 0, 6.85887710038768, 0, 7.04130653528599, 0, 7.21053924923296, 0, 7.36540280606647, 0, 7.50480693833876, 0, 7.62775630921048, 0, 7.73336233605698, 0, 7.82085387950912, 0, 7.88958661815778, 0, 7.93905094954024, 0, 7.96887828189528, 0, 7.97884560802865, 0, 7.96887828189528, 0, 7.93905094954024, 0, 7.88958661815778, 0, 7.82085387950912, 0, 7.73336233605698, 0, 7.62775630921048, 0, 7.50480693833876, 0, 7.36540280606647, 0, 7.21053924923296, 0, 7.04130653528599, 0, 6.85887710038768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        double logLik = q32Dt0005.calculateLogP();

        // if using JTransforms
        // int[] expectedHiLoIdxs4Transfer = new int[] { 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99 };
        // double[] expectedDsHiAtNodeInitialSp1 = new double[] { 0.308986942687903, 0.350566009871371, 0.396747087835907, 0.447890605896858, 0.504364398303888, 0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911, 2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839, 3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101, 4.59764281368466, 4.83941449038287, 5.08118112938378, 5.3217049979751, 5.55969772261993, 5.79383105522966, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599, 6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876, 7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024, 7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778, 7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647, 7.21053924923296, 7.04130653528599, 6.85887710038768, 6.66449205783599, 6.45944719335828, 6.24507866733522, 6.02274864309609, 5.79383105522965, 5.55969772261993, 5.32170499797509, 5.08118112938378, 4.83941449038287, 4.59764281368466, 4.35704354065101, 4.11872537439949, 3.88372109966426, 3.65298170778044, 3.42737184095615, 3.20766654683839, 2.9945493127149, 2.7886113289072, 2.59035191331783, 2.40018001393971, 2.21841669358911, 2.04529849127956, 1.88098154753774, 1.72554637653023, 1.57900316601788, 1.44129748672436, 1.31231629549353, 1.19189412137632, 1.07981933026376, 0.975840371583655, 0.879671919608544, 0.791000831787404, 0.709491856924628, 0.634793036713349, 0.566540754832024, 0.504364398303888, 0.447890605896858, 0.396747087835907, 0.350566009871371, 0.308986942687903 };
        // double[] expectedDsHiAtNodeInitialSp2 = new double[] { 0.000254946647636669, 0.000319674822138109, 0.000399835934138456, 0.000498849425801072, 0.000620828141157003, 0.000770703934841743, 0.000954372730824099, 0.0011788613551308, 0.00145251860604505, 0.00178523314354266, 0.00218868086879601, 0.00267660451529771, 0.00326512817532484, 0.00397310942785545, 0.00482253160451986, 0.0058389385158292, 0.00705191364734891, 0.00849560541101504, 0.0102092994868837, 0.0122380386022755, 0.0146332892566062, 0.0174536539009152, 0.0207656259132282, 0.0246443833694604, 0.0291746160933349, 0.0344513787810736, 0.0405809611459954, 0.0476817640292969, 0.0558851682975889, 0.0653363811239983, 0.0761952419644361, 0.08863696823876, 0.102852818461079, 0.119050648395517, 0.137455333812279, 0.158309031659599, 0.181871250031821, 0.208418696288452, 0.238244872152104, 0.271659384673712, 0.308986942687903, 0.350566009871371, 0.396747087835906, 0.447890605896858, 0.504364398303888, 0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911, 2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839, 3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101, 4.59764281368466, 4.83941449038287, 5.08118112938378, 5.32170499797509, 5.55969772261993, 5.79383105522965, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599, 6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876, 7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024, 7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778, 7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647, 7.21053924923296, 7.04130653528599, 6.85887710038768 };
        // Assert.assertArrayEquals(expectedDsHiAtNodeInitialSp1, Arrays.copyOfRange(esDsHiAtNodeInitial0[1], 0, 103), 1E-13);
        // Assert.assertArrayEquals(expectedDsHiAtNodeInitialSp2, Arrays.copyOfRange(esDsHiAtNodeInitial1[1], 0, 103), 1E-13);

        Assert.assertArrayEquals(expectedHiLoIdxs4Transfer, hiLoIdxs4Transfer);
        Assert.assertArrayEquals(expectedDsHiAtNodeInitialSp1, Arrays.copyOfRange(esDsHiAtNodeInitial0[1], 0, 256), 1E-13);
        Assert.assertArrayEquals(expectedDsHiAtNodeInitialSp2, Arrays.copyOfRange(esDsHiAtNodeInitial1[1], 0, 256), 1E-13);
        Assert.assertEquals(2.820698371976182, logLik, 1e-6);
    }

    /*
     * Checks log-likelihood and other internal quantities®®
     * for bifurcating tree (height = 0.025), with 32 quantitative trait
     * bins, and dtMax = 0.005.
     *
     * Same as previous test, but because bifurcating tree is taller, less
     * than half of the tree height is in high resolution, and more than half
     * is in low resolution
     */
    @Test
    public void testPruneBifTree0025Height32Bins() {

        // indices for high to low res transfer
        int[] hiLoIdxs4Transfer = q32BifTree0025HeightDt0005.getHiLoIdxs4Transfer();

        // let us grab the E's and D's before integration
        double[][] esDsHiAtNode0, esDsHiAtNode1;
        esDsHiAtNode0 = q32BifTree0025HeightDt0005.getEsDsAtNode(0, false); // sp1
        esDsHiAtNode1 = q32BifTree0025HeightDt0005.getEsDsAtNode(1, false); // sp2

        /*
         * we are going to have a look at (make a deep copy of) the initial D's
         * in the assert below because they are the starting point of everything
         */
        double[][] esDsHiAtNodeInitial0 = new double[esDsHiAtNode0.length][esDsHiAtNode0[0].length];
        esDsHiAtNodeInitial0[0] = Arrays.copyOf(esDsHiAtNode0[0], esDsHiAtNode0[0].length); // E
        esDsHiAtNodeInitial0[1] = Arrays.copyOf(esDsHiAtNode0[1], esDsHiAtNode0[1].length); // D

        double[][] esDsHiAtNodeInitial1 = new double[esDsHiAtNode1.length][esDsHiAtNode1[0].length];
        esDsHiAtNodeInitial1[0] = Arrays.copyOf(esDsHiAtNode1[0], esDsHiAtNode1[0].length); // E
        esDsHiAtNodeInitial1[1] = Arrays.copyOf(esDsHiAtNode1[1], esDsHiAtNode1[1].length); // D

        double logLik = q32BifTree0025HeightDt0005.calculateLogP();

        int[] expectedHiLoIdxs4Transfer = new int[] { 6, 14, 22, 30, 38, 46, 54, 62, 70, 78, 86, 94, 102, 110, 118, 126, 134, 142, 150, 158, 166, 174, 182, 190, 198 };
        double[] expectedDsHiAtNodeInitialSp1 = { 0.308986942687903, 0, 0.350566009871371, 0, 0.396747087835907, 0, 0.447890605896858, 0, 0.504364398303888, 0, 0.566540754832024, 0, 0.634793036713348, 0, 0.709491856924629, 0, 0.791000831787404, 0, 0.879671919608544, 0, 0.975840371583655, 0, 1.07981933026376, 0, 1.19189412137632, 0, 1.31231629549353, 0, 1.44129748672436, 0, 1.57900316601788, 0, 1.72554637653023, 0, 1.88098154753774, 0, 2.04529849127956, 0, 2.21841669358911, 0, 2.40018001393971, 0, 2.59035191331783, 0, 2.7886113289072, 0, 2.9945493127149, 0, 3.20766654683839, 0, 3.42737184095615, 0, 3.65298170778044, 0, 3.88372109966426, 0, 4.11872537439949, 0, 4.35704354065101, 0, 4.59764281368466, 0, 4.83941449038287, 0, 5.08118112938378, 0, 5.3217049979751, 0, 5.55969772261993, 0, 5.79383105522966, 0, 6.02274864309609, 0, 6.24507866733522, 0, 6.45944719335828, 0, 6.66449205783599, 0, 6.85887710038768, 0, 7.04130653528599, 0, 7.21053924923296, 0, 7.36540280606647, 0, 7.50480693833876, 0, 7.62775630921048, 0, 7.73336233605698, 0, 7.82085387950912, 0, 7.88958661815778, 0, 7.93905094954024, 0, 7.96887828189528, 0, 7.97884560802865, 0, 7.96887828189528, 0, 7.93905094954024, 0, 7.88958661815778, 0, 7.82085387950912, 0, 7.73336233605698, 0, 7.62775630921048, 0, 7.50480693833876, 0, 7.36540280606647, 0, 7.21053924923296, 0, 7.04130653528599, 0, 6.85887710038768, 0, 6.66449205783599, 0, 6.45944719335828, 0, 6.24507866733522, 0, 6.02274864309609, 0, 5.79383105522965, 0, 5.55969772261993, 0, 5.32170499797509, 0, 5.08118112938378, 0, 4.83941449038287, 0, 4.59764281368466, 0, 4.35704354065101, 0, 4.11872537439949, 0, 3.88372109966426, 0, 3.65298170778044, 0, 3.42737184095615, 0, 3.20766654683839, 0, 2.9945493127149, 0, 2.7886113289072, 0, 2.59035191331783, 0, 2.40018001393971, 0, 2.21841669358911, 0, 2.04529849127956, 0, 1.88098154753774, 0, 1.72554637653023, 0, 1.57900316601788, 0, 1.44129748672436, 0, 1.31231629549353, 0, 1.19189412137632, 0, 1.07981933026376, 0, 0.975840371583655, 0, 0.879671919608544, 0, 0.791000831787404, 0, 0.709491856924628, 0, 0.634793036713349, 0, 0.566540754832024, 0, 0.504364398303888, 0, 0.447890605896858, 0, 0.396747087835907, 0, 0.350566009871371, 0, 0.308986942687903, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
        double[] expectedDsHiAtNodeInitialSp2 = { 0.000254946647636669, 0, 0.000319674822138109, 0, 0.000399835934138456, 0, 0.000498849425801072, 0, 0.000620828141157003, 0, 0.000770703934841743, 0, 0.000954372730824099, 0, 0.0011788613551308, 0, 0.00145251860604505, 0, 0.00178523314354266, 0, 0.00218868086879601, 0, 0.00267660451529771, 0, 0.00326512817532484, 0, 0.00397310942785545, 0, 0.00482253160451986, 0, 0.0058389385158292, 0, 0.00705191364734891, 0, 0.00849560541101504, 0, 0.0102092994868837, 0, 0.0122380386022755, 0, 0.0146332892566062, 0, 0.0174536539009152, 0, 0.0207656259132282, 0, 0.0246443833694604, 0, 0.0291746160933349, 0, 0.0344513787810736, 0, 0.0405809611459954, 0, 0.0476817640292969, 0, 0.0558851682975889, 0, 0.0653363811239983, 0, 0.0761952419644361, 0, 0.08863696823876, 0, 0.102852818461079, 0, 0.119050648395517, 0, 0.137455333812279, 0, 0.158309031659599, 0, 0.181871250031821, 0, 0.208418696288452, 0, 0.238244872152104, 0, 0.271659384673712, 0, 0.308986942687903, 0, 0.350566009871371, 0, 0.396747087835906, 0, 0.447890605896858, 0, 0.504364398303888, 0, 0.566540754832024, 0, 0.634793036713348, 0, 0.709491856924629, 0, 0.791000831787404, 0, 0.879671919608544, 0, 0.975840371583655, 0, 1.07981933026376, 0, 1.19189412137632, 0, 1.31231629549353, 0, 1.44129748672436, 0, 1.57900316601788, 0, 1.72554637653023, 0, 1.88098154753774, 0, 2.04529849127956, 0, 2.21841669358911, 0, 2.40018001393971, 0, 2.59035191331783, 0, 2.7886113289072, 0, 2.9945493127149, 0, 3.20766654683839, 0, 3.42737184095615, 0, 3.65298170778044, 0, 3.88372109966426, 0, 4.11872537439949, 0, 4.35704354065101, 0, 4.59764281368466, 0, 4.83941449038287, 0, 5.08118112938378, 0, 5.32170499797509, 0, 5.55969772261993, 0, 5.79383105522965, 0, 6.02274864309609, 0, 6.24507866733522, 0, 6.45944719335828, 0, 6.66449205783599, 0, 6.85887710038768, 0, 7.04130653528599, 0, 7.21053924923296, 0, 7.36540280606647, 0, 7.50480693833876, 0, 7.62775630921048, 0, 7.73336233605698, 0, 7.82085387950912, 0, 7.88958661815778, 0, 7.93905094954024, 0, 7.96887828189528, 0, 7.97884560802865, 0, 7.96887828189528, 0, 7.93905094954024, 0, 7.88958661815778, 0, 7.82085387950912, 0, 7.73336233605698, 0, 7.62775630921048, 0, 7.50480693833876, 0, 7.36540280606647, 0, 7.21053924923296, 0, 7.04130653528599, 0, 6.85887710038768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

        // if using JTransforms
        // int[] expectedHiLoIdxs4Transfer = new int[] { 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99 };
        // double[] expectedDsHiAtNodeInitialSp1 = new double[] { 0.308986942687903, 0.350566009871371, 0.396747087835907, 0.447890605896858, 0.504364398303888, 0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911, 2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839, 3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101, 4.59764281368466, 4.83941449038287, 5.08118112938378, 5.3217049979751, 5.55969772261993, 5.79383105522966, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599, 6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876, 7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024, 7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778, 7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647, 7.21053924923296, 7.04130653528599, 6.85887710038768, 6.66449205783599, 6.45944719335828, 6.24507866733522, 6.02274864309609, 5.79383105522965, 5.55969772261993, 5.32170499797509, 5.08118112938378, 4.83941449038287, 4.59764281368466, 4.35704354065101, 4.11872537439949, 3.88372109966426, 3.65298170778044, 3.42737184095615, 3.20766654683839, 2.9945493127149, 2.7886113289072, 2.59035191331783, 2.40018001393971, 2.21841669358911, 2.04529849127956, 1.88098154753774, 1.72554637653023, 1.57900316601788, 1.44129748672436, 1.31231629549353, 1.19189412137632, 1.07981933026376, 0.975840371583655, 0.879671919608544, 0.791000831787404, 0.709491856924628, 0.634793036713349, 0.566540754832024, 0.504364398303888, 0.447890605896858, 0.396747087835907, 0.350566009871371, 0.308986942687903 };
        // double[] expectedDsHiAtNodeInitialSp2 = new double[] { 0.000254946647636669, 0.000319674822138109, 0.000399835934138456, 0.000498849425801072, 0.000620828141157003, 0.000770703934841743, 0.000954372730824099, 0.0011788613551308, 0.00145251860604505, 0.00178523314354266, 0.00218868086879601, 0.00267660451529771, 0.00326512817532484, 0.00397310942785545, 0.00482253160451986, 0.0058389385158292, 0.00705191364734891, 0.00849560541101504, 0.0102092994868837, 0.0122380386022755, 0.0146332892566062, 0.0174536539009152, 0.0207656259132282, 0.0246443833694604, 0.0291746160933349, 0.0344513787810736, 0.0405809611459954, 0.0476817640292969, 0.0558851682975889, 0.0653363811239983, 0.0761952419644361, 0.08863696823876, 0.102852818461079, 0.119050648395517, 0.137455333812279, 0.158309031659599, 0.181871250031821, 0.208418696288452, 0.238244872152104, 0.271659384673712, 0.308986942687903, 0.350566009871371, 0.396747087835906, 0.447890605896858, 0.504364398303888, 0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911, 2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839, 3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101, 4.59764281368466, 4.83941449038287, 5.08118112938378, 5.32170499797509, 5.55969772261993, 5.79383105522965, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599, 6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876, 7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024, 7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778, 7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647, 7.21053924923296, 7.04130653528599, 6.85887710038768 };
        // Assert.assertArrayEquals(expectedHiLoIdxs4Transfer, hiLoIdxs4Transfer);
        // Assert.assertArrayEquals(expectedDsHiAtNodeInitialSp1, Arrays.copyOfRange(esDsHiAtNodeInitial0[1], 0, 103), 1E-13);
        // Assert.assertArrayEquals(expectedDsHiAtNodeInitialSp2, Arrays.copyOfRange(esDsHiAtNodeInitial1[1], 0, 103), 1E-13);

        Assert.assertArrayEquals(expectedHiLoIdxs4Transfer, hiLoIdxs4Transfer);
        Assert.assertArrayEquals(expectedDsHiAtNodeInitialSp1, Arrays.copyOfRange(esDsHiAtNodeInitial0[1], 0, 256), 1E-13);
        Assert.assertArrayEquals(expectedDsHiAtNodeInitialSp2, Arrays.copyOfRange(esDsHiAtNodeInitial1[1], 0, 256), 1E-13);
        Assert.assertEquals(2.816105371976182, logLik, 1e-6);
    }

    /*
     * Checks log-likelihood and other internal quantities
     * for tree with 3 species (height = 0.02), with 32 quantitative trait
     * bins, and dtMax = 0.005.
     */
    @Test
    public void testPruneThreeSpTree32Bins() {
        double logLik = q32ThreeSpTreeDt0005.calculateLogP();

        Assert.assertEquals(0.124798371976182, logLik, 1e-6);
    }

    /*
     * The full likelihood must support trait-dependent extinction, beyond merely populating its
     * rate arrays. The reference comes from validation/r_scripts/QuaSSENonconstantMuReference.R.
     * Corrected diversitree commit 5af1bf7bef1742512f3bc84beaf0e6115fc08e42; C/R and wider grids agree.
     */
    @Test
    public void testTraitDependentExtinctionLikelihoodAgainstDiversitree() {
        QuaSSEDistribution distribution = smallDistribution(0.0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
        // drift, diffusion, left, right, xLo minimum, log likelihood (Observed prior).
        double[][] references = {
                {0, .001, 3, 3, -.60, 2.7770738626454392},
                {0, .004, 5, 5, -.58, 2.7709935261296716},
                {-.1, .001, 3, 3, -.60, 2.7770706614684051},
                {.1, .001, 3, 3, -.60, 2.7770763257023154},
                {-1, .001, 2, 3, -.61, 2.7626916689257088},
                {1, .001, 3, 2, -.61, 2.7631155125293141}
        };
        for (double[] reference : references) {
            QuaSSEDistribution fresh = smallDistribution(reference[0], reference[1], .005,
                    "(sp1:0.1,sp2:0.1);", "Observed");
            distribution.driftInput.get().setValue(reference[0]);
            distribution.diffusionInput.get().setValue(reference[1]);
            Assert.assertEquals(reference[5], fresh.calculateLogP(), 1e-9);
            Assert.assertEquals(reference[5], distribution.calculateLogP(), 1e-9);
            Assert.assertEquals((int) reference[2], distribution.getNLeftFlanks(true));
            Assert.assertEquals((int) reference[3], distribution.getNRightFlanks(true));
            Assert.assertEquals(reference[4], distribution.getXMinLo(), 1e-14);
            Assert.assertEquals(reference[4] - .0075, distribution.getXMinHi(), 1e-14);
            int[] transfer = distribution.getHiLoIdxs4Transfer();
            for (int i = 0; i < transfer.length; i++)
                Assert.assertEquals(distribution.getX(true)[i], distribution.getX(false)[transfer[i]/2], 1e-13);
        }
    }

    // References constrain final likelihoods, not every E/D entry or resized native resources.
    // Compare FFTs, or native integration against Java with FFTW; obsolete with one implementation.
    @Test
    public void testFftImplementationsAndReinitialization() {
        QuaSSEDistribution selected = smallDistribution(0, .001, .005, "(sp1:0.02,sp2:0.02);", "Observed");
        QuaSSEDistribution reference = smallDistribution(0, .001, .005, "(sp1:0.02,sp2:0.02);", "Observed");
        reference.integrationBackendInput.setValue("java", reference);
        String referenceFFT = "native".equals(System.getProperty("test.quasse.integration")) ? "fftw" : "sst";
        reference.fftBackendInput.setValue(referenceFFT, reference);
        reference.initAndValidate();
        for (double drift : new double[] {0, 1, -1}) {
            selected.driftInput.get().setValue(drift);
            reference.driftInput.get().setValue(drift);
            Assert.assertEquals(reference.calculateLogP(), selected.calculateLogP(), 1e-10);
            for (boolean low : new boolean[] {true, false}) {
                double[][][] expected = reference.getEsDs(low), actual = selected.getEsDs(low);
                for (int node = 0; node < expected.length; node++)
                    for (int dim = 0; dim < expected[node].length; dim++)
                        Assert.assertArrayEquals(expected[node][dim], actual[node][dim], 1e-10);
            }
        }
        // Reinitialization changes transform lengths, unlike parameter-dependent padding changes.
        for (QuaSSEDistribution distribution : new QuaSSEDistribution[] {reference, selected}) {
            distribution.nXbinsInput.get().setValue(256);
            distribution.highLowRatioInput.get().setValue(3);
            distribution.initAndValidate();
        }
        Assert.assertEquals(reference.calculateLogP(), selected.calculateLogP(), 1e-10);
        // Switching back releases native storage instead of retaining closed resources.
        selected.integrationBackendInput.setValue("java", selected);
        selected.initAndValidate();
        Assert.assertEquals(reference.calculateLogP(), selected.calculateLogP(), 1e-10);
        selected.integrationBackendInput.setValue("invalid", selected);
        Assert.assertThrows(IllegalArgumentException.class, selected::initAndValidate);
        selected.integrationBackendInput.setValue("native", selected);
        selected.fftBackendInput.setValue("sst", selected);
        Assert.assertThrows(IllegalArgumentException.class, selected::initAndValidate);
    }

    // Build the small reference model with independent inputs, so comparisons cannot share caches.
    private QuaSSEDistribution smallDistribution(double drift, double diffusion, double tc,
                                                 String newick, String prior) {
        Tree tree = new TreeParser(newick, false, false, true, 0);

        RealParameter traits = new RealParameter();
        traits.initByName("value", Arrays.asList(0.0, 0.1), "keys", "sp1 sp2");

        ConstantLinkFn speciation = new ConstantLinkFn();
        speciation.initByName("yV", new RealParameter(new Double[] { 0.15 }));

        LogisticFunction extinction = new LogisticFunction();
        extinction.initByName(
                "curveYBaseValue", new RealParameter(new Double[] { 0.01 }),
                "curveMaxY", new RealParameter(new Double[] { 0.08 }),
                "sigmoidMidpoint", new RealParameter(new Double[] { 0.05 }),
                "logisticGrowthRate", new RealParameter(new Double[] { 20.0 }));

        NormalCenteredAtObservedLinkFn tipLink = new NormalCenteredAtObservedLinkFn();
        tipLink.initByName("quTraits", traits,
                "sdNormalQuTrValue", new RealParameter(new Double[] { 0.05 }));

        QuaSSEDistribution distribution = newDistribution();
        distribution.initByName(
                "dtMax", new RealParameter(new Double[] { 0.005 }),
                "dynDt", new BooleanParameter(new Boolean[] { true }),
                "tc", new RealParameter(new Double[] { tc }),
                "nX", new IntegerParameter(new Integer[] { 128 }),
                "dX", new RealParameter(new Double[] { 0.01 }),
                "xMid", new RealParameter(new Double[] { 0.0 }),
                "flankWidthScaler", new RealParameter(new Double[] { 10.0 }),
                "hiLoRatio", new IntegerParameter(new Integer[] { 4 }),
                "drift", new RealParameter(new Double[] { drift }),
                "diffusion", new RealParameter(new Double[] { diffusion }),
                "q2mLambda", speciation,
                "q2mMu", extinction,
                "tree", tree,
                "q2d", tipLink,
                "priorProbAtRootType", prior);

        return distribution;
    }

    // Rejecting a proposal leaves unsnapshotted grids and link scalars at proposed values. A later
    // change to a different parameter must recover clean restored inputs, not merely the logP snapshot.
    @Test
    public void testGridRejectionLifecycle() {
        QuaSSEDistribution distribution = smallDistribution(0, .001, .1, "(sp1:0.09,sp2:0.09);", "Observed");
        RealParameter drift = distribution.driftInput.get();
        RealParameter diffusion = distribution.diffusionInput.get();
        RealParameter lambda = ((ConstantLinkFn) distribution.q2mLambdaInput.get()).yValueInput.get();
        RealParameter mu = ((LogisticFunction) distribution.q2mMuInput.get()).curveYBaseValueInput.get();
        NormalCenteredAtObservedLinkFn tip = (NormalCenteredAtObservedLinkFn) distribution.q2dInput.get();
        RealParameter sd = tip.sdNormalQuTrValueInput.get();
        RealParameter traits = tip.quTraitsInput.get();
        Tree tree = distribution.treeInput.get();
        State state = new State();
        state.initByName("stateNode", Arrays.asList(drift, diffusion, lambda, mu, sd, traits, tree));
        state.initialise();
        CompoundDistribution posterior = new CompoundDistribution();
        posterior.initByName("distribution", distribution);
        state.setPosterior(posterior);
        double accepted = state.robustlyCalcPosterior(posterior);
        for (double proposedDiffusion : new double[] {.004, .0011, 1e100}) {
            state.store(1);
            diffusion.setValue(proposedDiffusion);
            drift.setValue(1.0);
            mu.setValue(.02);
            sd.setValue(.07);
            traits.setValue(0, .01);
            tree.getRoot().setHeight(.11);
            state.storeCalculationNodes();
            state.checkCalculationNodesDirtiness();
            double proposed = posterior.calculateLogP();
            if (proposedDiffusion < 1) {
                QuaSSEDistribution fresh = smallDistribution(1, proposedDiffusion, .1,
                        "(sp1:0.11,sp2:0.11);", "Observed");
                ((ConstantLinkFn) fresh.q2mLambdaInput.get()).yValueInput.get().setValue(lambda.getValue());
                ((LogisticFunction) fresh.q2mMuInput.get()).curveYBaseValueInput.get().setValue(.02);
                NormalCenteredAtObservedLinkFn freshTip = (NormalCenteredAtObservedLinkFn) fresh.q2dInput.get();
                freshTip.sdNormalQuTrValueInput.get().setValue(.07);
                freshTip.quTraitsInput.get().setValue(0, .01);
                Assert.assertEquals(fresh.calculateLogP(), proposed, 1e-12);
            } else {
                Assert.assertEquals(Double.NEGATIVE_INFINITY, proposed, 0.0);
            }
            state.restore();
            state.restoreCalculationNodes();
            state.setEverythingDirty(false);
            Assert.assertEquals(accepted, distribution.getCurrentLogP(), 0.0);

            // Do not recalculate after rejection: only lambda is dirty in this next proposal.
            state.store(2);
            lambda.setValue(lambda.getValue() + .01);
            state.storeCalculationNodes();
            state.checkCalculationNodesDirtiness();
            QuaSSEDistribution fresh = smallDistribution(0, .001, .1, "(sp1:0.09,sp2:0.09);", "Observed");
            ((ConstantLinkFn) fresh.q2mLambdaInput.get()).yValueInput.get().setValue(lambda.getValue());
            accepted = posterior.calculateLogP();
            Assert.assertEquals(fresh.calculateLogP(), accepted, 1e-12);
            state.acceptCalculationNodes();
            state.setEverythingDirty(false);
        }
    }

    // Kernel validity is per resolution and representation, not a global dt/dirty hint. Compare
    // non-forced calls with forced construction through raw/FFT and parameter A-B-A transitions.
    @Test
    public void testKernelValidity() {
        QuaSSEDistribution distribution = smallDistribution(0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
        QuaSSEDistribution fresh = smallDistribution(0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
        for (double[] parameters : new double[][] {{0, .001}, {.1, .001}, {0, .001}, {0, .0011}, {0, .001}}) {
            distribution.driftInput.get().setValue(parameters[0]);
            fresh.driftInput.get().setValue(parameters[0]);
            distribution.diffusionInput.get().setValue(parameters[1]);
            fresh.diffusionInput.get().setValue(parameters[1]);
            for (boolean jtransforms : new boolean[] {false, true, false}) {
                for (boolean doFFT : new boolean[] {false, true, false, true}) {
                    for (double dt : new double[] {.005, .002, .005}) {
                        for (boolean low : new boolean[] {true, false, true}) {
                            distribution.populatefY(dt, false, doFFT, low, jtransforms);
                            fresh.populatefY(dt, true, doFFT, low, jtransforms);
                            Assert.assertArrayEquals(fresh.getfY(low), distribution.getfY(low), 0.0);
                            Assert.assertArrayEquals(fresh.getfftFY(low), distribution.getfftFY(low), 0.0);
                        }
                    }
                }
            }
        }
    }

    // A reused distribution must match a fresh one across padding changes, including equal-size
    // grids with exchanged left/right support. Invalid proposals must not mutate the usable grid.
    @Test
    public void testDynamicGridAgainstFreshCalculation() {
        QuaSSEDistribution distribution = smallDistribution(0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
        for (double[] parameters : new double[][] {{0, .004}, {1, .001}, {-1, .001}, {.1, .001}, {0, .001}}) {
            distribution.driftInput.get().setValue(parameters[0]);
            distribution.diffusionInput.get().setValue(parameters[1]);
            QuaSSEDistribution fresh = smallDistribution(parameters[0], parameters[1], .005,
                    "(sp1:0.1,sp2:0.1);", "Observed");
            Assert.assertEquals(fresh.calculateLogP(), distribution.calculateLogP(), 1e-12);
            Assert.assertEquals(fresh.getNLeftFlanks(true), distribution.getNLeftFlanks(true));
            Assert.assertEquals(fresh.getNRightFlanks(true), distribution.getNRightFlanks(true));
            Assert.assertArrayEquals(fresh.getX(true), distribution.getX(true), 0.0);
        }
        double[] ruler = distribution.getX(true);
        for (double invalid : new double[] {Double.NaN, Double.POSITIVE_INFINITY, 1e100}) {
            distribution.driftInput.get().setValue(invalid);
            Assert.assertEquals(Double.NEGATIVE_INFINITY, distribution.calculateLogP(), 0.0);
            Assert.assertSame(ruler, distribution.getX(true));
        }
        distribution.driftInput.get().setValue(0.0);
        for (double invalid : new double[] {-1, 0, Double.NaN, Double.POSITIVE_INFINITY, 1e100}) {
            distribution.diffusionInput.get().setValue(invalid);
            Assert.assertEquals(Double.NEGATIVE_INFINITY, distribution.calculateLogP(), 0.0);
            Assert.assertSame(ruler, distribution.getX(true));
        }
        distribution.diffusionInput.get().setValue(.001);
        Assert.assertEquals(2.7770738626454392, distribution.calculateLogP(), 1e-9);
        distribution.dtMaxInput.get().setValue(.01);
        Assert.assertThrows(IllegalArgumentException.class, distribution::calculateLogP);
    }

    // Padding must cover short-step support, not only dtMax, and allow valid one-bin strips.
    // Small-drift reference likelihoods cannot distinguish the interior maximum of the opposing side.
    @Test
    public void testAllTimestepPadding() {
        for (double drift : new double[] {-2, 2}) {
            QuaSSEDistribution distribution = smallDistribution(0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
            distribution.dtMaxInput.get().setValue(1.0);
            distribution.dXBinInput.get().setValue(.1);
            distribution.flankWidthScalerInput.get().setValue(3.0);
            distribution.driftInput.get().setValue(drift);
            distribution.diffusionInput.get().setValue(1.0);
            distribution.initAndValidate();
            Assert.assertEquals(drift > 0 ? 50 : 12, distribution.getNLeftFlanks(true));
            Assert.assertEquals(drift > 0 ? 12 : 50, distribution.getNRightFlanks(true));
        }
        QuaSSEDistribution one = smallDistribution(0, .00001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
        Assert.assertEquals(1, one.getNLeftFlanks(true));
        Assert.assertTrue(Double.isFinite(one.calculateLogP()));
    }

    // Preserve fixed-dt truncation and normalization even when a nonempty segment takes no steps.
    // Full likelihood references use dynamic dt; replace this check when fixed-dt behavior is repaired.
    @Test
    public void testFixedStepSegments() {
        for (boolean low : new boolean[]{true, false}) {
            for (double length : new double[]{0, .002, .012}) {
                QuaSSEDistribution selected = smallDistribution(0, .001, .005,
                        "(sp1:0.02,sp2:0.02);", "Observed");
                QuaSSEDistribution reference = smallDistribution(0, .001, .005,
                        "(sp1:0.02,sp2:0.02);", "Observed");
                reference.integrationBackendInput.setValue("java", reference);
                reference.initAndValidate();
                double[][] expected = reference.getEsDs(low)[0], actual = selected.getEsDs(low)[0];
                // Make normalization observable independently of the initial tip density's integral.
                for (int i = 0; i < expected[1].length; i += 2) {
                    expected[1][i] *= 3;
                    actual[1][i] *= 3;
                }
                if (length != 0) {
                    reference.populatefY(.005, false, true, low, false);
                    for (int step = 0; step < Math.floor(length / .005); ++step)
                        reference.doIntegrateInPlace(0, .005, low);
                    reference.normalizeDs(0, low, false);
                }
                selected.integrateLength(0, actual, selected.getScratchAtNode(0, low),
                        length, false, .005, low, false, false);
                for (int d = 0; d < expected.length; ++d)
                    Assert.assertArrayEquals(expected[d], actual[d], 1e-10);
            }
        }
    }

    // Root crossings must resize the prior, with time == tc consistently coarse and empty segments
    // acting as identities. Fixed-tree reference values never exercise these transitions.
    @Test
    public void testResolutionTransitions() {
        for (String prior : Arrays.asList("Flat", "Observed")) {
            QuaSSEDistribution distribution = smallDistribution(0, .001, .1, "(sp1:0.1,sp2:0.1);", prior);
            for (double height : new double[] {.09, .1, .11, .1, .09}) {
                distribution.treeInput.get().getRoot().setHeight(height);
                QuaSSEDistribution fresh = smallDistribution(0, .001, .1,
                        "(sp1:" + height + ",sp2:" + height + ");", prior);
                double expected = fresh.calculateLogP();
                Assert.assertTrue(Double.isFinite(expected));
                Assert.assertEquals(expected, distribution.calculateLogP(), 1e-12);
            }
            QuaSSEDistribution zero = smallDistribution(0, .001, 0, "(sp1:0.0,sp2:0.0);", prior);
            Assert.assertTrue(Double.isFinite(zero.calculateLogP()));
        }
    }

    // A supplied array has no defined moving-grid interpretation; fail explicitly rather than silently
    // using a built-in prior. Replace this check when coordinate-aware supplied priors are implemented.
    @Test
    public void testSuppliedRootPriorRejected() {
        QuaSSEDistribution distribution = smallDistribution(0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
        distribution.priorProbsAtRootInput.setValue(new RealParameter(new Double[] {1.0}), distribution);
        IllegalArgumentException error = Assert.assertThrows(IllegalArgumentException.class, distribution::initAndValidate);
        Assert.assertTrue(error.getMessage().contains("givenPriorProbsAtRoot"));
    }

    /*
     * A full likelihood call must depend only on current inputs, not arrays mutated by an earlier
     * traversal or a shared link's cache. Independent node/rate caches could replace this check.
     */
    @Test
    public void testRepeatedLikelihoodCalculation() {
        double first = q32ThreeSpTreeDt0005.calculateLogP();
        double second = q32ThreeSpTreeDt0005.calculateLogP();

        Assert.assertEquals(first, second, 1e-12);

        QuaSSEDistribution a = smallDistribution(0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
        QuaSSEDistribution b = smallDistribution(0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
        b.q2mMuInput.setValue(a.q2mMuInput.get(), b);
        b.initAndValidate();
        RealParameter mu = ((LogisticFunction) a.q2mMuInput.get()).curveYBaseValueInput.get();
        for (double value : new double[] {.02, .03, .02}) {
            mu.setValue(value);
            mu.setEverythingDirty(false);
            QuaSSEDistribution fresh = smallDistribution(0, .001, .005, "(sp1:0.1,sp2:0.1);", "Observed");
            ((LogisticFunction) fresh.q2mMuInput.get()).curveYBaseValueInput.get().setValue(value);
            double expected = fresh.calculateLogP();
            Assert.assertEquals(expected, a.calculateLogP(), 1e-12);
            Assert.assertEquals(expected, b.calculateLogP(), 1e-12);
        }
    }

    /*
     * BEAST acceptance/rejection must preserve both link values and the cached likelihood scalar.
     * One-shot checks and manual parameter round trips do not exercise dependency tracking.
     */
    @Test
    public void testLikelihoodParameterRoundTrip() {
        RealParameter deathRate = new RealParameter(new Double[] { 0.03 });
        ConstantLinkFn deathLink = new ConstantLinkFn();
        deathLink.initByName("yV", deathRate);
        QuaSSEDistribution distribution = newDistribution();
        distribution.initByName("dtMax", dt0005Rp, "dynDt", dynDtbpTrue,
                "tc", tc0005Rp,
                "nX", nXbins32Ip, "dX", dxBin001Rp, "xMid", xMid00Rp,
                "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", deathLink,
                "tree", threeSpTreeHeight002,
                "q2d", nfn3Sp,
                "priorProbAtRootType", rootPriorType);

        RealParameter unrelated = new RealParameter();
        unrelated.initByName("value", "0.0");
        Normal normal = new Normal();
        normal.initAndValidate();
        Prior prior = new Prior();
        prior.initByName("x", unrelated, "distr", normal);
        CompoundDistribution posterior = new CompoundDistribution();
        posterior.initByName("distribution", Arrays.asList(distribution, prior));
        State state = new State();
        state.initByName("stateNode", Arrays.asList(deathRate, unrelated));
        state.initialise();
        state.setPosterior(posterior);
        double initial = state.robustlyCalcPosterior(posterior);
        double initialLikelihood = distribution.getCurrentLogP();

        state.store(1);
        deathRate.setValue(0.06);
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
        double changed = posterior.calculateLogP();
        Assert.assertTrue(distribution.isDirtyCalculation());
        state.restore();
        state.restoreCalculationNodes();
        state.setEverythingDirty(false);
        Assert.assertEquals("restored scalar", initialLikelihood, distribution.getCurrentLogP(), 0.0);
        Assert.assertEquals("recalculated restored inputs", initialLikelihood, distribution.calculateLogP(), 1e-12);

        Assert.assertNotEquals(initial, changed, 1e-8);
        state.store(2);
        unrelated.setValue(1.0);
        state.storeCalculationNodes();
        state.checkCalculationNodesDirtiness();
        Assert.assertFalse(distribution.isDirtyCalculation());
        Assert.assertEquals("unchanged likelihood reuse", initial - .5, posterior.calculateLogP(), 1e-12);
        state.acceptCalculationNodes();
        state.setEverythingDirty(false);
    }

    /*
     * Checks log-likelihood and other internal quantities
     * for tree with 15 species, with 1024 quantitative trait
     * bins.
     *
     * This is shows decisively that this implementation
     * is correct.
     */
    @Test
    public void testPruneFifteenSpTree1024Bins() {
        double logLik = q32FifteenSp.calculateLogP();

        // Original rounded R reference, corrected by -2*log(dx) with this model's dx=0.01027592.
        Assert.assertEquals(-61.27245 - 2*Math.log(0.01027592), logLik, 1e-5);
    }

    /*
     * Checks that the number of quantitative character bins is a power of 2.
     * If it isn't, an exception must be thrown
     */
    @Test(expected = RuntimeException.class)
    public void testPowerOf2() {

        QuaSSEDistribution q48Dt001 = newDistribution();
        q48Dt001.initByName("dtMax", dt001Rp, "dynDt", dynDtbpTrue,
                "tc", tc100Rp,
                "nX", nXbins48Ip, "dX", dxBin001Rp, "xMid", xMid00Rp, "flankWidthScaler", flankWidthScaler10Rp, "hiLoRatio", hiLoRatioIp,
                "drift", driftRp, "diffusion", diffusionRp0001,
                "q2mLambda", lfn, "q2mMu", cfn,
                "tree", bifTreeHeight001,
                "q2d", nfn2Sp,
                "priorProbAtRootType", rootPriorType);
    }
}
