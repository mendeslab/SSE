package SSE;

import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.evolution.tree.Node;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class QuaSSEDistribution extends QuaSSEProcess {

    final public Input<LinkFn> q2mLambdaInput = new Input<>("q2mLambda", "Function converting quantitative trait into lambda parameter.", Input.Validate.REQUIRED);
    final public Input<LinkFn> q2mMuInput = new Input<>("q2mMu", "Function converting quantitative trait into mu parameter.", Input.Validate.REQUIRED);
    final public Input<LinkFn> q2dInput = new Input<>("q2d", "Function converting quantitative trait into initial D values.", Input.Validate.REQUIRED);

    // dimension state
    protected int nDimensions = 2;
    protected int nDimensionsE = 1;
    protected int nDimensionsD = 1;

    // macroevo state (for propagate in t)
    protected double[] birthRatesLo, deathRatesLo, birthRatesHi, deathRatesHi; // macroevol parameters
    protected LinkFn q2mLambda, q2mMu, q2d;

    // state that matters directly for calculateLogP
    protected double[][][] esDsLo, esDsHi; // first dimension are nodes, second is Es and Ds, third is each E (or D) along X ruler
    protected double[][][] fftBufferEsDsLo, fftBufferEsDsHi; // same as above, to be used in propagate in X (during convolution)
    protected double[][][] scratchLo, scratchHi;
    protected double[] logNormalizationFactors;

    @Override
    public void initAndValidate() {

        super.initAndValidate(); // read in all dimension-related stuff, populates fYLo and fYHi

        checkDimensions();

        int nNodes = tree.getNodeCount();

        logNormalizationFactors = new double[nNodes];

        q2mLambda = q2mLambdaInput.get();
        birthRatesLo = new double[nUsefulXbinsLo];
        birthRatesHi = new double[nUsefulXbinsHi];

        q2mMu = q2mMuInput.get();
        deathRatesLo = new double[nUsefulXbinsLo];
        deathRatesHi = new double[nUsefulXbinsHi];

        q2d = q2dInput.get();

        scratchLo = new double[nNodes][nDimensions][2 * nXbinsLo]; // for real and complex part after FFT using JTransforms
        scratchHi = new double[nNodes][nDimensions][2 * nXbinsHi]; // for real and complex part after FFT using JTransforms

        populateMacroevolParams(true);

        int nDimensionsFFT = nDimensions;
        initializeEsDs(nNodes, nDimensionsFFT, nXbinsLo, nXbinsHi);

        populateTipsEsDs(nDimensionsFFT, nXbinsHi, true, false);

    }

    private void checkDimensions() {
        if (!((nXbinsLo & nXbinsLo-1) == 0)) throw new RuntimeException("Number of quantitative character bins must be a power of 2. It was " + nXbinsLo + ". Exiting...");
        if ((nXbinsLo - nLeftNRightFlanksLo[1] - (nLeftNRightFlanksLo[0] + nLeftNRightFlanksLo[1])) <= nLeftNRightFlanksLo[0]) throw new RuntimeException("Left and right flanking bins were too many, leaving no useful bins between. Exiting...");
        if (nLeftNRightFlanksLo[0] == 0 || nLeftNRightFlanksLo[0] == 1) throw new RuntimeException("Left and right flank bins were set to 0 (invalid). Exiting...");
    }

    @Override
    public void populateMacroevolParams(boolean ignoreRefresh) {
        birthRatesLo = q2mLambda.getY(xLo, birthRatesLo, ignoreRefresh);
        birthRatesHi = q2mLambda.getY(xHi, birthRatesHi, ignoreRefresh);
        deathRatesLo = q2mMu.getY(xLo, deathRatesLo, ignoreRefresh);
        deathRatesHi = q2mMu.getY(xHi, deathRatesHi, ignoreRefresh);
    }

    @Override
    public void initializeEsDs(int nNodes, int nDimensionsFFT, int nXbinsLo, int nXbinsHi) {
        esDsLo = new double[nNodes][nDimensionsFFT][2 * nXbinsLo]; // 2 * for real and complex part after FFT
        esDsHi = new double[nNodes][nDimensionsFFT][2 * nXbinsHi];

        fftBufferEsDsLo = new double[nNodes][nDimensionsFFT][2 * nXbinsLo]; // 2 * for real and complex part after FFT
        fftBufferEsDsHi = new double[nNodes][nDimensionsFFT][2 * nXbinsHi];
    }

    @Override
    public void populateTipsEsDs(int nDimensionsFFT, int nXbinsHi, boolean ignoreRefresh, boolean jtransforms) {
        for (Node tip: tree.getExternalNodes()) {
            String tipName = tip.getID();
            int nodeIdx = tip.getNr();

            for (int i=0; i < nDimensionsFFT; i++) {
                // E's
                if (i == 0) {
                    // high res
                    for (int j=0; j < nXbinsHi; j++) esDsHi[nodeIdx][i][j] = 0.0; // E's = 1.0 - sampling.f

                    // low res (just for debugging)
                    for (int j=0; j < nXbinsLo; j++) esDsLo[nodeIdx][i][j] = 0.0; // E's = 1.0 - sampling.f
                }

                // D's (initialize both just in case tip is non-contemporaneous and t > tc)
                else {
                    // high res
                    esDsHi[nodeIdx][i] = q2d.getY(xHi, esDsHi[nodeIdx][i], nLeftNRightFlanksHi, tipName, ignoreRefresh);

                    // low res (just for debugging)
                    esDsLo[nodeIdx][i] = q2d.getY(xLo, esDsLo[nodeIdx][i], nLeftNRightFlanksLo, tipName, ignoreRefresh);

                    if (!jtransforms) {
                        SSEUtils.everyOtherExpandInPlace(esDsHi[nodeIdx][i]);
                        SSEUtils.everyOtherExpandInPlace(esDsLo[nodeIdx][i]);
                    }
                }
            }

            double dsSumLo = 0.0;
            double dsSumHi = 0.0;
            for (double d: esDsHi[nodeIdx][1]) dsSumHi += d;
            for (double d: esDsLo[nodeIdx][1]) dsSumLo += d;
            if (dsSumHi == 0.0 || dsSumLo == 0.0)
                throw new RuntimeException("ERROR: Initial D's were all 0.0 for species " + tipName + ". " +
                        "This can happen when the trait value standard deviation" +
                        "(of the normal centered at the observed trait value)" +
                        "is too small relative to trait value bin sizes." +
                        " Try increasing the standard deviation.");
        }
    }

    @Override
    public void startRecursionAtRootNode(Node rootNode, boolean forceRecalcKernel, boolean jtransforms) {
        int rootIdx = rootNode.getNr();
        double rootHeight = rootNode.getHeight();

        // start recursion
        processInternalNode(rootNode, forceRecalcKernel, jtransforms);

        // we're done pruning, let's deal with the prior probs at root now
        double[][] esDsAtRoot;
        double dxAtRightRes;
        int nXBinsRightRes, nUsefulXBinsRightRes;
        if (rootHeight > tc) {
            esDsAtRoot = esDsLo[rootIdx];
            dxAtRightRes = dXbin;
            nXBinsRightRes = nXbinsLo;
            nUsefulXBinsRightRes = nUsefulXbinsLo;
            got2LowRes = true;
        } else {
            esDsAtRoot = esDsHi[rootIdx];
            dxAtRightRes = dXbin / hiLoRatio;
            nXBinsRightRes = nXbinsHi;
            nUsefulXBinsRightRes = nUsefulXbinsHi;
            got2LowRes = false;
        }

        if (priorProbsAtRoot == null) priorProbsAtRoot = new double[nUsefulXBinsRightRes]; // if root prior probs not given, initialize it!

        populatePriorProbAtRoot(esDsAtRoot[1], dxAtRightRes, nXBinsRightRes, nUsefulXBinsRightRes, rootPriorType, jtransforms); // ok, prior probs are set!
    }

    @Override
    public void processInternalNode(Node aNode, boolean forceRecalcKernel, boolean jtransforms) {

        // debugging
        // System.out.println("\nDoing node " + aNode.getID());

        // recur if internal node or sampled ancestor
        if (!aNode.isLeaf()) {
            for (Node childNode: aNode.getChildren()) {
                processInternalNode(childNode, forceRecalcKernel, jtransforms); // recur
            }

            // after recursion, prepare initial conditions for integrating this node
            mergeChildrenNodes(aNode, jtransforms);
        }

        // unless we're at the root, now we integrate this branch
        if (!aNode.isRoot()) processBranch(aNode, forceRecalcKernel, jtransforms);
    }

    // main branch function
    @Override
    public void processBranch(Node aNode, boolean forceRecalcKernel, boolean jtransforms) {
        int nodeIdx = aNode.getNr();

        // dealing with branch lengths
        double startTime = aNode.getHeight(); // we're going backwards in time, toward the root
        double endTime = aNode.getParent().getHeight();
        double branchLength2Integrate = endTime - startTime;

        /*
         * Sorting out E's and D's prior to integration
         * (1) Grabbing esDs and scratch depending on resolution
         * (2) Calculating normalization factor
         * (3) Storing normalization factor to unnormalize at the very end of pruning
         * (4) Normalize DsThis is done only once if node is already at low, rather than
         * at every loop of the while block below
         */
        double[][] esDsAtNode, scratchAtNode;

        // Debugging
        // System.out.println("Log-normalization factor at start of branch = " + logNormalizationFactors[nodeIdx]);
//        // option 0: (rare) branch starts at tc (all low-res)
//        if (startTime == tc) {
//            for (int ithDim=0; ithDim<nDimensions; ithDim++)
//                SSEUtils.hiToLoTransferInPlace(esDsHi[nodeIdx][ithDim], esDsLo[nodeIdx][ithDim], hiLoIdxs4Transfer);
//            esDsAtNode = esDsLo[nodeIdx];
//            scratchAtNode = scratchLo[nodeIdx];
//        }

        // option 1: entire branch might already be > tc (all low-res)
        // note that here the high to low res transfer should have already happened
        if (startTime > tc) {
            esDsAtNode = esDsLo[nodeIdx];
            scratchAtNode = scratchLo[nodeIdx];

            // debugging
            // System.out.println("OPTION 1: E's size should be 2 * low-res = " + esDsAtNode[0].length);
            // System.out.println("Unnormalized D's at low-res = " + Arrays.toString(esDsAtNode[1]));

            // normalize and record normalization factor before integration
            // logNormalizationFactors[nodeIdx] = normalizeDs(esDsAtNode[1], dXbin); // normalize D's and returns factor, which we record
            logNormalizationFactors[nodeIdx] = normalizeDs(nodeIdx, true, jtransforms); // normalize D's and returns factor, which we record

            // now integrate whole branch
            integrateLength(nodeIdx, esDsAtNode, scratchAtNode, branchLength2Integrate, dynamicallyAdjustDt, dtMax, true, forceRecalcKernel, jtransforms);
        }

        // option 2: entire branch < tc (all high-res)
        else if ((startTime + branchLength2Integrate) < tc) {
            esDsAtNode = esDsHi[nodeIdx];
            scratchAtNode = scratchHi[nodeIdx];

            // debugging
            // System.out.println("OPTION 2: E's size should be 4 * low-res = " + esDsAtNode[0].length);
            // System.out.println("Unnormalized D's at low-res = " + Arrays.toString(esDsAtNode[1]));

            // normalize and record normalization factor before integration
            // logNormalizationFactors[nodeIdx] = normalizeDs(esDsAtNode[1], dXbin/hiLoRatio); // normalize D's and returns factor, which we record
            logNormalizationFactors[nodeIdx] = normalizeDs(nodeIdx, false, jtransforms); // normalize D's and returns factor, which we record

            // now integrate whole branch (again, normalization and adding to logNormalizationFactors inside)
            integrateLength(nodeIdx, esDsAtNode, scratchAtNode, branchLength2Integrate, dynamicallyAdjustDt, dtMax, false, forceRecalcKernel, jtransforms);
        }

        // option 3: tc happens inside branch (tip-end part in high-res, root-end part in low-res)
        else {
            esDsAtNode = esDsHi[nodeIdx]; // copying address, so whatever happens to esDsAtNode, happens to esDsHi[nodeIdx]
            scratchAtNode = scratchHi[nodeIdx];

            // debugging
            // System.out.println("OPTION 3: E's size should be 2 * 4 * low-res = " + esDsAtNode[0].length);
            // System.out.println("Unnormalized D's at high-res = " + Arrays.toString(esDsAtNode[1]));

            // normalize and record normalization factor before integration
            // logNormalizationFactors[nodeIdx] = normalizeDs(esDsAtNode[1], dXbin/hiLoRatio); // normalize D's and returns factor, which we record

            // System.out.println("calling normalize before integration");
            logNormalizationFactors[nodeIdx] = normalizeDs(nodeIdx, false, jtransforms); // normalize D's and returns factor, which we record

            // high res part
            double lenHi = tc - startTime;

            // debugging
            // System.out.println("high-res part, lenHi = " + lenHi);
            // System.out.println("calling normalize inside integration in high res");
            integrateLength(nodeIdx, esDsAtNode, scratchAtNode, lenHi, dynamicallyAdjustDt, dtMax, false, forceRecalcKernel, jtransforms); // normalization and adding to logNormalizationFactors inside

            // normalize and record normalization factor before integration at low res
            // logNormalizationFactors[nodeIdx] += normalizeDs(esDsAtNode[1], dXbin/hiLoRatio);

            // debugging
            // System.out.println("Log-normalization factors after high-res part = " + logNormalizationFactors[nodeIdx]);

            // transferring high-res esDs to low-res esDs
            for (int ithDim=0; ithDim<nDimensions; ithDim++)
                SSEUtils.hiToLoTransferInPlace(esDsHi[nodeIdx][ithDim], esDsLo[nodeIdx][ithDim], hiLoIdxs4Transfer, jtransforms);

            esDsAtNode = esDsLo[nodeIdx]; // copying address, so whatever happens to esDsAtNode, happens to esDsHi[nodeIdx]
            scratchAtNode = scratchLo[nodeIdx];

            // low res part
            double lenLo = branchLength2Integrate - lenHi;

            // debugging
            // System.out.println("low-res part, lenLo = " + lenLo);
            // System.out.println("calling normalize after integration in low res");
            integrateLength(nodeIdx, esDsAtNode, scratchAtNode, lenLo, dynamicallyAdjustDt, dtMax, true, forceRecalcKernel, jtransforms);
        }
    }

    @Override
    public void integrateLength(
    		int nodeIdx,
    		double[][] esDsAtNode,
    		double[][] scratchAtNode,
    		double aLength,
    		boolean dynamicallyAdjust,
    		double maxDt,
    		boolean lowRes,
    		boolean forceRecalcKernel,
    		boolean jtransforms) {

        // option of jtransforms, and jni

        // dealing with dt
        double dt, nIntervals;
        boolean dtChanged = false;
        // double nonIntegratedDt;

        if (dynamicallyAdjust) {
            nIntervals = Math.ceil(aLength / maxDt);
            dt = aLength / nIntervals; // dynamically adjusted dt
            dtChanged = true;
        } else {
            nIntervals = Math.floor(aLength / maxDt);
            dt = maxDt;
            // nonIntegratedDt = aLength % maxDt; // maybe later we do something with this when we fix dt at dtMax
        }

        // updating fY if necessary
        // remove next two later
        populatefY(dt, forceRecalcKernel, dtChanged, true, lowRes, false);

        // integrating!
        for (int i=0; i<nIntervals; i++) {
            doIntegrateInPlace(nodeIdx, dt, lowRes);
        }

        // normalization that happens in make.pde.quasse.fftR
        // System.out.println("Prior to normalization inside integrateLength (lowRes=" + lowRes + ") = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("After normalization inside integrateLength (lowRes=" + lowRes + ") = " + Arrays.toString(esDsAtNode[1]));

        // if (lowRes) logNormalizationFactors[nodeIdx] += normalizeDs(esDsAtNode[1], dXbin);
        // else logNormalizationFactors[nodeIdx] += normalizeDs(esDsAtNode[1], dXbin / hiLoRatio);

        // debugging
        logNormalizationFactors[nodeIdx] += normalizeDs(nodeIdx, lowRes, jtransforms);
        // System.out.println("    log-normalization factor total, for node " + nodeIdx + " = " + logNormalizationFactors[nodeIdx]);

    }

    @Override
    public double normalizeDs(int nodeIdx, boolean lowRes, boolean jtransforms) {
        double dx;
        double[] dsAtNode;
        if (lowRes) {
            dsAtNode = esDsLo[nodeIdx][1];
            dx = dXbin;
        } else {
            dsAtNode = esDsHi[nodeIdx][1];
            dx = dXbin / hiLoRatio;
        }

        double normalizationFactorFromDs;
        if (!jtransforms) {
            normalizationFactorFromDs = SSEUtils.calculateNormalizationFactor(dsAtNode, dx);

            // debugging
            // System.out.println("D's prior to normalization inside normalizeDs = " + Arrays.toString(dsAtNode));

            for (int i=0; i < dsAtNode.length; i += 2) {
                dsAtNode[i] /= normalizationFactorFromDs;
            }
        } else {
            normalizationFactorFromDs = SSEUtils.calculateNormalizationFactorJTransforms(dsAtNode, dx);

            // debugging
            // System.out.println("D's prior to normalization inside normalizeDs = " + Arrays.toString(dsAtNode));

            // just first half (real part)
            for (int i=0; i < (dsAtNode.length/2); i++) {
                dsAtNode[i] /= normalizationFactorFromDs;
                // there will be an additional step of normalization when in option 3
                // as a result of integrating the tip-end part in high-res
            }
        }

        // debugging
        // System.out.println("    normalization factor inside normalizeDs, for node " + nodeIdx + " = " + normalizationFactorFromDs);
        // System.out.println("    log-normalization factor inside normalizeDs, for node " + nodeIdx + " = " + Math.log(normalizationFactorFromDs));

        return Math.log(normalizationFactorFromDs); // the equivalent of ans[[1]] in diversitree (or ans.hi[[1]])
    }

    @Override
    public void mergeChildrenNodes(Node aNode, boolean jtransforms) {
        int nodeIdx = aNode.getNr();
        double nodeHeight = aNode.getHeight();

        List<Node> childrenNode = aNode.getChildren();
        Node leftChild = childrenNode.get(0);
        int leftChildIdx = leftChild.getNr();
        Node rightChild = childrenNode.get(1);
        int rightChildIdx = rightChild.getNr();

        double[][][] esDs;
        double[] birthRatesAtRightRes;
        int nUsefulXbinAtRightRes;
        if (nodeHeight < tc) {
            esDs = esDsHi;
            nUsefulXbinAtRightRes = nUsefulXbinsHi;
            birthRatesAtRightRes = birthRatesHi;
        } else {
            esDs = esDsLo;
            nUsefulXbinAtRightRes = nUsefulXbinsLo;
            birthRatesAtRightRes = birthRatesLo;
        }

        // debugging
        // System.out.println("\nMerging children");
        // System.out.println("E's before merging = " + Arrays.toString(esDs[leftChildIdx][0]));
        // System.out.println("D's left before merging = " + Arrays.toString(esDs[leftChildIdx][1]));
        // System.out.println("D's right before merging = " + Arrays.toString(esDs[rightChildIdx][1]));
        // System.out.println("Birth rates = " + Arrays.toString(birthRatesAtRightRes));

        // proper merging
        for (int i=0, j=0; i<(esDs[nodeIdx][0].length/2); ++i, j+=2) {
            if (i < nUsefulXbinAtRightRes) {
                if (jtransforms) {
                    esDs[nodeIdx][0][i] = esDs[leftChildIdx][0][i]; // E's
                    esDs[nodeIdx][1][i] = esDs[leftChildIdx][1][i] * esDs[rightChildIdx][1][i] * birthRatesAtRightRes[i]; // D's
                } else {
                    esDs[nodeIdx][0][j] = esDs[leftChildIdx][0][j]; // E's
                    esDs[nodeIdx][1][j] = esDs[leftChildIdx][1][j] * esDs[rightChildIdx][1][j] * birthRatesAtRightRes[i]; // D's
                }
            }
            // beyond useful continuous trait bins
            else {
                if (jtransforms) {
                    esDs[nodeIdx][0][i] = 0.0; // rest is 0's
                    esDs[nodeIdx][1][i] = 0.0;
                } else {
                    esDs[nodeIdx][0][j] = 0.0; // rest is 0's
                    esDs[nodeIdx][1][j] = 0.0;
                }
            }
        }
    }

    /*
     * Math-y methods start below
     */
    @Override
    public void populatefY(double aDt, boolean forceRecalcKernel, boolean dtChanged, boolean doFFT, boolean lowRes, boolean jtransforms) {
        super.populatefY(aDt, forceRecalcKernel, dtChanged, doFFT, lowRes, jtransforms);
    }

    @Override
    public void populatePriorProbAtRoot(double[] dsAtRoot, double dxAtRightRes, int nXBinsAtRightRes, int nUsefulXBinsAtRightRes, String aRootPriorType, boolean jtransforms) {
        // making sure things have the right dimensions (dividing dsAtRoot by two because it carries real and imaginary parts)
        if (providedPriorAtRoot && got2LowRes && ((priorProbsAtRoot.length / hiLoRatio) != nUsefulXbinsLo))
            throw new RuntimeException("ERROR: You need to specify " + nUsefulXbinsHi + " prior probability values, but you specified " + priorProbsAtRoot.length + ". Exiting...");
        else if (providedPriorAtRoot && !got2LowRes && (priorProbsAtRoot.length != nUsefulXbinsHi)) {
            throw new RuntimeException("ERROR: The number of prior probabilities did not match the number of discrete quantitative bins at high resolution. Exiting...");
        }

        // now we populate prior probs
        if (aRootPriorType.equals(FLAT)) {
            for (int i=0; i<nUsefulXBinsAtRightRes; ++i) {
                priorProbsAtRoot[i] = 1.0 / ((nXBinsAtRightRes - 1) * dxAtRightRes);
            }
        }

        else if (aRootPriorType.equals(OBS)) {
            double sumOfDsAtRoot = 0.0;
            for (int i=0, j=0; i<nXBinsAtRightRes; ++i, j+=2) {
                if (jtransforms) sumOfDsAtRoot += dsAtRoot[i];
                else sumOfDsAtRoot += dsAtRoot[j];
            }

            // debugging
            // System.out.println("priorProbsAtRoot.length = " + priorProbsAtRoot.length);
            // System.out.println("sumOfDsAtRoot = " + sumOfDsAtRoot);
            // System.out.println("dxAtRightRes = " + dxAtRightRes);

            for (int i=0, j=0; i<nUsefulXBinsAtRightRes; ++i, j+=2) {
                if (jtransforms) priorProbsAtRoot[i] = dsAtRoot[i] / (sumOfDsAtRoot * dxAtRightRes);
                else priorProbsAtRoot[i] = dsAtRoot[j] / (sumOfDsAtRoot * dxAtRightRes);
            }
        }

        else { throw new RuntimeException("ERROR: You specified an invalid prior probability distribution for the root. Exiting..."); }
    }

    public void populatePriorProbAtRootJTransforms(double[] dsAtRoot, double dxAtRightRes, int nXBinsAtRightRes, int nUsefulXBinsAtRightRes, String aRootPriorType) {
        // making sure things have the right dimensions (dividing dsAtRoot by two because it carries real and imaginary parts)
        if (providedPriorAtRoot && got2LowRes && ((priorProbsAtRoot.length / hiLoRatio) != nUsefulXbinsLo))
            throw new RuntimeException("ERROR: You need to specify " + nUsefulXbinsHi + " prior probability values, but you specified " + priorProbsAtRoot.length + ". Exiting...");
        else if (providedPriorAtRoot && !got2LowRes && (priorProbsAtRoot.length != nUsefulXbinsHi)) {
            throw new RuntimeException("ERROR: The number of prior probabilities did not match the number of discrete quantitative bins at high resolution. Exiting...");
        }

        // now we populate prior probs
        if (aRootPriorType.equals(FLAT)) {
            for (int i=0; i<nUsefulXBinsAtRightRes; ++i) {
                priorProbsAtRoot[i] = 1.0 / ((nXBinsAtRightRes - 1) * dxAtRightRes);
            }
        }
        else if (aRootPriorType.equals(OBS)) {
            double sumOfDsAtRoot = 0.0;
            for (int i=0; i<nXBinsAtRightRes; ++i) sumOfDsAtRoot += dsAtRoot[i];

            // debugging
            // System.out.println("priorProbsAtRoot.length = " + priorProbsAtRoot.length);
            // System.out.println("sumOfDsAtRoot = " + sumOfDsAtRoot);
            // System.out.println("dxAtRightRes = " + dxAtRightRes);

            for (int i=0; i<nUsefulXBinsAtRightRes; ++i) {
                priorProbsAtRoot[i] = dsAtRoot[i] / (sumOfDsAtRoot * dxAtRightRes);
            }
        }
        else { throw new RuntimeException("ERROR: You specified an invalid prior probability distribution for the root. Exiting..."); }
    }

    @Override
    public void doIntegrateInPlace(int nodeIdx, double aDt, boolean lowRes) {

        double[][] esDsAtNode, fftBufferEsDsAtNode;
        double[][] scratchAtNode;
        if (lowRes) {
            esDsAtNode = esDsLo[nodeIdx];
            fftBufferEsDsAtNode = fftBufferEsDsLo[nodeIdx];
            scratchAtNode = scratchLo[nodeIdx];
        } else {
            esDsAtNode = esDsHi[nodeIdx];
            fftBufferEsDsAtNode = fftBufferEsDsHi[nodeIdx];
            scratchAtNode = scratchHi[nodeIdx];
        }

        // debugging
        // System.out.println("esAtNode before propagate in t = " + Arrays.toString(esDsAtNode[0]));
        // System.out.println("dsAtNode before propagate in t = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode[0] before propagate in t = " + Arrays.toString(scratchAtNode[0]));
        // System.out.println("scratchAtNode[1] before propagate in t = " + Arrays.toString(scratchAtNode[1]));

        // integrate over birth and death events (low or high resolution inside)
        boolean jtransforms = false;
        propagateTInPlace(esDsAtNode, scratchAtNode, aDt, lowRes, jtransforms);

        // debugging
        // System.out.println("esAtNode after propagate in t and before x = " + Arrays.toString(esDsAtNode[0]));
        // System.out.println("dsAtNode after propagate in t and before x = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode[0] after propagate in t and before x = " + Arrays.toString(scratchAtNode[0]));
        // System.out.println("scratchAtNode[1] after propagate in t and before x = " + Arrays.toString(scratchAtNode[1]));

        // make normal kernel and FFTs it
        // populatefY(aDt, true, true, lowRes);

        // integrate over diffusion of quantitative trait
        // debugging
        // System.out.println("D's length before propagating in X = " + esDsAtNode[1].length);
        propagateXInPlace(esDsAtNode, fftBufferEsDsAtNode, scratchAtNode, lowRes);

        // debugging
        // System.out.println("esAtNode after propagate in t and x = " + Arrays.toString(esDsAtNode[0]));
        // System.out.println("dsAtNode after propagate in t and x = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratchAtNode[0] after propagate in t and x = " + Arrays.toString(scratchAtNode[0]));
        // System.out.println("scratchAtNode[1] after propagate in t and x = " + Arrays.toString(scratchAtNode[1]));
    }

    @Override
    public void propagateTInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double dt, boolean lowRes, boolean jtranforms) {
        if (jtranforms) {
            // grab scratch, dt and nDimensions from QuaSSEDistribution state
            // if (lowRes) SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesLo, deathRatesLo, dt, nUsefulXbinsLo, nDimensionsD);
            // else SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesHi, deathRatesHi, dt, nUsefulXbinsHi, nDimensionsD);
            if (lowRes)
                SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesLo, deathRatesLo, dt, nUsefulXbinsLo, nDimensionsD);
            else
                SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesHi, deathRatesHi, dt, nUsefulXbinsHi, nDimensionsD);
        } else {
            if (lowRes)
                SSEUtils.propagateEandDinTQuaSSEInPlaceSSTJavaFftService(esDsAtNode, scratchAtNode, birthRatesLo, deathRatesLo, dt, nUsefulXbinsLo, nDimensionsD);
            else
                SSEUtils.propagateEandDinTQuaSSEInPlaceSSTJavaFftService(esDsAtNode, scratchAtNode, birthRatesHi, deathRatesHi, dt, nUsefulXbinsHi, nDimensionsD);
        }
    }

    // JTransforms version
//    @Override
//    public void propagateXInPlace(int nodeIdx, double[][] esDsAtNode, double[][] scratchAtNode, boolean lowRes) {
//
//        // debugging
//        // System.out.println("esDsAtNode[1] = " + Arrays.toString(esDsAtNode[1]));
//        // System.out.println("scratch[1] = " + Arrays.toString(scratch[1]));
//
//        // grab dt and nDimensions from state
//        if (lowRes) SSEUtils.propagateEandDinXQuaLike(esDsAtNode, scratchAtNode, fYLo, nXbinsLo, nLeftNRightFlanksLo[0], nLeftNRightFlanksLo[1], nDimensionsE, nDimensionsD, fftForEandDLo);
//        else SSEUtils.propagateEandDinXQuaLike(esDsAtNode, scratchAtNode, fYHi, nXbinsHi, nLeftNRightFlanksHi[0], nLeftNRightFlanksHi[1], nDimensionsE, nDimensionsD, fftForEandDHi);
//    }

    @Override
    public void propagateXInPlace(double[][] esDsAtNode, double[][] fftBufferEsDsAtNode, double[][] scratchAtNode, boolean lowRes) {

        // debugging
        // System.out.println("esDsAtNode[1] = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratch[1] = " + Arrays.toString(scratch[1]));

        // grab dt and nDimensions from state
        if (lowRes) {
            SSEUtils.propagateEandDinXQuaLikeSSTJavaFftService(esDsAtNode, fftBufferEsDsAtNode, fftFYLo, scratchAtNode, nXbinsLo, nLeftNRightFlanksLo[0], nLeftNRightFlanksLo[1], nDimensionsE, nDimensionsD, jffts);
        }
        else {
            SSEUtils.propagateEandDinXQuaLikeSSTJavaFftService(esDsAtNode, fftBufferEsDsAtNode, fftFYHi, scratchAtNode, nXbinsHi, nLeftNRightFlanksHi[0], nLeftNRightFlanksHi[1], nDimensionsE, nDimensionsD, jffts);
        }
    }

    public void propagateXInPlaceJTransforms(double[][] esDsAtNode, double[][] scratchAtNode, boolean lowRes) {

        // debugging
        // System.out.println("esDsAtNode[1] = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratch[1] = " + Arrays.toString(scratch[1]));

        // grab dt and nDimensions from state
        if (lowRes) SSEUtils.propagateEandDinXQuaLike(esDsAtNode, scratchAtNode, fYLo, nXbinsLo, nLeftNRightFlanksLo[0], nLeftNRightFlanksLo[1], nDimensionsE, nDimensionsD, fftForEandDLo);
        else SSEUtils.propagateEandDinXQuaLike(esDsAtNode, scratchAtNode, fYHi, nXbinsHi, nLeftNRightFlanksHi[0], nLeftNRightFlanksHi[1], nDimensionsE, nDimensionsD, fftForEandDHi);
    }

    @Override
    public double getLogPFromRelevantObjects(double[][] esDsAtRoot, double sumOfLogNormalizationFactors, double[] birthRates, double dXAtRightRes, boolean jtransforms) {
        double[] esAtRoot = esDsAtRoot[0];
        double[] dsAtRoot = esDsAtRoot[1];

        // debugging
        // System.out.println("\nProcessing quantitites at root");
        // System.out.println("priorProbsAtRoot = " + Arrays.toString(priorProbsAtRoot));
        // System.out.println("esAtRoot = " + Arrays.toString(esAtRoot));
        // System.out.println("dsAtRoot = " + Arrays.toString(dsAtRoot));
        // System.out.println("esAtRoot.length = " + dsAtRoot.length);
        // System.out.println("priorProbsAtRoot.length = " + priorProbsAtRoot.length);
        // System.out.println("birthRates = " + Arrays.toString(birthRates));
        // System.out.println("birthRates.length = " + birthRates.length);

        // conditioning on survival of lineages
        double denomSumForConditioning = 0.0;

        /*
         * priorProbsAtRoot.length should be = nUsefulXBins at the right resolution
         * unlike the esDs array, we do not interdigitate values, which is why we
         * index priorProbsAtRoot with 'i', but index 'esAtRoot' and 'dsAtRoot'
         * below with 'j'
         */
        for (int i=0, j=0; i<priorProbsAtRoot.length; ++i, j+=2) {
            if (jtransforms) denomSumForConditioning += (priorProbsAtRoot[i] * birthRates[i] * Math.pow(1 - esAtRoot[i], 2));
            else denomSumForConditioning += (priorProbsAtRoot[i] * birthRates[i] * Math.pow(1 - esAtRoot[j], 2));
        }

        for (int i=0, j=0; i<priorProbsAtRoot.length; ++i, j+=2) {
            if (jtransforms) dsAtRoot[i] = dsAtRoot[i] / denomSumForConditioning * dXAtRightRes;
            else dsAtRoot[j] = dsAtRoot[j] / denomSumForConditioning * dXAtRightRes;
        }

        // debugging
        // System.out.println("Denominator sum = " + denomSumForConditioning);
        // System.out.println("dsAtRoot after conditioning on survival = " + Arrays.toString(dsAtRoot));

        // incorporating prior probabilities at root
        double prodSumPriorDs = 0.0;
        for (int i=0, j=0; i<priorProbsAtRoot.length; ++i, j+=2) {
            if (jtransforms) prodSumPriorDs += (priorProbsAtRoot[i] * dsAtRoot[i]);
            else prodSumPriorDs += (priorProbsAtRoot[i] * dsAtRoot[j]);
        }

        double thisLogLik = Math.log(prodSumPriorDs * dXAtRightRes) + sumOfLogNormalizationFactors; // denormalizing by adding sumOfLogNormalizationFactors

        return thisLogLik;
    }

    @Override
    public double calculateLogP() {

        // refreshing parameters
        populateMacroevolParams(false);

        // refreshing tree and root
        tree = treeInput.get();
        Node rootNode = tree.getRoot();
        int rootIdx = rootNode.getNr();
        
        // refreshing qu trait parameters
        drift = driftInput.get().getValue();
        diffusion = diffusionInput.get().getValue();
        
        boolean forceRecalcKernel = false;
        if (driftInput.get().somethingIsDirty() || diffusionInput.get().somethingIsDirty())
        	forceRecalcKernel = true;

        // start recursion for likelihood calculation
        boolean jtransforms = false;
        startRecursionAtRootNode(rootNode, forceRecalcKernel, jtransforms);

        // debugging
        // System.out.println("Log-normalization factors = " + Arrays.toString(logNormalizationFactors));

        // put it all together
        double sumOfLogNormalizationFactors = 0.0;
        for (double d: logNormalizationFactors) sumOfLogNormalizationFactors += d;
        double[][] esDsAtRootAtRightRes;
        double[] birthRatesAtRightRes;
        double dxAtRightRes;
        if (got2LowRes) {
            esDsAtRootAtRightRes = esDsLo[rootIdx]; // low res
            birthRatesAtRightRes = birthRatesLo;
            dxAtRightRes = dXbin;
        }
        else {
            esDsAtRootAtRightRes = esDsHi[rootIdx]; // high res
            birthRatesAtRightRes = birthRatesHi;
            dxAtRightRes = dXbin / hiLoRatio;
        }

        // debugging
        // System.out.println("logNormalizationFactors = " + Arrays.toString(logNormalizationFactors));

        // apply prior within
        double myLogP = getLogPFromRelevantObjects(esDsAtRootAtRightRes, sumOfLogNormalizationFactors, birthRatesAtRightRes, dxAtRightRes, jtransforms);

        return myLogP;
    }


    /*
     * Getters and setters
     */
    public double[] getLambda(boolean lowRes) {
        if (lowRes) return birthRatesLo;
        else return birthRatesHi;
    }

    public double[] getMu(boolean lowRes) {
        if (lowRes) return deathRatesLo;
        else return deathRatesHi;
    }

    public double[] getfY(boolean lowRes) {
        if (lowRes) return fYLo;
        else return fYHi;
    }

    public double[] getfftFY(boolean lowRes) {
        if (lowRes) return fftFYLo;
        else return fftFYHi;
    }

    public double[][][] getEsDs(boolean lowRes) {
        if (lowRes) return esDsLo;
        else return esDsHi;
    }

    public double[][][] getScratch(boolean lowRes) {
        if (lowRes) return scratchLo;
        else return scratchHi;
    }

    public double[][] getEsDsAtNode(int nodeIdx, boolean lowRes) {
        if (lowRes) return esDsLo[nodeIdx];
        return esDsHi[nodeIdx];
    }

    public double[][] getScratchAtNode(int nodeIdx, boolean lowRes) {
        if (lowRes) return scratchLo[nodeIdx];
        else return scratchHi[nodeIdx];
    }

    public int getNUsefulTraitBins(boolean lowRes) {
        if (lowRes) return nUsefulXbinsLo;
        else return nUsefulXbinsHi;
    }

    // for debugging
    public void setEsDsAtNodeElementAtDim(int nodeIdx, int dimIdx, int eleIdx, double val, boolean lowRes) {
        if (lowRes) esDsLo[nodeIdx][dimIdx][eleIdx] = val;
        else esDsHi[nodeIdx][dimIdx][eleIdx] = val;
    }

    // for debugging
    public void setGot2LowRes() {
        got2LowRes = true;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {

    }
}
