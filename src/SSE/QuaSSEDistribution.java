package SSE;

import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Node;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class QuaSSEDistribution extends QuaSSEProcess {

    public final Input<String> integrationBackendInput = new Input<>("integrationBackend",
            "Integration implementation: java or native T/X segments (requires fftBackend=fftw).", "java");
    private QuaSSENativeIntegrator nativeLo, nativeHi, nativeELo, nativeEHi;
    public final Input<RealParameter> fossilSamplingRateInput = new Input<>("fossilSamplingRate",
            "Constant non-removing fossil sampling rate; defaults to zero.");
    public final Input<RealParameter> presentSamplingProbabilityInput = new Input<>("presentSamplingProbability",
            "Probability of sampling a living lineage; defaults to one.");
    private double fossilSamplingRate, presentSamplingProbability = 1;
    public final Input<Boolean> conditionOnRootInput = new Input<>("conditionOnRoot",
            "Start at a sampled bifurcation (default true), rather than a single-lineage origin.", true);
    public final Input<RealParameter> originInput = new Input<>("origin", "Origin age, required only when conditionOnRoot=false.");
    public final Input<Boolean> conditionOnSamplingInput = new Input<>("conditionOnSampling",
            "Condition on sampled descendants of either kind; default true if neither sampling flag is supplied.");
    public final Input<Boolean> conditionOnRhoSamplingInput = new Input<>("conditionOnRhoSampling",
            "Condition on sampled living descendants; mutually exclusive with conditionOnSampling.");
    private boolean conditionOnRoot = true, conditionOnSampling = true, conditionOnRhoSampling;
    private double startingAge;
    private boolean hasFossils, hasTerminalFossils, hasBifurcations;
    // One reusable probability workspace per resolution: terminal-fossil E, then living-sample Q.
    private double[][] fossilELo, fossilEHi, scratchELo, scratchEHi, fftBufferELo, fftBufferEHi;
    public final Input<RealParameter> driftInput = new Input<>("drift", "Trait drift for an automatic grid.");
    public final Input<RealParameter> diffusionInput = new Input<>("diffusion", "Diffusion variance rate for an automatic grid.");
    private QuaSSEGrid automaticGrid;

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

        // Release any existing or newly allocated FFT resources if initialization fails.
        try {
            initializeGrid();
            super.initAndValidate(); // read dimensions and allocate FFT resources
            if (!refreshSamplingParameters())
                throw new IllegalArgumentException("QuaSSE sampling inputs must be scalar, with finite "
                        + "fossilSamplingRate >= 0 and presentSamplingProbability in [0,1].");
            resolveConditioning();
            validateSamplingTree();
            if (!refreshStartingAge())
                throw new IllegalArgumentException("Origin must be finite, scalar, and at least the tree root age.");
            String integration = integrationBackendInput.get();
            if (!"java".equals(integration) && !"native".equals(integration))
                throw new IllegalArgumentException("integrationBackend must be java or native.");
            if ("native".equals(integration)) {
                if (!"fftw".equals(fftBackendInput.get()))
                    throw new IllegalArgumentException("Native integration requires fftBackend=fftw.");
                nativeLo = new QuaSSENativeIntegrator(nXbinsLo, nDimensionsE + nDimensionsD);
                nativeHi = new QuaSSENativeIntegrator(nXbinsHi, nDimensionsE + nDimensionsD);
            }
            if (hasFossils || conditionOnRhoSampling) {
                fossilELo = new double[1][2 * nXbinsLo];
                fossilEHi = new double[1][2 * nXbinsHi];
                scratchELo = new double[1][2 * nXbinsLo];
                scratchEHi = new double[1][2 * nXbinsHi];
                fftBufferELo = new double[1][2 * nXbinsLo];
                fftBufferEHi = new double[1][2 * nXbinsHi];
                if (nativeLo != null) {
                    nativeELo = new QuaSSENativeIntegrator(nXbinsLo, 1);
                    nativeEHi = new QuaSSENativeIntegrator(nXbinsHi, 1);
                }
            }
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

            populateMacroevolParams();

            int nDimensionsFFT = nDimensions;
            initializeEsDs(nNodes, nDimensionsFFT, nXbinsLo, nXbinsHi);

            populateTipsEsDs(nDimensionsFFT, nXbinsHi, true, false);
        } catch (RuntimeException | Error failure) {
            close();
            throw failure;
        }
    }

    // Resolve configuration once; omission of both sampling flags preserves historical behavior.
    // Supplying only the rho flag must not implicitly enable the other condition.
    private void resolveConditioning() {
        conditionOnRoot = conditionOnRootInput.get();
        Boolean any = conditionOnSamplingInput.get(), living = conditionOnRhoSamplingInput.get();
        conditionOnSampling = any == null ? living == null : any;
        conditionOnRhoSampling = Boolean.TRUE.equals(living);
        if (conditionOnSampling && conditionOnRhoSampling)
            throw new IllegalArgumentException("Only one sampling condition may be true.");
        if (conditionOnRoot ? originInput.get() != null : originInput.get() == null)
            throw new IllegalArgumentException("Supply origin exactly when conditionOnRoot=false.");
        if (!dynamicallyAdjustDt && (!conditionOnRoot || conditionOnRhoSampling))
            throw new IllegalArgumentException("Origin and rho conditioning require dynDt=true.");
    }

    // An origin proposal changes the interval, not the initialized grid or its prior support.
    private boolean refreshStartingAge() {
        if (conditionOnRoot) {
            startingAge = tree.getRoot().getHeight();
            return true;
        }
        RealParameter origin = originInput.get();
        if (origin == null || origin.getDimension() != 1) return false;
        startingAge = origin.getValue();
        return Double.isFinite(startingAge) && startingAge >= tree.getRoot().getHeight();
    }

    // Root/rho conditioning needs a living observation on each daughter, not just anywhere.
    private boolean hasLivingDescendant(Node node) {
        if (node.isLeaf()) return node.getHeight() == 0 && !node.isDirectAncestor();
        for (Node child : node.getChildren()) if (hasLivingDescendant(child)) return true;
        return false;
    }

    // Read live parameters on every evaluation, including clean values restored after rejection.
    // Null inputs retain the legacy defaults and remain distinguishable for the dynDt requirement.
    private boolean refreshSamplingParameters() {
        RealParameter psi = fossilSamplingRateInput.get(), rho = presentSamplingProbabilityInput.get();
        if ((psi != null && psi.getDimension() != 1) || (rho != null && rho.getDimension() != 1))
            return false;
        fossilSamplingRate = psi == null ? 0 : psi.getValue();
        presentSamplingProbability = rho == null ? 1 : rho.getValue();
        return Double.isFinite(fossilSamplingRate) && fossilSamplingRate >= 0
                && Double.isFinite(presentSamplingProbability)
                && presentSamplingProbability >= 0 && presentSamplingProbability <= 1;
    }

    // BEAST represents an ancestor observation by an exactly zero-length leaf, not a unary node.
    // Validate that encoding before using binary-node helpers; repeat after tree proposals.
    private void validateSamplingTree() {
        hasFossils = hasTerminalFossils = hasBifurcations = false;
        boolean hasLivingSample = false;
        for (Node node : tree.getNodesAsArray()) {
            if (!Double.isFinite(node.getHeight()) || node.getHeight() < 0
                    || (!node.isRoot() && node.getHeight() > node.getParent().getHeight()))
                throw new IllegalArgumentException("Invalid age or branch length at node " + node.getID());
            if (node.isLeaf()) {
                if (node.getHeight() == 0 && !node.isDirectAncestor()) hasLivingSample = true;
                if (node.getHeight() > 0) {
                    hasFossils = true;
                    if (!node.isDirectAncestor()) hasTerminalFossils = true;
                }
                if (node.isDirectAncestor() && node.getHeight() == 0)
                    throw new IllegalArgumentException("Ancestor observations must have positive age: " + node.getID());
                if (q2dInput.get() instanceof NormalCenteredAtObservedLinkFn normal
                        && !normal.quTraitsInput.get().getKeysList().contains(node.getID()))
                    throw new IllegalArgumentException("Missing trait observation for sample " + node.getID());
            } else {
                if (!node.isFake()) hasBifurcations = true;
                if (node.getChildCount() != 2)
                    throw new IllegalArgumentException("QuaSSE requires binary nodes: " + node.getID());
                if (node.getLeft().isDirectAncestor() && node.getRight().isDirectAncestor())
                    throw new IllegalArgumentException("Ambiguous pair of ancestor leaves at node " + node.getID());
            }
        }
        if (tree.getRoot().isLeaf() || (conditionOnRoot && tree.getRoot().isFake()))
            throw new IllegalArgumentException("QuaSSE requires an ordinary bifurcating root, without a stem.");
        if (conditionOnRoot && conditionOnRhoSampling
                && (!hasLivingDescendant(tree.getRoot().getLeft())
                || !hasLivingDescendant(tree.getRoot().getRight())))
            throw new IllegalArgumentException("Root/rho conditioning needs living samples on both sides.");
        if (!hasLivingSample)
            throw new IllegalArgumentException("QuaSSE requires a living sample at age zero; fossil-only trees "
                    + "need a present-age offset and are not supported.");
        if (!dynamicallyAdjustDt && (hasFossils || fossilSamplingRateInput.get() != null
                || presentSamplingProbabilityInput.get() != null))
            throw new IllegalArgumentException("QuaSSE fossil events and explicit sampling inputs require dynDt=true.");
    }

    // Own only the generated grid; explicit grids retain caller-controlled initialization.
    // Attaching through Input registers it for ordinary BEAST dependency discovery before MCMC starts.
    private void initializeGrid() {
        QuaSSEGrid supplied = gridInput.get();
        if (supplied == null || supplied == automaticGrid) {
            if (driftInput.get() == null || diffusionInput.get() == null)
                throw new IllegalArgumentException("Without an explicit grid, QuaSSE requires drift and diffusion.");
            if (!(q2dInput.get() instanceof NormalCenteredAtObservedLinkFn normal))
                throw new IllegalArgumentException("Automatic grids require NormalCenteredAtObservedLinkFn; "
                        + "supply an explicit grid for other links.");
            if (automaticGrid == null) automaticGrid = new QuaSSEGrid();
            automaticGrid.initByName("tree", treeInput.get(), "traits", normal.quTraitsInput.get(),
                    "drift", driftInput.get(), "diffusion", diffusionInput.get());
            gridInput.setValue(automaticGrid, this);
        } else {
            // Serialization may expose the generated grid and the same parameter references twice.
            if ((driftInput.get() != null && driftInput.get() != supplied.driftInput.get())
                    || (diffusionInput.get() != null && diffusionInput.get() != supplied.diffusionInput.get()))
                throw new IllegalArgumentException(
                        "Likelihood drift/diffusion must reference the explicit grid's parameters.");
        }
    }

    // Also called by parent initialization: release old resolution-specific resources before resizing.
    @Override
    public void close() {
        if (nativeLo != null) nativeLo.close();
        if (nativeHi != null) nativeHi.close();
        if (nativeELo != null) nativeELo.close();
        if (nativeEHi != null) nativeEHi.close();
        nativeLo = nativeHi = nativeELo = nativeEHi = null;
        fossilELo = fossilEHi = scratchELo = scratchEHi = fftBufferELo = fftBufferEHi = null;
        super.close();
    }

    @Override
    public void populateMacroevolParams() {
        // Refresh once, then populate both rulers. Another likelihood may share the link and have
        // already refreshed its scalars; that must not prevent refreshing this likelihood's arrays.
        q2mLambda.refreshParams();
        q2mMu.refreshParams();
        birthRatesLo = q2mLambda.getY(xLo, birthRatesLo, true);
        birthRatesHi = q2mLambda.getY(xHi, birthRatesHi, true);
        deathRatesLo = q2mMu.getY(xLo, deathRatesLo, true);
        deathRatesHi = q2mMu.getY(xHi, deathRatesHi, true);
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
                    // Living boundary E(x,0) = 1 − ρ; fossil boundaries are set during evaluation.
                    Arrays.fill(esDsHi[nodeIdx][i], 0);
                    Arrays.fill(esDsLo[nodeIdx][i], 0);
                    if (tip.getHeight() == 0) {
                        int stride = jtransforms ? 1 : 2;
                        for (int j = 0; j < nUsefulXbinsHi; j++)
                            esDsHi[nodeIdx][i][stride * j] = 1 - presentSamplingProbability;
                        for (int j = 0; j < nUsefulXbinsLo; j++)
                            esDsLo[nodeIdx][i][stride * j] = 1 - presentSamplingProbability;
                    }
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

            // Only the grid at the observation's age is used; an unused coarse density may underflow.
            boolean lowRes = tip.getHeight() >= tc;
            double dsSum = 0;
            for (double d : (lowRes ? esDsLo : esDsHi)[nodeIdx][1]) dsSum += d;
            if (!(dsSum > 0) || !Double.isFinite(dsSum))
                throw new ArithmeticException("Invalid observation density for " + tipName
                        + " on the " + (lowRes ? "coarse" : "fine")
                        + " grid; check observation SD and trait spacing.");
        }
    }

    // Initialize actual sample events after support checks, keeping scalar factors in log space.
    // Virtual leaves store only g; their ψ factor and continuing D are handled at the parent event.
    private void initializeSampleEvents() {
        for (Node tip : tree.getExternalNodes()) {
            if (tip.isDirectAncestor()) continue;
            int nodeIdx = tip.getNr();
            if (tip.getHeight() == 0) {
                logNormalizationFactors[nodeIdx] = Math.log(presentSamplingProbability);
                continue;
            }
            // E is independent of the observed subtree. Restart at the present for every fossil,
            // using the same steps and transport as pruning, but without likelihood normalization.
            boolean lowRes = tip.getHeight() >= tc;
            double[] e = probabilityAtAge(tip.getHeight(), fossilSamplingRate);
            double[][] rows = (lowRes ? esDsLo : esDsHi)[nodeIdx];
            System.arraycopy(e, 0, rows[0], 0, e.length);
            // A terminal non-removing fossil contributes ψ g E: its continuation is unobserved.
            for (int i = 0; i < (lowRes ? nUsefulXbinsLo : nUsefulXbinsHi); ++i)
                rows[1][2 * i] *= e[2 * i];
            logNormalizationFactors[nodeIdx] = Math.log(fossilSamplingRate);
        }
    }

    // E and Q share the same equation: Q uses zero sampling rate because fossils do not remove
    // lineages. Restart from 1-rho; transfer at tc before evaluating an endpoint exactly at tc.
    private double[] probabilityAtAge(double age, double samplingRate) {
        Arrays.fill(fossilELo[0], 0);
        Arrays.fill(fossilEHi[0], 0);
        boolean lowRes = age >= tc;
        double[][] initial = tc == 0 ? fossilELo : fossilEHi;
        int useful = tc == 0 ? nUsefulXbinsLo : nUsefulXbinsHi;
        for (int i = 0; i < useful; ++i) initial[0][2 * i] = 1 - presentSamplingProbability;
        if (tc > 0) {
            integrateSegment(fossilEHi, scratchEHi, fftBufferEHi, Math.min(age, tc),
                    true, dtMax, false, false, nativeEHi, samplingRate);
            if (lowRes)
                SSEUtils.hiToLoTransferInPlace(fossilEHi[0], fossilELo[0], hiLoIdxs4Transfer, false);
        }
        if (lowRes)
            integrateSegment(fossilELo, scratchELo, fftBufferELo, age - tc,
                    true, dtMax, true, false, nativeELo, samplingRate);
        return (lowRes ? fossilELo : fossilEHi)[0];
    }

    @Override
    public void startRecursionAtRootNode(Node rootNode, boolean forceRecalcKernel, boolean jtransforms) {
        int rootIdx = rootNode.getNr();
        double rootHeight = conditionOnRoot ? rootNode.getHeight() : startingAge;

        // start recursion
        processInternalNode(rootNode, forceRecalcKernel, jtransforms);
        // A zero-length stem is exactly the identity; it adds no event or normalization.
        if (!conditionOnRoot && startingAge > rootNode.getHeight())
            processInterval(rootIdx, rootNode.getHeight(), startingAge, forceRecalcKernel, jtransforms);

        // we're done pruning, let's deal with the prior probs at root now
        double[][] esDsAtRoot;
        double dxAtRightRes;
        int nXBinsRightRes, nUsefulXBinsRightRes;
        if (rootHeight >= tc) {
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

        if (priorProbsAtRoot == null || priorProbsAtRoot.length != nUsefulXBinsRightRes)
            priorProbsAtRoot = new double[nUsefulXBinsRightRes];

        if (conditionOnRoot)
            populatePriorProbAtRoot(esDsAtRoot[1], dxAtRightRes, nXBinsRightRes, nUsefulXBinsRightRes,
                    rootPriorType, jtransforms);
        else {
            // Origin weights are densities on useful bins. Common scaling of D cancels in
            // Observed weights; these are comparison weights, not a data-independent prior.
            double sum = 0;
            for (int i = 0; i < nUsefulXBinsRightRes; ++i) {
                priorProbsAtRoot[i] = rootPriorType.equals(FLAT) ? 1 : esDsAtRoot[1][jtransforms ? i : 2 * i];
                sum += priorProbsAtRoot[i];
            }
            if (!(sum > 0) || !Double.isFinite(sum))
                throw new ArithmeticException("Invalid origin trait weights.");
            for (int i = 0; i < nUsefulXBinsRightRes; ++i) priorProbsAtRoot[i] /= sum * dxAtRightRes;
        }
    }

    @Override
    public void processInternalNode(Node aNode, boolean forceRecalcKernel, boolean jtransforms) {

        // debugging
        // System.out.println("\nDoing node " + aNode.getID());

        // recur if internal node or sampled ancestor
        if (!aNode.isLeaf()) {
            if (aNode.isFake()) {
                Node continuing = aNode.getNonDirectAncestorChild();
                processInternalNode(continuing, forceRecalcKernel, jtransforms);
                double[][][] rows = aNode.getHeight() >= tc ? esDsLo : esDsHi;
                int useful = aNode.getHeight() >= tc ? nUsefulXbinsLo : nUsefulXbinsHi;
                double[][] parent = rows[aNode.getNr()], child = rows[continuing.getNr()];
                double[] observation = rows[aNode.getDirectAncestorChild().getNr()][1];
                // An ancestor contributes ψ g Dcontinuing, with no λ or extra E factor.
                // Copy E unchanged; the virtual leaf has no branch calculation or likelihood scale.
                for (int i = 0; i < useful; ++i) {
                    int j = jtransforms ? i : 2 * i;
                    parent[0][j] = child[0][j];
                    parent[1][j] = child[1][j] * observation[j];
                }
                logNormalizationFactors[aNode.getNr()] += Math.log(fossilSamplingRate);
            } else {
                for (Node childNode: aNode.getChildren()) {
                    processInternalNode(childNode, forceRecalcKernel, jtransforms); // recur
                }

                // after recursion, prepare initial conditions for integrating this node
                mergeChildrenNodes(aNode, jtransforms);
            }
        }

        // unless we're at the root, now we integrate this branch
        if (!aNode.isRoot()) processBranch(aNode, forceRecalcKernel, jtransforms);
    }

    // main branch function
    @Override
    public void processBranch(Node aNode, boolean forceRecalcKernel, boolean jtransforms) {
        processInterval(aNode.getNr(), aNode.getHeight(), aNode.getParent().getHeight(),
                forceRecalcKernel, jtransforms);
    }

    // Propagate one node row backwards over an age interval, retaining the existing resolution
    // transition and accumulating each normalization in that row's scale exactly once.
    private void processInterval(int nodeIdx, double startTime, double endTime,
            boolean forceRecalcKernel, boolean jtransforms) {
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
        if (startTime >= tc) {
            esDsAtNode = esDsLo[nodeIdx];
            scratchAtNode = scratchLo[nodeIdx];

            // debugging
            // System.out.println("OPTION 1: E's size should be 2 * low-res = " + esDsAtNode[0].length);
            // System.out.println("Unnormalized D's at low-res = " + Arrays.toString(esDsAtNode[1]));

            // normalize and record normalization factor before integration
            // logNormalizationFactors[nodeIdx] = normalizeDs(esDsAtNode[1], dXbin); // normalize D's and returns factor, which we record
            logNormalizationFactors[nodeIdx] += normalizeDs(nodeIdx, true, jtransforms); // normalize D's and returns factor, which we record

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
            logNormalizationFactors[nodeIdx] += normalizeDs(nodeIdx, false, jtransforms); // normalize D's and returns factor, which we record

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
            logNormalizationFactors[nodeIdx] += normalizeDs(nodeIdx, false, jtransforms); // normalize D's and returns factor, which we record

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
    public void integrateLength(int nodeIdx, double[][] rows, double[][] scratch, double length,
            boolean dynamicallyAdjust, double maxDt, boolean lowRes, boolean forceRecalcKernel,
            boolean jtransforms) {
        // An empty interval is the identity, including branch ends exactly at tc.
        if (length == 0) return;
        integrateSegment(rows, scratch, (lowRes ? fftBufferEsDsLo : fftBufferEsDsHi)[nodeIdx],
                length, dynamicallyAdjust, maxDt, lowRes, forceRecalcKernel, lowRes ? nativeLo : nativeHi, fossilSamplingRate);
        logNormalizationFactors[nodeIdx] += normalizeDs(nodeIdx, lowRes, jtransforms);
    }

    // Advance explicit E/D or E-only storage; the caller alone owns likelihood scaling.
    // This shares timestep selection, kernels and Strang ordering without inventing a D row for E.
    private void integrateSegment(double[][] rows, double[][] scratch, double[][] transformed,
            double length, boolean dynamicallyAdjust, double maxDt, boolean lowRes, boolean force,
            QuaSSENativeIntegrator owner, double samplingRate) {
        if (length == 0) return;
        double intervals = dynamicallyAdjust ? Math.ceil(length / maxDt) : Math.floor(length / maxDt);
        double dt = dynamicallyAdjust ? length / intervals : maxDt;
        // Compatibility: legacy fixed-dt mode drops the remainder. Sampling inputs require dynDt=true;
        // repairing the old mode's timestep convention is a separate numerical change.
        populatefY(dt, force, true, lowRes, false);
        if (owner != null && Double.isFinite(intervals) && intervals >= 0
                && intervals <= Integer.MAX_VALUE && intervals == Math.rint(intervals)) {
            integrateNativeSegment(rows, owner, dt, (int) intervals, lowRes, samplingRate);
        } else {
            // Compatibility: preserve the existing int-loop limitation outside JNI's step-count range.
            for (int i = 0; i < intervals; ++i)
                integrateStep(rows, scratch, transformed, dt, lowRes, owner, samplingRate);
        }
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

        double normalizationFactorFromDs = jtransforms
                ? SSEUtils.calculateNormalizationFactorJTransforms(dsAtNode, dx)
                : SSEUtils.calculateNormalizationFactor(dsAtNode, dx);
        if (!(normalizationFactorFromDs > 0) || !Double.isFinite(normalizationFactorFromDs))
            throw new ArithmeticException("Invalid D normalization at node " + nodeIdx
                    + " on the " + (lowRes ? "coarse" : "fine") + " grid.");
        if (!jtransforms) {

            // debugging
            // System.out.println("D's prior to normalization inside normalizeDs = " + Arrays.toString(dsAtNode));

            for (int i=0; i < dsAtNode.length; i += 2) {
                dsAtNode[i] /= normalizationFactorFromDs;
            }
        } else {

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
    public void populatefY(double aDt, boolean forceRecalcKernel, boolean doFFT, boolean lowRes, boolean jtransforms) {
        super.populatefY(aDt, forceRecalcKernel, doFFT, lowRes, jtransforms);
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

    // Copy current inputs into the supplied native owner; owners cache storage, not model state.
    private void integrateNativeSegment(double[][] rows, QuaSSENativeIntegrator owner,
            double dt, int steps, boolean lowRes, double samplingRate) {
        int[] flanks = lowRes ? nLeftNRightFlanksLo : nLeftNRightFlanksHi;
        owner.integrateSegment(rows, lowRes ? birthRatesLo : birthRatesHi,
                lowRes ? deathRatesLo : deathRatesHi, lowRes ? fftFYLo : fftFYHi,
                samplingRate, dt, steps, flanks[0], flanks[1]);
    }

    // Single-step entry point retained for callers using node storage and a prepared full-dt kernel.
    @Override
    public void doIntegrateInPlace(int nodeIdx, double dt, boolean lowRes) {
        integrateStep((lowRes ? esDsLo : esDsHi)[nodeIdx], (lowRes ? scratchLo : scratchHi)[nodeIdx],
                (lowRes ? fftBufferEsDsLo : fftBufferEsDsHi)[nodeIdx], dt, lowRes, lowRes ? nativeLo : nativeHi, fossilSamplingRate);
    }

    // Strang splitting: reaction half-steps surround transport with the prepared full-dt kernel.
    private void integrateStep(double[][] rows, double[][] scratch, double[][] transformed,
            double dt, boolean lowRes, QuaSSENativeIntegrator owner, double samplingRate) {
        if (owner != null) {
            integrateNativeSegment(rows, owner, dt, 1, lowRes, samplingRate);
            return;
        }
        propagateReaction(rows, scratch, dt / 2, lowRes, false, samplingRate);
        propagateXInPlace(rows, transformed, scratch, lowRes);
        propagateReaction(rows, scratch, dt / 2, lowRes, false, samplingRate);
    }

    @Override
    public void propagateTInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double dt, boolean lowRes, boolean jtranforms) {
        propagateReaction(esDsAtNode, scratchAtNode, dt, lowRes, jtranforms, fossilSamplingRate);
    }

    // Pass the rate explicitly so Q propagation cannot mutate the data-likelihood model.
    private void propagateReaction(double[][] esDsAtNode, double[][] scratchAtNode, double dt,
            boolean lowRes, boolean jtranforms, double samplingRate) {
        if (jtranforms) {
            // grab scratch, dt and nDimensions from QuaSSEDistribution state
            // if (lowRes) SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesLo, deathRatesLo, samplingRate, dt, nUsefulXbinsLo, esDsAtNode.length - 1);
            // else SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesHi, deathRatesHi, samplingRate, dt, nUsefulXbinsHi, esDsAtNode.length - 1);
            if (lowRes)
                SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesLo, deathRatesLo, samplingRate, dt, nUsefulXbinsLo, esDsAtNode.length - 1);
            else
                SSEUtils.propagateEandDinTQuaSSEInPlace(esDsAtNode, scratchAtNode, birthRatesHi, deathRatesHi, samplingRate, dt, nUsefulXbinsHi, esDsAtNode.length - 1);
        } else {
            if (lowRes)
                SSEUtils.propagateEandDinTQuaSSEInPlaceSSTJavaFftService(esDsAtNode, scratchAtNode, birthRatesLo, deathRatesLo, samplingRate, dt, nUsefulXbinsLo, esDsAtNode.length - 1);
            else
                SSEUtils.propagateEandDinTQuaSSEInPlaceSSTJavaFftService(esDsAtNode, scratchAtNode, birthRatesHi, deathRatesHi, samplingRate, dt, nUsefulXbinsHi, esDsAtNode.length - 1);
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
//        if (lowRes) SSEUtils.propagateEandDinXQuaLike(esDsAtNode, scratchAtNode, fYLo, nXbinsLo, nLeftNRightFlanksLo[0], nLeftNRightFlanksLo[1], 1, esDsAtNode.length - 1, fftForEandDLo);
//        else SSEUtils.propagateEandDinXQuaLike(esDsAtNode, scratchAtNode, fYHi, nXbinsHi, nLeftNRightFlanksHi[0], nLeftNRightFlanksHi[1], 1, esDsAtNode.length - 1, fftForEandDHi);
//    }

    @Override
    public void propagateXInPlace(double[][] esDsAtNode, double[][] fftBufferEsDsAtNode, double[][] scratchAtNode, boolean lowRes) {

        QuaSSENativeIntegrator nativeX = esDsAtNode.length == nDimensions ? (lowRes ? nativeLo : nativeHi) : null;
        if (nativeX != null) {
            int[] flanks = lowRes ? nLeftNRightFlanksLo : nLeftNRightFlanksHi;
            nativeX.propagateX(esDsAtNode, lowRes ? fftFYLo : fftFYHi, flanks[0], flanks[1]);
            return;
        }

        // debugging
        // System.out.println("esDsAtNode[1] = " + Arrays.toString(esDsAtNode[1]));
        // System.out.println("scratch[1] = " + Arrays.toString(scratch[1]));

        // Use the same resolution's prepared transform for E/D as for its Gaussian kernel.
        if (lowRes) {
            SSEUtils.propagateEandDinXQuaSSE(esDsAtNode, fftBufferEsDsAtNode, fftFYLo, scratchAtNode,
                    nXbinsLo, nLeftNRightFlanksLo[0], nLeftNRightFlanksLo[1], 1, esDsAtNode.length - 1, fftLo);
        }
        else {
            SSEUtils.propagateEandDinXQuaSSE(esDsAtNode, fftBufferEsDsAtNode, fftFYHi, scratchAtNode,
                    nXbinsHi, nLeftNRightFlanksHi[0], nLeftNRightFlanksHi[1], 1, esDsAtNode.length - 1, fftHi);
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
        if (!conditionOnRoot || !conditionOnSampling)
            return conditionedLogP(esDsAtRoot, sumOfLogNormalizationFactors, birthRates,
                    dXAtRightRes, jtransforms);
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

        // Root-prior values are densities: the conditioning integral is the sum times dx.
        // Integrate before dividing D; multiplying D by dx afterwards adds an erroneous dx^2.
        denomSumForConditioning *= dXAtRightRes;
        if (!(denomSumForConditioning > 0) || !Double.isFinite(denomSumForConditioning))
            throw new ArithmeticException("Invalid root-conditioning integral on the "
                    + (got2LowRes ? "coarse" : "fine") + " grid.");
        for (int i=0, j=0; i<priorProbsAtRoot.length; ++i, j+=2) {
            if (jtransforms) dsAtRoot[i] /= denomSumForConditioning;
            else dsAtRoot[j] /= denomSumForConditioning;
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

        // Branch normalization cannot catch disjoint underflowed root partials: validate the
        // final product integral too, rather than reporting numerical loss as zero model support.
        if (!(prodSumPriorDs > 0) || !Double.isFinite(prodSumPriorDs))
            throw new ArithmeticException("Invalid root likelihood integral on the "
                    + (got2LowRes ? "coarse" : "fine") + " grid.");

        double thisLogLik = Math.log(prodSumPriorDs * dXAtRightRes) + sumOfLogNormalizationFactors; // denormalizing by adding sumOfLogNormalizationFactors

        return thisLogLik;
    }

    // Integrate first, then condition: averaging pointwise ratios would change the trait model.
    // Root mode retains the starting birth-density factor; origin mode has no event at its start.
    private double conditionedLogP(double[][] rows, double scales, double[] birthRates,
            double dx, boolean jtransforms) {
        double[] failure = conditionOnRhoSampling ? probabilityAtAge(startingAge, 0) : rows[0];
        double numerator = 0, denominator = 0;
        for (int i = 0; i < priorProbsAtRoot.length; ++i) {
            int j = jtransforms ? i : 2 * i;
            // Probability workspaces use interleaved storage even for legacy packed-row callers.
            int probabilityIndex = conditionOnRhoSampling ? 2 * i : j;
            double success = conditionOnSampling || conditionOnRhoSampling ? 1 - failure[probabilityIndex] : 1;
            numerator += priorProbsAtRoot[i] * rows[1][j];
            denominator += priorProbsAtRoot[i] * (conditionOnRoot ? birthRates[i] * success * success : success);
        }
        numerator *= dx;
        denominator = !conditionOnRoot && !conditionOnSampling && !conditionOnRhoSampling ? 1 : denominator * dx;
        if (!(numerator > 0) || !Double.isFinite(numerator)
                || !(denominator > 0) || !Double.isFinite(denominator))
            throw new ArithmeticException("Invalid starting likelihood or conditioning integral.");
        return Math.log(numerator) - Math.log(denominator) + scales;
    }

    @Override
    public double calculateLogP() {

        tree = treeInput.get();
        if (!refreshSamplingParameters() || !refreshStartingAge() || !refreshGrid())
            return logP = Double.NEGATIVE_INFINITY;
        try {
            validateSamplingTree();
        } catch (IllegalArgumentException invalidTree) {
            return logP = Double.NEGATIVE_INFINITY;
        }
        if ((hasFossils && fossilSamplingRate == 0) || presentSamplingProbability == 0)
            return logP = Double.NEGATIVE_INFINITY;

        // refreshing parameters
        populateMacroevolParams();

        // Identifiable mathematical zeros precede any observation/product normalization.
        // An origin tree with only ancestral samples needs no birth; fake nodes are observations.
        if ((hasBifurcations && Arrays.stream(birthRatesLo).allMatch(x -> x == 0)
                && Arrays.stream(birthRatesHi).allMatch(x -> x == 0))
                || (hasTerminalFossils && presentSamplingProbability == 1
                && Arrays.stream(deathRatesLo).allMatch(x -> x == 0)
                && Arrays.stream(deathRatesHi).allMatch(x -> x == 0)))
            return logP = Double.NEGATIVE_INFINITY;

        // refreshing tree and root
        tree = treeInput.get();
        Node rootNode = tree.getRoot();
        int rootIdx = rootNode.getNr();
        
        // Pruning mutates every node partial in place, so a full calculation must rebuild them from tips.
        Arrays.fill(logNormalizationFactors, 0.0);
        for (int nodeIdx = 0; nodeIdx < esDsLo.length; nodeIdx++) {
            for (int dimension = 0; dimension < nDimensions; dimension++) {
                Arrays.fill(esDsLo[nodeIdx][dimension], 0.0);
                Arrays.fill(esDsHi[nodeIdx][dimension], 0.0);
            }
        }
        // Refresh q2d's cached inputs; the forced call then repopulates arrays that were just cleared.
        q2d.refreshParams();
        populateTipsEsDs(nDimensions, nXbinsHi, true, false);
        initializeSampleEvents();
        
        boolean forceRecalcKernel = false;
        // Kernel keys track actual parameter values, including clean restored values after rejection.

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
        // Distribution owns the scalar snapshot used by CompoundDistribution and BEAST rejection.
        logP = getLogPFromRelevantObjects(esDsAtRootAtRightRes, sumOfLogNormalizationFactors, birthRatesAtRightRes, dxAtRightRes, jtransforms);

        return logP;
    }

    // Reconcile the live grid with current parameters, including clean values restored after rejection.
    // Validate before touching any live arrays; invalid proposals leave the last usable grid intact.
    private boolean refreshGrid() {
        if (dynamicDtInput.get().getValue() != dynamicallyAdjustDt)
            throw new IllegalArgumentException("QuaSSE dynDt must remain fixed after initialization.");
        double[] previousX = xLo;
        if (!refreshGridGeometry()) return false;
        if (previousX != xLo) {
            birthRatesLo = new double[nUsefulXbinsLo];
            birthRatesHi = new double[nUsefulXbinsHi];
            deathRatesLo = new double[nUsefulXbinsLo];
            deathRatesHi = new double[nUsefulXbinsHi];
            priorProbsAtRoot = null;
        }
        return true;
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
