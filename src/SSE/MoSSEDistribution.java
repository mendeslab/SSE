package SSE;

import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import org.shared.array.RealArray;

import java.util.List;
import java.util.Random;

public class MoSSEDistribution extends QuaSSEDistribution {

    final public Input<Alignment> alignInput = new Input<>("seqAlign", "Molecular alignment.", Input.Validate.REQUIRED);
    final public Input<SiteModelInterface> siteModelInput = new Input<>("siteModel", "Molecular site model.", Input.Validate.REQUIRED);

    private int nDimensions = 5;
    private Alignment algn;
    private SiteModel.Base siteModel;
    private SubstitutionModel substModel;
    private MoSSELikelihoodCore likelihoodCore;

    private int nSitePat;
    private double[] instRateMat; // flattened out matrix

    // state that matters for calculateLogP
    private double[][][] esDsLo, esDsHi; // first dimension are all nodes, second dimension are E's and D's, third are the bins
    private double[][] scratch;


    @Override
    public void initAndValidate() {

        // super will (1) populate macroevol parameters, ...
        super.initAndValidate();

        String what2check = "algn";
        algn = alignInput.get();
        areInputsOK(what2check); // checks if tree and MSA are OK

        what2check = "sitemodel";
        areInputsOK(what2check);
        siteModel = (SiteModel.Base) siteModelInput.get();

        nSitePat = algn.getPatternCount();

        int dtRes = 1;
        instRateMat = new double[4 * 4];

        int nDimensionsFFT = 5; // E, and 4 D's (one for each nuc)
        initializeEsDs(tree.getNodeCount(), nSitePat, nDimensionsFFT, nXbinsLo, nXbinsHi);
    }

    // needs extra dimension of nSitePat
    public void initializeEsDs(int nNodes, int nSitePat, int nDimensionsFFT, int nXbinsLo, int nXbinsHi) {
        esDsLo = new double[nNodes][nDimensionsFFT][nXbinsHi];
        esDsHi = new double[nNodes][nDimensionsFFT][nXbinsHi];
        // do stuff with nDimensionsFFT
    }

    @Override
    public void populateTipsEsDs(int nDimensionsFFT, int nUsefulXbinsHi, boolean ignoreRefresh, boolean jtransforms) {
        ;
    }

    @Override
    public void populatePriorProbAtRoot(double[] dsAtRoot, double dxAtRightRes, int nXBinsAtRightRes, int nUsefulXBinsAtRightRes, String aRootPriorType, boolean jtransforms) {
        ;
    }

    @Override
    public void processBranch(Node aNode, boolean forceRecalcKernel, boolean jtransforms) {

        // TODO stuff
        doMoSSEIntegrateInPlace(aNode.getNr(), esDsHi[0], fftBufferEsDsHi[0], scratch, 0.0, true, dtMax, true);
    }

    @Override
    public void processInternalNode(Node aNode, boolean forceRecalcKernel, boolean jtransforms) {
    }

    @Override
    public void startRecursionAtRootNode(Node rootNode, boolean forceRecalcKernel, boolean jtransforms) {

    }

    /*
     * Because we're dealing with transition probability matrix implementations,
     * we need start times, and need to deal with the first dt in a different way.
     * MoSSE will thus have its special 'doIntegrateInPlace'.
     */
    private void doMoSSEIntegrateInPlace(int nodeIdx, double[][] esDsAtNode, double[][] fftBufferEsDsAtNode, double[][] scratchAtNode, double startTime, boolean isFirstDt, double dt, boolean lowRes) {

        // integrate over birth and death events (either at low or high resolution inside)
        boolean jtransforms = false;
        this.propagateTInPlace(esDsAtNode, scratchAtNode, dt, lowRes, jtransforms);

        /*
         * For the step below, I am working on the MoSSELikelihoodCore
         * class.
         *
         * This class will never deal with partials at internal nodes as
         * LikelihoodCore does, because internal nodes are dealt with
         * by the MoSSE class, which just does the D * D product.
         *
         * Instead, this class will update a single partial (per bin) along
         * a branch, at every dt. P(dt,r) will be constant, and the
         * partial update function should take the partial as input and
         * just update it in place.
         *
         * Then I'd need two separate functions for updating partials:
         * (1) For the very first dt of a leaf (where missing data is dealt with)
         * (2) For all other dts
         */

        // MoSSE: integrate over substitution events
        propagateSubst(esDsAtNode, startTime, dtMax, isFirstDt, lowRes);

        // integrate over diffusion of substitution rate
        this.propagateXInPlace(esDsAtNode, fftBufferEsDsAtNode, scratchAtNode, lowRes);

        // return null;
    }

    @Override
    public void propagateTInPlace(double[][] esDsAtNode, double[][] scratchAtNode, double dt, boolean lowRes, boolean jtransforms) {
        super.propagateTInPlace(esDsAtNode, scratchAtNode, dt, lowRes, false);

        // return null;
    }

    @Override
    public void propagateXInPlace(double[][] esDsAtNode, double[][] fftBufferEsDs, double[][] scratchAtNode, boolean lowRes) {
        super.propagateXInPlace(esDsAtNode, fftBufferEsDs, scratchAtNode, lowRes);
    }

    public void propagateSubst(double[][] esDsAtNode, double startTime, double aDt, boolean isFirstDt, boolean lowRes) {

        double endTime = startTime + aDt;
        double[] xRuler;
        if (lowRes) xRuler = xLo;
        else xRuler = xHi;

            for (int i=0; i<nXbinsHi; i++) {
                double xValue = xRuler[i];
                double[] esDsAtNodeAtBin = esDsAtNode[i];
                substModel.getTransitionProbabilities(null, startTime, endTime, xValue, instRateMat); // P(dt,r), where r = x = subst rate
                likelihoodCore.updateMatrices(instRateMat); // set Q matrix inside likelihoodCore

                if (isFirstDt) likelihoodCore.updatePartialInPlaceFirstDt(esDsAtNodeAtBin, aDt); // first dt deals with missing data
                else likelihoodCore.updatePartialInPlace(esDsAtNodeAtBin, aDt);
            }

    }

    @Override
    public double calculateLogP() {

        populateMacroevolParams();

        // iterating over all unique site patterns (multiply by their frequencies)
//        for (int i=0; i<nSitePat; ++i) {
//            // process(root);
//        }

        return 0.0;
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

    private void areInputsOK(String what2Check) {

        // alignment should have same #taxa as tree, and only works with nucleotide alignments
        if (what2Check == "algn") {
            if (alignInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) throw new IllegalArgumentException("The number of nodes in the tree does not match the number of sequences. Exiting...");
            if (alignInput.get().getMaxStateCount() != 4) throw new IllegalArgumentException("MoSSE only works with multiple nucleotide alignments. Exiting...");
        }

        // site model should have specific type
        if (what2Check == "sitemodel" && (!(siteModelInput.get() instanceof SiteModel.Base))) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base. Exiting...");
        }
    }
}
