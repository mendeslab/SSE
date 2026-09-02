package mosse;

import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TraitSet;
import beast.base.evolution.tree.Tree;
import java.lang.UnsupportedOperationException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


/**
 * @author Kylie Chen
 */

@Description("Mosse likelihood class calculates the probability of sequence and trait data on a tree")
public class MosseTreeLikelihood extends TreeLikelihood {

    // trait data
    final public Input<List<TraitSet>> traitListInput = new Input<>("traits", "list of traits", new ArrayList<>());
    // tip model, species diversification model and trait model
    final public Input<MosseTipLikelihood> tipModelInput = new Input<>("tipModel", "model of tip probabilities", Input.Validate.REQUIRED);
    final public Input<Distribution> treeModelInput = new Input<>("treeModel", "species diversification model", Input.Validate.REQUIRED);
    // substitution rate parameters
    final public Input<RealParameter> startSubsRateInput = new Input<>("startSubsRate", "lower range for substitution rate", Input.Validate.REQUIRED);
    final public Input<RealParameter> endSubsRateInput = new Input<>("endSubsRate", "upper range for substitution rate", Input.Validate.REQUIRED);
    final public Input<IntegerParameter> numRateBinsInput = new Input<>("numRateBins", "number of bins for substitution rate", Input.Validate.REQUIRED);

    protected List<TraitSet> traits;
    protected MosseTipLikelihood tipModel;
    protected MosseDistribution treeModel;
    protected double startSubsRate;
    protected double endSubsRate;
    protected int numRateBins;

    @Override
    public void initAndValidate() {
        traits = traitListInput.get();
        tipModel = tipModelInput.get();
        treeModel = (MosseDistribution) treeModelInput.get();

        startSubsRate = startSubsRateInput.get().getValue();
        endSubsRate = endSubsRateInput.get().getValue();
        numRateBins = numRateBinsInput.get().getValue();

        // input checking
        if (dataInput.get().getTaxonCount() != treeInput.get().getLeafNodeCount()) {
            String message = String.format(
                    "The number of leaves in tree (%d) does not match the number of sequences (%d).",
                    treeInput.get().getLeafNodeCount(),
                    dataInput.get().getTaxonCount());
            throw new IllegalArgumentException(message);
        } else if (numRateBins <= 0) {
            throw new IllegalArgumentException("numRateBins input must be a positive integer");
        } else if (!(siteModelInput.get() instanceof SiteModel.Base)) {
            throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        } else if (branchRateModelInput.get() != null) {
            System.err.println("Ignoring clock model " + branchRateModelInput.get().getID());
        }

        beagle = null;

        int nodeCount = treeInput.get().getNodeCount();
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(dataInput.get().getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        // remove requirement for clock model
        branchRateModelInput.setRule(Input.Validate.OPTIONAL);
        branchRateModel = null;
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        int stateCount = dataInput.get().getMaxStateCount();
        int patterns = dataInput.get().getPatternCount();

        // set likelihood core number of states and number of rate bins
        likelihoodCore = new MosseLikelihoodCore(stateCount, numRateBins);

        String className = getClass().getSimpleName();
        Alignment alignment = dataInput.get();

        // logging likelihood class
        Log.info.println(className + "(" + getID() + ") uses " + likelihoodCore.getClass().getSimpleName());
        Log.info.println("  " + alignment.toString(true));

        setProportionInvariant(m_siteModel.getProportionInvariant());
        m_siteModel.setPropInvariantIsCategory(false);
        if (getProportionInvariant() > 0) {
            calcConstantPatternIndices(patterns, stateCount);
        }

        // set up likelihood core and initialize partials
        this.initCore();

        patternLogLikelihoods = new double[patterns];
        m_fRootPartials = new double[patterns * stateCount * numRateBins];
        matrixSize = (stateCount + 1) * (stateCount + 1);
        probabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(probabilities, 1.0);

        if (dataInput.get().isAscertained) {
            useAscertainedSitePatterns = true;
        }
    }

    /**
     * set leaf partials using tip likelihood model *
     */
    @Override
    protected void setPartials(Node node, int patternCount) {
        if (node.isLeaf()) {
            Alignment data = dataInput.get();
            int states = data.getDataType().getStateCount();
            String taxonName = node.getID();
            // get traits
            double[] traitValues = new double[traits.size()];
            for (int i = 0; i < traits.size(); i++) {
                double traitValue = traits.get(i).getValue(taxonName);
                traitValues[i] = traitValue;
            }
            double[] partials = new double[patternCount * states * numRateBins];
            int k = 0;
            int taxonIndex = data.getTaxonIndex(node.getID());
            for (int patternIndex = 0; patternIndex < patternCount; patternIndex++) {
                double[] tipLikelihoods = tipModel.getTipLikelihoods(traitValues, numRateBins, startSubsRate, endSubsRate);
                int stateCount = data.getPattern(taxonIndex, patternIndex);
                    boolean[] stateSet = data.getStateSet(stateCount);
                    for (int state = 0; state < states; state++) {
                        if (stateSet[state] == true) {
                            // set likelihoods for nucleotide in data
                            for (int i = 0; i < numRateBins; i++) {
                                partials[k++] = tipLikelihoods[i];
                            }
                        } else {
                            // otherwise set likelihoods to zero
                            for (int i = 0; i < numRateBins; i++) {
                                partials[k++] = 0.0;
                            }
                        }
                    }
            }
            likelihoodCore.setNodePartials(node.getNr(), partials);

        } else {
            setPartials(node.getLeft(), patternCount);
            setPartials(node.getRight(), patternCount);
        }
    }

    /**
     * set leaf states (not applicable for this class, use setPartials instead)
     */
    @Override
    protected void setStates(Node node, int patternCount) {
        throw new UnsupportedOperationException();
    }

    protected void initCore() {
        final int nodeCount = treeInput.get().getNodeCount();
        likelihoodCore.initialize(
                nodeCount,
                dataInput.get().getPatternCount(),
                m_siteModel.getCategoryCount(),
                true,
                m_useAmbiguities.get()
        );

        final int extNodeCount = nodeCount / 2 + 1;
        final int intNodeCount = nodeCount / 2;

        // set up tip partials
        setPartials(treeInput.get().getRoot(), dataInput.get().getPatternCount());

        hasDirt = Tree.IS_FILTHY;
        for (int i = 0; i < intNodeCount; i++) {
            likelihoodCore.createNodePartials(extNodeCount + i);
        }
    }

    public double calculateLogP() {
        double deltaT = 0.001;
        double rate = 1.0;
        Node node = new Node();
        // BEAST 2.7 substitution models need not expose their state count; the alignment defines the matrix size.
        int numStates = dataInput.get().getMaxStateCount();
        double[] transitionMatrix = new double[numStates * numStates];
        substitutionModel.getTransitionProbabilities(node, 0, deltaT, rate, transitionMatrix);
        double[] vars = new double[1024];
        double[] lambda = new double[1024];
        double[] mu = new double[1024];
        treeModel.calculateBranchLogP(vars, lambda, mu, transitionMatrix);

        return 0.0;
    }

    @Override
    protected boolean requiresRecalculation() {
        hasDirt = Tree.IS_CLEAN;

        if (m_siteModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
        }
        if (tipModel.isDirtyCalculation()) {
            hasDirt = Tree.IS_DIRTY;
            return true;
            // update leaf partials from tip model
        }

        return treeInput.get().somethingIsDirty();
    }
}
