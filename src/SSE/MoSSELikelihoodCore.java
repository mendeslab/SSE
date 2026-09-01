package SSE;

/*
 * This class is the analog of 'LikelihoodCore' in beast.base.evolution.likelihood.
 *
 * Note that unlike 'LikelihoodCore' or any of its daughter classes such
 * as BeerLikelihoodCore, we do not manage all site patterns. We only look at
 * a single site pattern at a time. The loop over site patterns is in the class
 * calling 'MosseLikelihoodCore'.
 *
 * They key difference here is that we have 1 partial array per quantitative trait
 * bin (i.e., there are 'nXbins' partials per nucleotide site). And we don't update
 * these partials (from the substitution model perspective) at internal nodes, but
 * only along branches instead.
 *
 * This means one partial goes in, one partial 'dt' time units before comes out.
 *
 */
public class MoSSELikelihoodCore {

    private int nrOfStates;
    protected int matrixSize;

    protected double[] matrices;

    public MoSSELikelihoodCore(int nrOfStates, double[] matrices) {
        this.nrOfStates = nrOfStates;
        this.matrices = matrices;
    } // ctor

    public void updateMatrices(double[] newMatrices) {
        matrices = newMatrices;
    }

    // with 1 E and 4 D's, eSDsAtNodeAtBin should have dimension 5
    public void updatePartialInPlaceFirstDt(double[] esDsAtNodeAtBin, double dt) {
        ;
    }

    // with 1 E and 4 D's, eSDsAtNodeAtBin should have dimension 5
    public void updatePartialInPlace(double[] esDsAtNodeAtBin, double dt) {
        ;
    }
}
