package loggers;

import SSE.StateDependentSpeciationExtinctionProcess;
import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

import java.io.PrintStream;

public class StocCharMapLogger extends BEASTObject implements Loggable {

    public final Input<Tree> treeInput = new Input<>("tree",
            "Tree to log and carry out stochastic character mapping.", Input.Validate.REQUIRED);
    public final Input<StateDependentSpeciationExtinctionProcess> likelihoodInput = new Input<>("likelihood",
            "Model under which to carry out stochastic character mapping.", Input.Validate.REQUIRED);
    public final Input<Boolean> nielsenInput = new Input<>("nielsen",
            "If true, carries out Nielsen-like stochastic character mapping; otherwise maps along branches as well.",
            Input.Validate.REQUIRED);

    private Tree tree;
    private StateDependentSpeciationExtinctionProcess likelihood;
    private boolean isNielsen;

    // Retain the historical logger's single sampled history per logged tree.
    @Override
    public void initAndValidate() {
        tree = treeInput.get();
        likelihood = likelihoodInput.get();
        isNielsen = nielsenInput.get();
    }

    // Draw one history, attach its states to nodes, and emit the annotated tree in NEXUS tree syntax.
    @Override
    public void log(long sample, PrintStream out) {
        int[][] samples = likelihood.sampleStatesForTree(1, isNielsen);
        for (Node node : tree.getNodesAsArray()) {
            node.metaDataString = node.getID() + " state=" + samples[0][node.getNr()];
        }

        out.print("tree STATE_" + sample + " =");
        out.print(tree.getRoot().toSortedNewick(new int[1], true));
        out.print(";");
    }

    @Override
    public void init(PrintStream out) {
        tree.init(out);
    }

    @Override
    public void close(PrintStream out) {
        tree.close(out);
    }
}
