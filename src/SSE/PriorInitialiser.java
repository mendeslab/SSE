package SSE;

import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.util.Randomizer;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Delegates initialization to a sampleable prior, independently of state-node type.
 * Include all random dependencies in the prior; arbitrary overlapping density factors need not
 * define a directly sampleable joint distribution. Sampling and bounds handling belong to BEAST.
 */
@Description("Initialize state nodes by sampling their existing prior distribution.")
public class PriorInitialiser extends BEASTObject implements StateNodeInitialiser {
    public final Input<State> stateInput = new Input<>("state", "State containing the prior arguments",
            Input.Validate.REQUIRED);
    public final Input<Distribution> priorInput = new Input<>("prior", "Complete sampleable prior",
            Input.Validate.REQUIRED);

    private Random random;

    @Override
    public void initAndValidate() {
        // Resolve for validation only; no state values are changed here.
        argumentNodes();
    }

    // Resolve arguments, not upstream dependencies: fixed hyperparameters must not be reported.
    private List<StateNode> argumentNodes() {
        List<String> arguments = priorInput.get().getArguments();
        if (arguments == null) {
            throw new IllegalArgumentException("Prior must report its argument IDs");
        }
        List<StateNode> nodes = new ArrayList<>();
        for (String id : arguments) {
            StateNode match = null;
            if (id != null && !id.isBlank()) {
                for (StateNode node : stateInput.get().stateNodeInput.get()) {
                    if (id.equals(node.getID())) {
                        match = node;
                        break;
                    }
                }
            }
            if (match == null) {
                throw new IllegalArgumentException("Prior argument '" + id + "' is absent from the state");
            }
            if (!nodes.contains(match)) nodes.add(match);
        }
        return nodes;
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        stateNodes.addAll(argumentNodes());
    }

    // Reset each distribution once, then let BEAST order conditional draws and update the state.
    @Override
    public void initStateNodes() {
        // As in DirectSimulator: some samplers use Random, others use the global Randomizer.
        // Seed lazily and retain the generator so repeated initialization advances both streams.
        if (random == null) random = new Random(Randomizer.getSeed());
        List<BEASTInterface> objects = new ArrayList<>();
        priorInput.get().getPredecessors(objects);
        for (BEASTInterface object : objects) {
            if (object instanceof Distribution distribution) distribution.sampledFlag = false;
        }
        priorInput.get().sample(stateInput.get(), random);
    }
}
