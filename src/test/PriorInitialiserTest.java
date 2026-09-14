package test;

import SSE.PriorInitialiser;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.DirectSimulator;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.distribution.Exponential;
import beast.base.inference.distribution.Normal;
import beast.base.inference.distribution.Poisson;
import beast.base.inference.distribution.Prior;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.Assert.*;

public class PriorInitialiserTest {
    // Protect delegation, conditional ordering and flag reset without expensive QuaSSE runs.
    // Replace this coverage if BEAST supplies the initializer itself.
    @Test
    public void delegatesRepeatedConditionalDraws() {
        RealParameter mean = new RealParameter();
        mean.initByName("value", "2");
        mean.setID("mean");
        RealParameter values = new RealParameter();
        values.initByName("value", "0 0");
        values.setID("values");
        IntegerParameter count = new IntegerParameter();
        count.initByName("value", "1");
        count.setID("count");
        State state = new State();
        state.initByName("stateNode", mean, "stateNode", values, "stateNode", count);
        state.initialise();

        Exponential exponential = new Exponential();
        exponential.initByName("mean", "2");
        Prior meanPrior = new Prior();
        meanPrior.initByName("x", mean, "distr", exponential);
        Normal normal = new Normal();
        normal.setID("conditionalNormal");
        normal.initByName("mean", mean, "sigma", "0.1");
        Prior valuesPrior = new Prior();
        valuesPrior.initByName("x", values, "distr", normal);
        Poisson poisson = new Poisson();
        poisson.initByName("lambda", "5");
        Prior countPrior = new Prior();
        countPrior.initByName("x", count, "distr", poisson);
        CompoundDistribution prior = new CompoundDistribution();
        prior.initByName("distribution", valuesPrior, "distribution", meanPrior, "distribution", countPrior);
        PriorInitialiser initializer = new PriorInitialiser();
        initializer.initByName("state", state, "prior", prior);
        List<StateNode> reported = new ArrayList<>();
        initializer.getInitialisedStateNodes(reported);
        assertEquals(List.of(values, mean, count), reported);
        assertEquals(2, mean.getValue(), 0); // Configuration must not draw.

        // Create the thread-local generator before reseeding: BEAST seeds first use differently.
        Randomizer.getSeed();
        Randomizer.setSeed(127);
        double[][] draws = new double[2][];
        for (int i = 0; i < draws.length; i++) {
            initializer.initStateNodes();
            draws[i] = new double[] {mean.getValue(), values.getValue(0), values.getValue(1), count.getValue()};
        }
        assertNotEquals(draws[0][0], draws[1][0], 0);

        Randomizer.setSeed(127);
        Random random = new Random(Randomizer.getSeed());
        DirectSimulator simulator = new DirectSimulator();
        for (double[] expected : draws) {
            mean.setValue(1000.0); // Conditional sampling must replace this before drawing values.
            simulator.clearSampledFlags(prior);
            prior.sample(state, random);
            assertArrayEquals(expected,
                    new double[] {mean.getValue(), values.getValue(0), values.getValue(1), count.getValue()}, 0);
            assertTrue(Math.abs(values.getValue(0) - mean.getValue()) < 1);
        }
        Randomizer.setSeed(127);
        initializer = new PriorInitialiser();
        initializer.initByName("state", state, "prior", prior);
        initializer.initStateNodes();
        assertEquals(draws[0][0], mean.getValue(), 0);
        Randomizer.setSeed(128);
        initializer.initStateNodes();
        assertNotEquals(draws[0][0], mean.getValue(), 0);
    }

    // Missing state arguments must not silently retain placeholders. This check is specific to
    // initializer configuration and can go when equivalent validation is provided by BEAST.
    @Test
    public void rejectsMissingArgument() {
        RealParameter value = new RealParameter("1");
        value.setID("missing");
        Prior prior = new Prior();
        prior.initByName("x", value, "distr", new Normal());
        PriorInitialiser initializer = new PriorInitialiser();
        initializer.stateInput.setValue(new State(), initializer);
        initializer.priorInput.setValue(prior, initializer);
        IllegalArgumentException error = assertThrows(IllegalArgumentException.class, initializer::initAndValidate);
        assertTrue(error.getMessage().contains("missing"));
    }
}
