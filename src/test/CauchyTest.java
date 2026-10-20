package test;

import SSE.Cauchy;
import beast.base.inference.parameter.RealParameter;
import org.junit.Assert;
import org.junit.Test;

public class CauchyTest {
    // Analytic quartiles/densities protect parameterization and offset handling; parameter
    // changes ensure the wrapper does not cache stale values after MCMC proposals.
    @Test
    public void distribution() throws Exception {
        RealParameter location = new RealParameter("2"), scale = new RealParameter("3");
        Cauchy prior = new Cauchy();
        prior.initByName("location", location, "scale", scale, "offset", 1.0);
        Assert.assertEquals(-Math.log(3 * Math.PI), prior.logDensity(3), 1e-14);
        Assert.assertEquals(-Math.log(6 * Math.PI), prior.logDensity(6), 1e-14);
        // BEAST's public CDF omits offset; the QuaSSE integrator adjusts the argument itself.
        Assert.assertEquals(.25, prior.getDistribution().cumulativeProbability(-1), 1e-14);
        Assert.assertEquals(6, prior.inverseCumulativeProbability(.75), 1e-14);
        Assert.assertTrue(Double.isNaN(prior.getMean()));
        location.setValue(-2.0); scale.setValue(2.0);
        Assert.assertEquals(-Math.log(2 * Math.PI), prior.logDensity(-1), 1e-14);
        location.setValue(2.0); scale.setValue(3.0);
        Assert.assertEquals(-Math.log(3 * Math.PI), prior.logDensity(3), 1e-14);
        for (String invalid : new String[]{"0", "-1", "NaN", "Infinity", "1 2"}) {
            Cauchy bad = new Cauchy();
            bad.locationInput.setValue(new RealParameter("0"), bad);
            bad.scaleInput.setValue(new RealParameter(invalid), bad);
            Assert.assertThrows(IllegalArgumentException.class, bad::initAndValidate);
        }
    }
}
