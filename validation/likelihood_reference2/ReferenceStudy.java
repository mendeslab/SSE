import SSE.*;
import beast.base.inference.*;
import beast.base.parser.XMLParser;
import beast.pkgmgmt.BEASTClassLoader;
import java.nio.file.*;
import java.io.PrintWriter;
import java.util.*;

/** Experimental uniform-grid likelihood measurement; no production algorithms are replaced. */
public class ReferenceStudy extends QuaSSEDistribution {
    private double smallest = Double.POSITIVE_INFINITY, largest = 0, ratio = Double.POSITIVE_INFINITY;
    private long requests, completed, below95, below999;

    // Measure the normalized real-space kernel used by production, before any FFT overwrites it.
    @Override
    public void populatefY(double dt, boolean force, boolean transform, boolean low, boolean jt) {
        requests++;
        smallest = Math.min(smallest, dt);
        largest = Math.max(largest, dt);
        super.populatefY(dt, force, transform, low, jt);
        double[] w = low ? fYLo : fYHi;
        int n = low ? nXbinsLo : nXbinsHi;
        int[] pad = low ? nLeftNRightFlanksLo : nLeftNRightFlanksHi;
        double dx = low ? dXbin : dXbin / hiLoRatio, mean = 0, variance = 0;
        for (int j = -pad[0]; j <= pad[1]; j++) mean += w[2 * (j < 0 ? n + j : j)] * j * dx;
        for (int j = -pad[0]; j <= pad[1]; j++) {
            double delta = j * dx - mean;
            variance += w[2 * (j < 0 ? n + j : j)] * delta * delta;
        }
        double retention = variance / (grid.getDiffusion() * dt);
        ratio = Math.min(ratio, retention);
        completed++;
        if (retention < .95) below95++;
        if (retention < .999) below999++;
    }

    // Parse the frozen biological inputs; vary only numerical controls and evaluate once.
    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.ROOT);
        Path output = Path.of(args[0]);
        int n = Integer.parseInt(args[1]);
        double divisor = Double.parseDouble(args[2]);
        double widthFactor = Double.parseDouble(args[3]), support = Double.parseDouble(args[4]);
        String mode = args.length > 5 ? args[5] : "reference";
        double paddingDivisor = args.length > 6 ? Double.parseDouble(args[6]) : 62.5;
        if (!Double.isFinite(divisor) || !Double.isFinite(paddingDivisor)
                || paddingDivisor <= 0 || divisor < paddingDivisor)
            throw new IllegalArgumentException("Padding timestep must cover the integration timestep");
        if (!mode.equals("reference") && !mode.equals("diagnostic"))
            throw new IllegalArgumentException("Expected reference or diagnostic mode");
        boolean diagnostic = mode.equals("diagnostic");
        BEASTClassLoader.initServices();
        BEASTClassLoader.addServices("version.xml");
        BEASTClassLoader.addService("beast.base.core.BEASTInterface", "ReferenceStudy", "biogeo");
        String xml = Files.readString(output.resolve("input.xml"))
                .replace("spec=\"QuaSSEDistribution\"", "spec=\"ReferenceStudy\"");
        MCMC run = (MCMC) new XMLParser().parseFragment(xml, true);
        CompoundDistribution posterior = (CompoundDistribution) run.posteriorInput.get();
        try (ReferenceStudy lik = (ReferenceStudy) posterior.pDistributions.get().get(1)) {
            QuaSSEGrid g = lik.gridInput.get();
            double height = g.treeInput.get().getRoot().getHeight();
            double dx = g.getDx() * g.getBins(true) * widthFactor / n;
            g.nXInput.setValue(n, g);
            g.dXInput.setValue(dx, g);
            g.hiLoRatioInput.setValue(1, g);
            // tc=0 avoids an identity grid transfer and extra segmentation on uniform grids.
            g.tcInput.setValue(0., g);
            g.flankWidthScalerInput.setValue(support, g);
            g.dtMaxInput.setValue(height / paddingDivisor, g);
            g.initAndValidate();
            lik.initAndValidate();
            // Diagnostic only: permit positive variance loss, retaining normalization/finite checks.
            // The existing input cannot permit exactly zero variance; report that failure separately.
            lik.minimumKernelVarianceRatioInput.setValue(diagnostic ? Double.MIN_VALUE : .999, lik);
            lik.dtMax = height / divisor;
            double[] x = g.copyX(true);
            int[] pad = g.copyPadding(true);
            double start = System.nanoTime(), logL = Double.NaN;
            String status = "ok", failure = "";
            try {
                // BEAST invalidation forces evaluation; the stored value below excludes priors.
                run.startStateInput.get().robustlyCalcPosterior(posterior);
                logL = lik.getCurrentLogP();
                if (!Double.isFinite(logL)) throw new IllegalStateException("Nonfinite likelihood");
            } catch (QuaSSEKernelException e) {
                status = diagnostic ? (e.getMessage().contains("ratio=0.0,")
                        ? "zero_variance_rejected" : "kernel_failure") : "guard_failed";
                failure = e.getMessage().replace('\t', ' ').replace('\n', ' ');
            } catch (IllegalStateException e) {
                if (!"Nonfinite likelihood".equals(e.getMessage())) throw e;
                status = "nonfinite_likelihood";
                failure = e.getMessage();
            }
            double seconds = (System.nanoTime() - start) * 1e-9;
            if (!Arrays.equals(x, g.copyX(true)) || !Arrays.equals(pad, g.copyPadding(true))
                    || lik.dtMax != height / divisor)
                throw new AssertionError("Geometry or timestep changed during evaluation");
            try (PrintWriter out = new PrintWriter(output.resolve("result.tsv").toFile())) {
                out.println("n\tdivisor\twidth\tsupport\tdx\tdtMax\tlogL\tseconds\tminDt\tmaxDt\tminVariance"
                        + "\txMin\txMax\tleftPadding\trightPadding\tstatus\tfailure"
                        + "\tmode\trequests\tcompleted\tbelow95\tbelow999\tpaddingDivisor");
                out.printf("%d\t%s\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.6f\t%.17g\t%.17g\t%.17g"
                                + "\t%.17g\t%.17g\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s%n",
                        n, divisor, widthFactor, support, dx, lik.dtMax, logL, seconds,
                        lik.smallest, lik.largest, lik.ratio, x[0], x[x.length - 1], pad[0], pad[1], status, failure,
                        mode, lik.requests, lik.completed, lik.below95, lik.below999, paddingDivisor);
            }
            System.out.printf("n=%d H/%s width=%g support=%g %s logL=%.13f seconds=%.3f%n",
                    n, divisor, widthFactor, support, status, logL, seconds);
        }
    }
}
