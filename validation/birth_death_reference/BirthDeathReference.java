import SSE.QuaSSEDistribution;
import SSE.QuaSSEGrid;
import beast.base.inference.CompoundDistribution;
import beast.base.inference.MCMC;
import beast.base.parser.XMLParser;
import beast.pkgmgmt.BEASTClassLoader;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/** Evaluate an existing constant-rate, fixed-origin model without starting MCMC. */
public class BirthDeathReference extends QuaSSEDistribution {
    private double minKernelRatio = Double.POSITIVE_INFINITY;

    // Report discretization of the actual Gaussian kernels, without changing production behavior.
    @Override
    public void populatefY(double dt, boolean force, boolean transform, boolean low, boolean jt) {
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
        minKernelRatio = Math.min(minKernelRatio, variance / (grid.getDiffusion() * dt));
    }

    // Require a unique element rather than guessing which of several model components to use.
    private static Element element(Document doc, String xpath) throws Exception {
        var nodes = (org.w3c.dom.NodeList) XPathFactory.newInstance().newXPath()
                .evaluate(xpath, doc, XPathConstants.NODESET);
        if (nodes.getLength() != 1) throw new IllegalArgumentException("Expected one element: " + xpath);
        return (Element) nodes.item(0);
    }

    private static Element id(Document doc, String name) throws Exception {
        return element(doc, "//*[@id='" + name + "']");
    }

    // This diagnostic deliberately supports the coverage-study template, not arbitrary BEAST models.
    private static void expect(Element e, String attribute, String value) {
        if (!e.getAttribute(attribute).equals(value))
            throw new IllegalArgumentException(e.getTagName() + " " + e.getAttribute("id")
                    + ": expected " + attribute + "=" + value);
    }

    // Export the original Newick and prior settings; R independently computes the branching ages.
    private static void inputs(Document doc, Path output) throws Exception {
        Element lik = id(doc, "likelihood");
        expect(lik, "spec", "QuaSSEDistribution");
        for (String[] pair : new String[][] {{"conditionOnRoot", "false"}, {"origin", "@origin"},
                {"conditionOnSampling", "false"}, {"conditionOnRhoSampling", "true"},
                {"fossilSamplingRate", "0"}, {"presentSamplingProbability", "1"}, {"tree", "@observedTree"}})
            expect(lik, pair[0], pair[1]);
        expect(id(doc, "grid"), "drift", "0");
        expect(id(doc, "grid"), "diffusion", "@diffusionVar");
        expect(id(doc, "grid"), "origin", "@origin");
        expect(id(doc, "lambda"), "expression", "diversification 1 turnover - /");
        expect(id(doc, "mu"), "expression", "turnover lambda *");
        expect(id(doc, "diffusionVar"), "expression", "diffusionSD diffusionSD *");
        if (id(doc, "prior").getElementsByTagName("distribution").getLength() != 3)
            throw new IllegalArgumentException("Expected only the three independent parameter priors");
        Element start = element(doc, "//*[@id='likelihood']/startTraitDistribution");
        expect(start, "spec", "Normal");
        if (!Double.isFinite(Double.parseDouble(start.getAttribute("mean")))
                || !(Double.parseDouble(start.getAttribute("sigma")) > 0))
            throw new IllegalArgumentException("Expected a fixed Normal starting-trait distribution");
        for (String[] pair : new String[][] {{"q2mLambda", "@lambda"}, {"q2mMu", "@mu"}}) {
            Element link = element(doc, "//*[@id='likelihood']/" + pair[0]);
            expect(link, "spec", "ConstantLinkFn");
            expect(link, "yV", pair[1]);
        }
        var states = id(doc, "state").getElementsByTagName("stateNode");
        if (states.getLength() != 3) throw new IllegalArgumentException("Expected three scalar state nodes");
        for (String name : new String[] {"diversification", "turnover", "diffusionSD"})
            element(doc, "//*[@id='state']/stateNode[@id='" + name + "']");
        Element d = element(doc, "//*[@id='prior']/distribution[@x='@diversification']/distr");
        Element r = element(doc, "//*[@id='prior']/distribution[@x='@turnover']/distr");
        Element sd = element(doc, "//*[@id='prior']/distribution[@x='@diffusionSD']/distr");
        for (Element e : new Element[] {d, sd}) {
            expect(e, "spec", "LogNormalDistributionModel");
            expect(e, "meanInRealSpace", "false");
        }
        expect(r, "spec", "Beta");
        Files.writeString(output.resolve("tree.nwk"), id(doc, "observedTree").getAttribute("newick") + "\n");
        try (PrintWriter out = new PrintWriter(output.resolve("inputs.tsv").toFile())) {
            out.println("name\tvalue");
            out.println("origin\t" + id(doc, "origin").getAttribute("value"));
            out.println("d_meanlog\t" + d.getAttribute("M"));
            out.println("d_sdlog\t" + d.getAttribute("S"));
            out.println("r_alpha\t" + r.getAttribute("alpha"));
            out.println("r_beta\t" + r.getAttribute("beta"));
            out.println("sd_meanlog\t" + sd.getAttribute("M"));
            out.println("sd_sdlog\t" + sd.getAttribute("S"));
        }
    }

    // Each row starts with a freshly parsed model. Grid padding retains its original timestep
    // while the integration step is refined, so time and padding effects can be distinguished.
    private static void evaluate(Document doc, Path points, Path output) throws Exception {
        id(doc, "likelihood").setAttribute("spec", "BirthDeathReference");
        try (PrintWriter out = new PrintWriter(output.toFile())) {
            out.println("point\tconfig\td\tr\tsd\tlogL\tnX\tdx\tdt\txMin\txMax\tpadLeft\tpadRight"
                    + "\tminKernelRatio\tseconds");
            var rows = Files.readAllLines(points);
            for (String row : rows.subList(1, rows.size())) {
                String[] f = row.split("\t");
                id(doc, "diversification").setAttribute("value", f[2]);
                id(doc, "turnover").setAttribute("value", f[3]);
                id(doc, "diffusionSD").setAttribute("value", f[4]);
                StringWriter xml = new StringWriter();
                TransformerFactory.newInstance().newTransformer().transform(new DOMSource(doc), new StreamResult(xml));
                MCMC run = (MCMC) new XMLParser().parseFragment(xml.toString(), true);
                CompoundDistribution posterior = (CompoundDistribution) run.posteriorInput.get();
                try (BirthDeathReference lik = (BirthDeathReference) posterior.pDistributions.get().get(1)) {
                    QuaSSEGrid g = lik.gridInput.get();
                    int resolution = Integer.parseInt(f[5]), width = Integer.parseInt(f[7]);
                    double timeFactor = Double.parseDouble(f[6]), originalDt = g.getDtMax();
                    if (resolution != 1 || width != 1) {
                        g.nXInput.setValue(g.getBins(true) * resolution * width, g);
                        g.dXInput.setValue(g.getDx() / resolution, g);
                        g.xMidInput.setValue(g.getMidpoint(), g);
                        g.initAndValidate();
                        lik.initAndValidate();
                    }
                    lik.dtMax = originalDt * timeFactor;
                    double start = System.nanoTime();
                    run.startStateInput.get().robustlyCalcPosterior(posterior);
                    double logL = lik.getCurrentLogP();
                    if (!Double.isFinite(logL)) throw new ArithmeticException("Nonfinite likelihood: " + row);
                    if (lik.dtMax != originalDt * timeFactor) throw new IllegalStateException("Timestep overwritten");
                    double[] x = g.copyX(true);
                    int[] pad = g.copyPadding(true);
                    out.printf("%s\t%s\t%s\t%s\t%s\t%.17g\t%d\t%.17g\t%.17g\t%.17g\t%.17g"
                                    + "\t%d\t%d\t%.8g\t%.5f%n",
                            f[0], f[1], f[2], f[3], f[4], logL, g.getBins(true), g.getDx(), lik.dtMax,
                            x[0], x[x.length - 1], pad[0], pad[1], lik.minKernelRatio, (System.nanoTime() - start) / 1e9);
                    out.flush();
                    System.out.println("Completed point " + f[0] + " / " + f[1]);
                }
            }
        }
    }

    // Inputs mode requires no BEAST initialization; evaluation registers the local package explicitly.
    public static void main(String[] args) throws Exception {
        Locale.setDefault(Locale.ROOT);
        if (args.length != 4 && args.length != 6)
            throw new IllegalArgumentException("inputs XML OUTPUT_DIR SSE_DIR | evaluate XML POINTS OUTPUT SSE_DIR BEAST_DIR");
        var factory = DocumentBuilderFactory.newInstance();
        factory.setFeature("http://apache.org/xml/features/disallow-doctype-decl", true);
        Document doc = factory.newDocumentBuilder().parse(Path.of(args[1]).toFile());
        if (args[0].equals("inputs")) {
            inputs(doc, Path.of(args[2]));
        } else if (args[0].equals("evaluate")) {
            BEASTClassLoader.initServices();
            for (String file : new String[] {"BEAST.base.version.xml", "BEAST.app.version.xml"})
                BEASTClassLoader.addServices(Path.of(args[5], "lib/packages", file).toString());
            BEASTClassLoader.addServices(Path.of(args[4], "version.xml").toString());
            BEASTClassLoader.addService("beast.base.core.BEASTInterface", "BirthDeathReference", "biogeo");
            evaluate(doc, Path.of(args[2]), Path.of(args[3]));
        } else throw new IllegalArgumentException("Unknown mode: " + args[0]);
    }
}
