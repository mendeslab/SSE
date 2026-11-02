import argparse
import csv
from pathlib import Path
import tempfile
import unittest
import xml.etree.ElementTree as ET

import generate_quasse_xml as converter


# These small end-to-end checks protect joins and template preservation, which BEAST parsing
# alone cannot verify. They remain useful while simulation tables supply inference inputs.
class ConverterTest(unittest.TestCase):
    # Write two replicates in different table orders, with latent and unobserved traits as decoys.
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.directory = Path(self.temp.name)
        self.args = argparse.Namespace(
            simulation_dir=self.directory, output_dir=self.directory / "xml",
            template=Path(__file__).resolve().parents[1] / "examples/QuaSSE_constant_origin_template.xml",
            bind=["origin=origin_age", "observationSD=observation_sd"], chain_length=3, log_every=1,
            bins=4096, dt_max=.01, fft_backend="sst", integration_backend="java")
        self.write("trees_reconstructed.tsv", ["@replicate", "@element", "trees"],
                   [["8", "", "(a:1,b:1);"], ["2", "", "single;"]])
        self.write("trees_traits.tsv", ["@replicate", "@element", "node", "alive", "sampled", "trait", "observed_trait"],
                   [["2", "", "single", "True", "True", "99", "-.2"],
                    ["8", "", "b", "True", "True", "99", ".4"],
                    ["8", "", "ancestor", "False", "False", "99", ""],
                    ["8", "", "a", "True", "True", "99", ".1"]])
        self.write("combined.csv", ["@replicate", "origin_age", "observation_sd", "diversification"],
                   [["2", "4", ".1", "99"], ["8", "5", ".2", "99"]])

    # Keep the test input tables readable while using the same CSV quoting as exported data.
    def write(self, name, columns, rows):
        with (self.directory / name).open("w", newline="") as stream:
            writer = csv.writer(stream, delimiter="\t" if name.endswith("tsv") else ",")
            writer.writerow(columns)
            writer.writerows(rows)

    # Check values, singleton encoding, and reuse with a changed prior and birth link.
    def test_data_and_model_preservation(self):
        root = ET.parse(self.args.template).getroot()
        root.find(".//distr[@M='-0.65']").set("M", "-0.7")
        birth = root.find(".//q2mLambda")
        birth.attrib = dict(spec="LogisticFunction", curveYBaseValue="@lambda", curveMaxY="1.0",
                            sigmoidMidpoint="0", logisticGrowthRate="2")
        self.args.template = self.directory / "model.xml"
        ET.ElementTree(root).write(self.args.template)
        self.assertEqual(2, converter.generate(self.args))
        self.assertEqual({"replicate_002.xml", "replicate_008.xml"},
                         {path.name for path in self.args.output_dir.iterdir()})
        output = ET.parse(self.args.output_dir / "replicate_008.xml").getroot()
        ids = converter.index_xml(output)
        self.assertEqual("5.0", ids["origin"].get("value"))
        self.assertEqual("a b", ids["observations"].get("keys"))
        self.assertEqual("0.1 0.4", ids["observations"].get("value"))
        self.assertEqual("0.52", ids["diversification"].get("value"))
        self.assertEqual("replicate_008.log", ids["trace"].get("fileName"))
        self.assertIsNotNone(output.find(".//distr[@M='-0.7']"))
        self.assertEqual(birth.attrib, output.find(".//q2mLambda").attrib)
        singleton = ET.parse(self.args.output_dir / "replicate_002.xml").getroot()
        self.assertEqual("single;", converter.index_xml(singleton)["observedTree"].get("newick"))
        with self.assertRaisesRegex(ValueError, "existing output"):
            converter.generate(self.args)

    # Invalid mappings must fail before any XMLs are written, including attempts to bind truth
    # into estimated parameters. No likelihood test can detect these data-selection errors.
    def test_bad_inputs(self):
        for bindings in (["diversification=diversification"], ["missing=origin_age"],
                         ["origin=origin_age", "origin=observation_sd"]):
            with self.subTest(bindings=bindings):
                self.args.bind = bindings
                with self.assertRaises(ValueError):
                    converter.generate(self.args)
                self.assertFalse(self.args.output_dir.exists())
        traits_path = self.directory / "trees_traits.tsv"
        traits = traits_path.read_text()
        for bad_traits in (traits.replace("True\tTrue", "False\tTrue"),
                           traits.replace("\tb\t", "\ta\t"), traits.replace("8\t", "9\t")):
            with self.subTest(traits=bad_traits), self.assertRaises(ValueError):
                traits_path.write_text(bad_traits)
                converter.read_datasets(self.directory, [])
        traits_path.write_text(traits)
        for newick in ("(a:1,c:1);", "(a:1,b:2);", "(a:-1,b:1);"):
            with self.subTest(tree=newick), self.assertRaises(ValueError):
                converter.validate_tree(newick, {"a": .1, "b": .2})
        with self.assertRaisesRegex(ValueError, "multiple-element"):
            converter.replicate_id({"@replicate": "1", "@element": "2"})


if __name__ == "__main__":
    unittest.main()
