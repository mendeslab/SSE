#!/usr/bin/env python3
r"""Fill a BEAST XML template from fossil-free PhyloJunction simulation outputs.

Requires DendroPy (available in the PhyloJunction Python environment). The template uses IDs
observedTree, observations, mcmc, grid, likelihood and trace. observedTree contains a taxonset;
observations is a keyed RealParameter. Other model elements are preserved.

Example, from SSE (numerical and chain settings here are only an execution example):
  ../PhyloJunction/.venv/bin/python scripts/generate_quasse_xml.py \
      ../quasse_fixed_origin_100 examples/QuaSSE_constant_origin_template.xml /path/to/xmls \
      --bind origin=origin_age --bind observationSD=observation_sd \
      --chain-length 100 --log-every 10 --bins 4096 --dt-max 0.01

Run each XML from its result directory, supplying a distinct BEAST seed per replicate. The output
XMLs contain their data and use relative log filenames. No manifest or subdirectory is created.
"""

import argparse
import copy
import csv
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET

import dendropy


# Read named columns without silently accepting duplicate headers or incomplete records.
def read_rows(path, delimiter, required):
    with path.open(newline="") as stream:
        reader = csv.DictReader(stream, delimiter=delimiter)
        headers = reader.fieldnames or []
        if len(headers) != len(set(headers)) or not set(required) <= set(headers):
            raise ValueError(f"{path}: expected unique columns including {', '.join(required)}")
        rows = list(reader)
    if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError(f"{path}: empty table or incomplete record")
    return rows


# Keep original replicate IDs, permitting only characters that can safely appear in filenames.
def replicate_id(row):
    replicate = row["@replicate"]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", replicate):
        raise ValueError(f"Invalid replicate ID: {replicate!r}")
    if row.get("@element", ""):
        raise ValueError(f"Replicate {replicate}: multiple-element output is not supported")
    return replicate


# Validate the reconstructed tree and observations while retaining the original Newick text.
def validate_tree(newick, observations):
    tree = dendropy.Tree.get(data=newick, schema="newick", preserve_underscores=True)
    labels = [leaf.taxon.label for leaf in tree.leaf_node_iter() if leaf.taxon is not None]
    if len(labels) != len(list(tree.leaf_node_iter())) or len(labels) != len(set(labels)):
        raise ValueError("Tree tips must have unique labels")
    if set(labels) != set(observations):
        raise ValueError("Tree tip labels and sampled observations do not match")
    if any(re.search(r"[\s,]", label) for label in labels):
        raise ValueError("Trait keys cannot contain whitespace or commas")
    distances = {}
    for node in tree.preorder_node_iter():
        length = node.edge_length
        if node.parent_node is not None and length is None:
            raise ValueError("Every non-root branch needs a length")
        if length is not None and (not math.isfinite(length) or length < 0):
            raise ValueError("Branch lengths must be finite and nonnegative")
        if node.parent_node is None:
            distances[node] = 0.0
        else:
            distances[node] = distances[node.parent_node] + length
        if not node.is_leaf() and len(node.child_nodes()) != 2:
            raise ValueError("Reconstructed trees must be binary (or a single tip)")
    heights = [distances[node] for node in tree.leaf_node_iter()]
    # Floating-point summation during simulation can leave tiny differences between living tips.
    if max(heights) - min(heights) > 1e-8:
        raise ValueError("Fossil-free trees must be ultrametric within 1e-8 time units")
    return labels, max(heights)


# Join observations and trees by replicate, reading scalar truth only for explicit fixed bindings.
def read_datasets(directory, columns):
    trees = read_rows(directory / "trees_reconstructed.tsv", "\t", ["@replicate", "@element", "trees"])
    traits = read_rows(directory / "trees_traits.tsv", "\t",
                       ["@replicate", "@element", "node", "alive", "sampled", "observed_trait"])
    datasets = {}
    for row in trees:
        replicate = replicate_id(row)
        if replicate in datasets:
            raise ValueError(f"Duplicate tree for replicate {replicate}")
        datasets[replicate] = {"newick": row["trees"], "observations": {}, "fixed": {}}
    seen = set()
    for row in traits:
        replicate = replicate_id(row)
        if replicate not in datasets:
            raise ValueError(f"Trait replicate {replicate} has no tree")
        seen.add(replicate)
        if row["sampled"] not in ("True", "False") or row["alive"] not in ("True", "False"):
            raise ValueError(f"Replicate {replicate}: expected True/False alive and sampled flags")
        if row["sampled"] == "False":
            continue
        if row["alive"] != "True":
            raise ValueError(f"Replicate {replicate}: fossil observations are not supported")
        observations = datasets[replicate]["observations"]
        if row["node"] in observations:
            raise ValueError(f"Replicate {replicate}: duplicate sampled node {row['node']}")
        value = float(row["observed_trait"])
        if not math.isfinite(value):
            raise ValueError(f"Replicate {replicate}: nonfinite observation")
        observations[row["node"]] = value
    if seen != set(datasets):
        raise ValueError("Tree and trait replicate sets differ")
    if columns:
        rows = read_rows(directory / "combined.csv", ",", ["@replicate", *sorted(set(columns))])
        seen = set()
        for row in rows:
            replicate = replicate_id(row)
            if replicate not in datasets or replicate in seen:
                raise ValueError(f"Unexpected or duplicate scalar replicate {replicate}")
            seen.add(replicate)
            fixed = {column: float(row[column]) for column in columns}
            if not all(math.isfinite(value) for value in fixed.values()):
                raise ValueError(f"Replicate {replicate}: nonfinite fixed value")
            datasets[replicate]["fixed"] = fixed
        if seen != set(datasets):
            raise ValueError("Tree and scalar replicate sets differ")
    for replicate, dataset in datasets.items():
        try:
            dataset["labels"], dataset["height"] = validate_tree(dataset["newick"], dataset["observations"])
        except (ValueError, dendropy.utility.error.DataParseError) as error:
            raise ValueError(f"Replicate {replicate}: {error}") from error
    return datasets


# Resolve IDs once per XML and reject ambiguity before using them as editing targets.
def index_xml(root):
    elements = {}
    for element in root.iter():
        name = element.get("id")
        if name is not None:
            if name in elements:
                raise ValueError(f"Duplicate XML ID: {name}")
            elements[name] = element
    return elements


# Copy and fill only the documented data/execution fields, leaving the model in the template.
def fill_template(template, dataset, bindings, stem, args):
    root = copy.deepcopy(template)
    ids = index_xml(root)
    for name in ("observedTree", "observations", "mcmc", "grid", "likelihood", "trace"):
        if name not in ids:
            raise ValueError(f"Template is missing ID {name}")
    states = list(root.iter("state"))
    state_ids = {node.get("id", node.get("idref")) for state in states for node in state}
    state_ids.update(reference.lstrip("@") for state in states
                     for reference in state.get("stateNode", "").split())
    for target, column in bindings.items():
        element = ids.get(target)
        if (element is None or not element.get("spec", "").endswith("RealParameter")
                or target in state_ids or element.get("dimension", "1") != "1" or len(element)):
            raise ValueError(f"Binding target {target} must be a fixed scalar RealParameter outside the state")
        element.text = None
        element.set("value", repr(dataset["fixed"][column]))
    tree = ids["observedTree"]
    tree.set("newick", dataset["newick"])
    tree.set("adjustTipHeights", "false")
    tree.set("threshold", "1e-8")
    taxa = tree.find("taxonset")
    if taxa is None or taxa.get("idref"):
        raise ValueError("observedTree must contain an inline taxonset")
    old_taxa_ids = {node.get("id") for node in taxa}
    for child in list(taxa):
        taxa.remove(child)
    for label in dataset["labels"]:
        if label in ids and label not in old_taxa_ids:
            raise ValueError(f"Taxon label {label} conflicts with a template XML ID")
        ET.SubElement(taxa, "taxon", id=label, spec="beast.base.evolution.alignment.Taxon")
    observations = ids["observations"]
    observations.text = None
    observations.set("keys", " ".join(dataset["labels"]))
    observations.set("value", " ".join(repr(dataset["observations"][label]) for label in dataset["labels"]))
    observations.set("dimension", str(len(dataset["labels"])))
    origin_ref = ids["grid"].get("origin", "")
    if origin_ref.startswith("@"):
        origin = ids[origin_ref[1:]]
        age = float(origin.get("value", origin.text))
        if not math.isfinite(age) or age < dataset["height"]:
            raise ValueError("Origin age must be finite and at least the reconstructed root height")
    ids["mcmc"].set("chainLength", str(args.chain_length))
    ids["grid"].set("nX", str(args.bins))
    if args.dt_max is not None:
        ids["grid"].set("dtMax", repr(args.dt_max))
    ids["likelihood"].set("fftBackend", args.fft_backend)
    ids["likelihood"].set("integrationBackend", args.integration_backend)
    ids["trace"].set("fileName", f"{stem}.log")
    ids["trace"].set("logEvery", str(args.log_every))
    ET.indent(root, space="  ")
    return ET.tostring(root, encoding="utf-8", xml_declaration=True) + b"\n"


# Validate every replicate and destination before creating the requested output files.
def generate(args):
    bindings = {}
    for binding in args.bind:
        target, separator, column = binding.partition("=")
        if not separator or not target or not column or target in bindings:
            raise ValueError(f"Expected a unique XML_ID=COLUMN binding, got {binding!r}")
        bindings[target] = column
    if args.chain_length < 1 or not 1 <= args.log_every <= args.chain_length:
        raise ValueError("Require chain-length >= log-every >= 1")
    if args.bins < 2 or args.bins & (args.bins - 1):
        raise ValueError("bins must be a power of two greater than one")
    if args.dt_max is not None and (not math.isfinite(args.dt_max) or args.dt_max <= 0):
        raise ValueError("dt-max must be finite and positive")
    template = ET.parse(args.template, parser=ET.XMLParser(target=ET.TreeBuilder(insert_comments=True))).getroot()
    datasets = read_datasets(args.simulation_dir, bindings.values())
    outputs = {}
    for replicate, dataset in datasets.items():
        suffix = replicate.zfill(3) if replicate.isdecimal() else replicate
        path = args.output_dir / f"replicate_{suffix}.xml"
        if path in outputs or path.exists():
            raise ValueError(f"Duplicate or existing output: {path}")
        try:
            outputs[path] = fill_template(template, dataset, bindings, path.stem, args)
        except (ValueError, KeyError, TypeError) as error:
            raise ValueError(f"Replicate {replicate}: {error}") from error
    args.output_dir.mkdir(parents=True, exist_ok=True)
    for path, content in outputs.items():
        with path.open("xb") as stream:
            stream.write(content)
    return len(outputs)


# Keep model choices in XML; command-line options control data bindings and execution only.
def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("simulation_dir", type=Path)
    parser.add_argument("template", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--bind", action="append", default=[], metavar="XML_ID=COLUMN")
    parser.add_argument("--chain-length", type=int, required=True)
    parser.add_argument("--log-every", type=int, required=True)
    parser.add_argument("--bins", type=int, required=True)
    parser.add_argument("--dt-max", type=float,
                        help="override the template timestep; if absent from both, use SSE defaults")
    parser.add_argument("--fft-backend", choices=("sst", "fftw"), default="sst")
    parser.add_argument("--integration-backend", choices=("java", "native"), default="java")
    args = parser.parse_args()
    try:
        count = generate(args)
    except (OSError, ValueError, ET.ParseError) as error:
        parser.exit(1, f"{error}\n")
    print(f"Wrote {count} XML files to {args.output_dir}")


if __name__ == "__main__":
    main()
