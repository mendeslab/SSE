"""Run the fixed 233-taxon accuracy study, or regenerate its figures with --plot-only."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILD = ROOT / "build/likelihood-reference2"
TABLE = HERE / "measurements.tsv"
DIVISORS = (62.5, 125, 250, 500, 1000, 2000, 4000, 8000)
LIBRARIES = {"Strang": ROOT / "build/gcc-16/libsse_quasse.so",
             "Original": ROOT / "build/strang/baseline/libsse_quasse.so"}
ORIGINAL_HASH = "9bcb4ca968390566e93d791b5da91f643329edd04e70c8b50d8e69e3b3e17959"


# Explicit settings keep the study small: method, bins, timestep divisor, width factor, support SDs.
def settings():
    rows = [("Strang", n, d, 1, 8) for n in (256, 512, 1024, 2048, 4096, 8192) for d in DIVISORS]
    rows += [("Strang", 16384, d, 1, 8) for d in (62.5, 125, 250, 8000, 16000, 32000)]
    rows += [("Original", 4096, d, 1, 8) for d in DIVISORS]
    return rows + [("Strang", 16384, 8000, 2, 8), ("Strang", 8192, 8000, 1, 10)]


def key(row):
    return (row["method"], int(row["n"]), float(row["divisor"]), float(row["width"]), float(row["support"]))


# The table has one JSON comment containing provenance, followed by ordinary tab-separated rows.
def read_table():
    if not TABLE.exists():
        return {}, []
    with TABLE.open() as stream:
        metadata = json.loads(next(stream).removeprefix("# "))
        return metadata, list(csv.DictReader(stream, delimiter="\t"))


# Hash a set of compiled files without filling the report directory with individual manifests.
def combined_hash(paths):
    digest = hashlib.sha256()
    for path in sorted(paths):
        digest.update(str(path.relative_to(ROOT) if path.is_relative_to(ROOT) else path).encode())
        digest.update(b"\0")
        digest.update(path.read_bytes())
    return digest.hexdigest()


# Resume only with unchanged input and executable code; compile outside the report folder.
def prepare(previous):
    beast = Path(os.environ["BEAST2_HOME"])
    input_xml = ROOT / "examples/QuaSSE_233_primates_MCMC.xml"
    files = {"input": input_xml, "driver": HERE / "ReferenceStudy.java", "runner": Path(__file__),
             **LIBRARIES}
    hashes = {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in files.items()}
    hashes["java"] = combined_hash(list((ROOT / "build/SSE").rglob("*.class")))
    hashes["dependencies"] = combined_hash(list((ROOT / "lib").glob("*.jar")) +
                                          [ROOT / "version.xml", beast / "lib/packages/BEAST.base.jar",
                                           beast / "lib/launcher.jar"])
    if hashes["Original"] != ORIGINAL_HASH:
        raise SystemExit("Archived original-method library does not match the source-checked binary.")
    # The runner may extend the matrix or change plots; input, driver and numerical binaries must match.
    if previous and any(previous["sha256"].get(name) != digest
                        for name, digest in hashes.items() if name != "runner"):
        raise SystemExit("Input or code changed. Preserve the existing table before starting a fresh study.")
    metadata = previous or {"revision": subprocess.check_output(
        ["jj", "log", "-r", "@", "--no-graph", "-T", "commit_id"], text=True).strip(),
        "sha256": hashes, "paddingDivisor": 62.5}
    BUILD.mkdir(parents=True, exist_ok=True)
    (BUILD / "input.xml").write_bytes(input_xml.read_bytes())
    cp = (f"{BUILD}:build:lib/*:{beast}/lib/packages/BEAST.base.jar:{beast}/lib/launcher.jar:"
          "build/test-services/BEAST.base/lib:build/test-services/biogeo/lib")
    subprocess.run(["javac", "-cp", cp, "-d", str(BUILD), str(HERE / "ReferenceStudy.java")], check=True)
    return metadata, cp


# Each subprocess releases native memory. Completed rows are appended immediately for resuming.
def measure(budget):
    os.chdir(ROOT)
    deadline = time.monotonic() + budget
    metadata, rows = read_table()
    done = {key(row) for row in rows}
    if done == set(settings()):
        print(f"All {len(settings())} measurements already saved; regenerating figures only.")
        return
    metadata, cp = prepare(metadata)
    for setting in settings():
        if setting in done:
            continue
        if time.monotonic() >= deadline:
            break
        method, n, divisor, width, support = setting
        name = f"{method}-n{n}-d{divisor:g}-w{width}-s{support}"
        directory = BUILD / name
        directory.mkdir(exist_ok=True)
        (directory / "input.xml").write_bytes((BUILD / "input.xml").read_bytes())
        mode = "diagnostic" if n <= 2048 else "reference"
        cmd = ["java", "-Xmx6g", "--enable-native-access=ALL-UNNAMED",
               f"-Djava.library.path={LIBRARIES[method].parent}", "-cp", cp, "ReferenceStudy",
               str(directory), str(n), str(divisor), str(width), str(support), mode, "62.5"]
        print(f"Starting {name} ({len(rows) + 1}/{len(settings())})", flush=True)
        start = time.monotonic()
        try:
            with (directory / "run.log").open("w") as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True,
                               timeout=max(.001, deadline - time.monotonic()))
        except subprocess.TimeoutExpired:
            print("Budget reached; completed measurements are saved.", flush=True)
            break
        elapsed = time.monotonic() - start
        with (directory / "result.tsv").open() as stream:
            row, = csv.DictReader(stream, delimiter="\t")
        row = {"method": method, "kind": "main" if width == 1 and support == 8 else "boundary",
               **row, "wallSeconds": f"{elapsed:.6f}"}
        if row["status"] == "ok":
            assert int(row["completed"]) == 464 and int(row["requests"]) == 464
        with TABLE.open("a") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(row), delimiter="\t", lineterminator="\n")
            if not rows:
                stream.write("# " + json.dumps(metadata, separators=(",", ":")) + "\n")
                writer.writeheader()
            writer.writerow(row)
        rows.append(row)
        done.add(setting)
        print(f"  {row['status']}: logL={row['logL']}, seconds={row['seconds']}, "
              f"minimum variance retention={row['minVariance']}", flush=True)
    print(f"Saved {len(rows)}/{len(settings())} measurements.", flush=True)


# Draw connected measurements; markers identify kernel variance loss, not likelihood accuracy.
def curve(ax, rows, xkey, label):
    rows = sorted(rows, key=lambda r: float(r[xkey]))
    line, = ax.plot([float(r[xkey]) for r in rows], [r["error"] for r in rows], label=label)
    for row in rows:
        retention = float(row["minVariance"])
        marker = "x" if retention < .95 else "^" if retention < .999 else "o"
        ax.scatter(float(row[xkey]), row["error"], color=line.get_color(), marker=marker, s=32)


# Pair each actual numerical coordinate with its generating condition on the same log-axis tick.
def condition_ticks(ax, rows, spatial=False):
    xkey, setting = ("dx", "n") if spatial else ("dtMax", "divisor")
    ticks = sorted({(float(row[xkey]), float(row[setting])) for row in rows})
    labels = [f"{int(value)} bins\n{x:.4g}" if spatial else f"H/{value:g}\n{x:.4g}"
              for x, value in ticks]
    ax.set_xticks([x for x, _ in ticks], labels, fontsize=8)
    ax.tick_params(axis="x", which="minor", labelbottom=False)
    ax.set_xlabel("Uniform grid spacing Δx\nBin count / numerical value" if spatial else
                  "Maximum timestep dtMax\nFraction of tree height / numerical value")


# Plot only saved data: the report and figures do not require BEAST or a native library.
def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    _, rows = read_table()
    references = [r for r in rows if key(r) == ("Strang", 16384, 32000, 1, 8) and r["status"] == "ok"]
    if not references:
        print("Reference not yet available; no error figures generated.")
        return
    reference, = references
    ref = float(reference["logL"])
    main = [r for r in rows if r["kind"] == "main" and r["status"] == "ok"]
    for row in main:
        row["error"] = abs(float(row["logL"]) - ref)
    positive = [r for r in main if r["error"] > 0]
    strang = [r for r in positive if r["method"] == "Strang"]
    figures = []
    fig, ax = plt.subplots(figsize=(13, 6), layout="constrained")
    for n in sorted({int(r["n"]) for r in strang}, reverse=True):
        curve(ax, [r for r in strang if int(r["n"]) == n], "dtMax", f"{n} bins")
    figures.append((fig, [ax], "error-vs-dt.png"))
    ax.set_title("Timestep sensitivity — Strang")
    ax.set_xscale("log")
    condition_ticks(ax, main)
    fig, ax = plt.subplots(figsize=(11, 6), layout="constrained")
    for divisor in (*DIVISORS, 16000):
        curve(ax, [r for r in strang if float(r["divisor"]) == divisor], "dx", f"H/{divisor:g}")
    figures.append((fig, [ax], "error-vs-dx.png"))
    ax.set_title("Spatial sensitivity — Strang")
    ax.set_xscale("log")
    condition_ticks(ax, main, spatial=True)
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), layout="constrained")
    comparison = [r for r in positive if int(r["n"]) == 4096]
    for method in ("Strang", "Original"):
        group = [r for r in comparison if r["method"] == method]
        curve(axes[0], group, "dtMax", method)
        curve(axes[1], group, "seconds", method)
    axes[0].set_xscale("log")
    condition_ticks(axes[0], comparison)
    axes[0].set_title("Timestep sensitivity — 4096 bins")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Seconds per likelihood evaluation")
    axes[1].set_title("Accuracy versus evaluation cost — 4096 bins")
    figures.append((fig, axes, "splitting-comparison.png"))
    for fig, axes, filename in figures:
        for ax in axes:
            ax.set_yscale("log")
            ax.set_ylabel("Absolute log-likelihood discrepancy")
            ax.grid(True, alpha=.25)
            ax.legend(fontsize=9)
        fig.savefig(HERE / filename, dpi=160)
        plt.close(fig)
    print(f"Figures use reference logL={ref:.13f}; its zero self-discrepancy is omitted.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--seconds", type=float, default=900, help="Total computation budget; default 900 seconds")
    args = parser.parse_args()
    if not args.plot_only:
        measure(args.seconds)
    plot()
