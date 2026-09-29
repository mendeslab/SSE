"""Run the fixed 233-taxon accuracy study, or regenerate its figures with --plot-only."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BUILD = ROOT / "build/likelihood-reference2"
TABLE = HERE / "measurements.tsv"
DIVISORS = (62.5, 125, 250, 500, 1000, 2000, 4000, 8000)
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
def prepare(previous, libraries):
    beast = Path(os.environ["BEAST2_HOME"])
    input_xml = ROOT / "examples/QuaSSE_233_primates_MCMC.xml"
    files = {"input": input_xml, "driver": HERE / "ReferenceStudy.java", "runner": Path(__file__),
             **libraries}
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
def measure(budget, native_library):
    os.chdir(ROOT)
    libraries = {"Strang": native_library, "Original": ROOT / "build/strang/baseline/libsse_quasse.so"}
    deadline = time.monotonic() + budget
    metadata, rows = read_table()
    done = {key(row) for row in rows}
    if done == set(settings()):
        print(f"All {len(settings())} measurements already saved; regenerating figures only.")
        return
    metadata, cp = prepare(metadata, libraries)
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
               f"-Djava.library.path={libraries[method].parent}", "-cp", cp, "ReferenceStudy",
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


# Prepare one set of logged points and edges for both renderers. Only immediate neighbors
# in the full study grid are connected, so missing, failed and zero-error points leave gaps.
def error_grid(rows):
    bins = sorted({n for method, n, _, width, support in settings() if method == "Strang"})
    divisors = sorted({d for method, _, d, width, support in settings() if method == "Strang"})
    points = {(int(r["n"]), float(r["divisor"])): r for r in rows
              if r["method"] == "Strang" and r["kind"] == "main" and r["status"] == "ok"
              and math.isfinite(r["error"]) and r["error"] > 0}
    coordinates = {key: (math.log10(float(r["dx"])), math.log10(float(r["dtMax"])),
                         math.log10(r["error"])) for key, r in points.items()}
    edges = []
    for i, n in enumerate(bins):
        for j, d in enumerate(divisors):
            neighbors = []
            if i + 1 < len(bins):
                neighbors.append((bins[i + 1], d))
            if j + 1 < len(divisors):
                neighbors.append((n, divisors[j + 1]))
            for neighbor in neighbors:
                if (n, d) in points and neighbor in points:
                    edges.append((coordinates[n, d], coordinates[neighbor]))
    return list(points.values()), list(coordinates.values()), edges


# Render the same measured points and adjacency in a static figure and optional offline HTML.
# Explicitly logged coordinates keep the 3D geometry identical between plotting libraries.
def plot_error_3d(rows, interactive=False):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
    from matplotlib.colors import Normalize

    points, coordinates, edges = error_grid(rows)
    labels = ("log₁₀ Δx", "log₁₀ dtMax", "log₁₀ |log L − log L_reference|")
    caption = ("Error is discrepancy from the numerical reference; Δt denotes dtMax.\n"
               "Lines connect adjacent measured settings; missing, nonfinite and zero-error points are omitted.")
    groups = [("o", "circle", "≥99.9%", []), ("^", "diamond", "95–99.9%", []),
              ("x", "x", "<95%", [])]
    for i, row in enumerate(points):
        retention = float(row["minVariance"])
        groups[0 if retention >= .999 else 1 if retention >= .95 else 2][3].append(i)
    norm = Normalize(min(p[2] for p in coordinates), max(p[2] for p in coordinates))
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(projection="3d")
    ax.add_collection3d(Line3DCollection(edges, colors="0.55", linewidths=.8, alpha=.65))
    for marker, _, name, indices in groups:
        if indices:
            x, y, z = zip(*(coordinates[i] for i in indices))
            ax.scatter(x, y, z, c=z, cmap="viridis", norm=norm, marker=marker, s=42,
                       depthshade=False, label=name)
    ax.set(xlabel=labels[0], ylabel=labels[1], zlabel=labels[2], title="Space and timestep sensitivity — Strang")
    ax.view_init(elev=25, azim=-130)
    ax.legend(title="Minimum variance retention", loc="upper left")
    fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="viridis"), ax=ax, shrink=.55,
                 pad=.13, label="log₁₀ absolute discrepancy")
    fig.subplots_adjust(left=.02, right=.88, bottom=.12, top=.94)
    fig.text(.5, .035, caption, ha="center", fontsize=10)
    fig.savefig(HERE / "error-3d.png", dpi=160)
    plt.close(fig)

    if not interactive:
        return
    import plotly.graph_objects as go

    lines = [[value for edge in edges for value in (edge[0][axis], edge[1][axis], None)]
             for axis in range(3)]
    figure = go.Figure(go.Scatter3d(x=lines[0], y=lines[1], z=lines[2], mode="lines",
                                  line=dict(color="#888", width=2), hoverinfo="skip", showlegend=False))
    for _, symbol, name, indices in groups:
        if not indices:
            continue
        x, y, z = zip(*(coordinates[i] for i in indices))
        details = [[int(points[i]["n"]), float(points[i]["dx"]), float(points[i]["divisor"]),
                    float(points[i]["dtMax"]), points[i]["error"], float(points[i]["minVariance"])]
                   for i in indices]
        figure.add_trace(go.Scatter3d(
            x=x, y=y, z=z, mode="markers", name=name, customdata=details,
            marker=dict(symbol=symbol, size=5, color=z, coloraxis="coloraxis"),
            hovertemplate="%{customdata[0]} bins<br>Δx: %{customdata[1]:.6g}"
                          "<br>Timestep: H/%{customdata[2]}<br>dtMax: %{customdata[3]:.6g}"
                          "<br>Absolute error: %{customdata[4]:.6g}"
                          "<br>Minimum variance retention: %{customdata[5]:.6g}<extra></extra>"))
    figure.update_layout(
        title="Space and timestep sensitivity — Strang", autosize=True,
        scene=dict(xaxis_title=labels[0], yaxis_title=labels[1], zaxis_title=labels[2],
                   aspectmode="cube", camera=dict(eye=dict(x=-1.5, y=-1.8, z=1.2))),
        coloraxis=dict(colorscale="Viridis", cmin=norm.vmin, cmax=norm.vmax,
                       colorbar=dict(title="log₁₀ error", len=.55, y=.4)),
        legend=dict(title="Minimum variance retention", x=0, y=1),
        margin=dict(l=0, r=60, b=100, t=60),
        annotations=[dict(text=caption.replace("\n", "<br>") +
                          "<br>Drag to rotate; scroll to zoom. Diamonds denote 95–99.9% retention.",
                          x=.5, y=-.12, xref="paper", yref="paper", showarrow=False)])
    figure.write_html(HERE / "error-3d.html", include_plotlyjs=True, full_html=True,
                      default_height="95vh", config=dict(displaylogo=False, scrollZoom=True))


# Plot only saved data: the report and figures do not require BEAST or a native library.
def plot(interactive=False):
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
    plot_error_3d(main, interactive)
    print(f"Figures use reference logL={ref:.13f}; its zero self-discrepancy is omitted.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--interactive", action="store_true", help="Also generate offline 3D HTML (requires Plotly)")
    parser.add_argument("--seconds", type=float, default=900, help="Total computation budget; default 900 seconds")
    parser.add_argument("--native-library", type=Path, default=ROOT / "build/native/libsse_quasse.so",
                        help="Current-method library; default build/native/libsse_quasse.so (ignored for --plot-only)")
    args = parser.parse_args()
    if args.interactive:
        try:
            import plotly.graph_objects
        except ImportError:
            parser.error("--interactive requires Plotly: install it with python3 -m pip install plotly")
    if not args.plot_only:
        # Java loads this fixed basename from its directory; hash the same file, retaining symlinks.
        if args.native_library.name != "libsse_quasse.so":
            parser.error("--native-library must point to a file named libsse_quasse.so")
        measure(args.seconds, args.native_library.expanduser().absolute())
    plot(args.interactive)
