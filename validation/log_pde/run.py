"""Reproduce the log-PDE study, or regenerate figures/report with --plot-only."""
import argparse
import csv
import gzip
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import make_interp_spline
from scipy.special import logsumexp, ndtr, roots_legendre

from solver import advance, derivatives, log_reference, merge, propagate

HERE = Path(__file__).resolve().parent
DIFFUSION = .2
DURATION = .2
# Fixed physical cases across all grids; off-grid means do not move with resolution.
CASES = {
    "narrow-centered": ([(0, 0, .08 ** 2)], 0),
    "narrow-offset": ([(0, .137, .08 ** 2)], .3),
    "broad": ([(0, .137, .8 ** 2)], .3),
    "overlapping": ([(np.log(.6), -.25, .3 ** 2), (np.log(.4), .35, .45 ** 2)], .3),
    "separated": ([(np.log(.8), -.85, .25 ** 2), (np.log(.2), .85, .35 ** 2)], -.3),
    "merge": ([(np.log(.6), -.25, .3 ** 2), (np.log(.4), .35, .45 ** 2)], .3),
}

MAIN_CASES = tuple(CASES)
for width in (.16, .32):
    CASES[f"width-{width}"] = ([(0, .137, width ** 2)], .3)


# Integrate the exponentiated log interpolant separately from propagation. Increasing
# Gauss-Legendre order checks quadrature error, not interpolation error between grid points.
def metrics(x, h, components, order):
    reference = log_reference(x, components)
    delta = h - reference
    active = reference >= max(reference) + np.log(1e-8)
    scale = max(np.max(reference), np.max(h))
    numerical, exact = np.exp(h - scale), np.exp(reference - scale)
    integral = np.trapezoid
    l1 = integral(abs(numerical - exact), x) / integral(exact, x)
    l2 = np.sqrt(integral((numerical - exact) ** 2, x) / integral(exact ** 2, x))
    interpolant = make_interp_spline(x, h, k=order + 1)
    log_integrals = []
    for count in (8, 16, 32):
        nodes, weights = roots_legendre(count)
        positions = (x[:-1, None] + x[1:, None]) / 2 + np.diff(x)[:, None] / 2 * nodes
        log_values = interpolant(positions) + np.log(weights * np.diff(x)[:, None] / 2)
        log_integrals.append(float(logsumexp(log_values)))
    w, mean, variance = np.asarray(components).T
    mass = ndtr((x[-1] - mean) / np.sqrt(variance)) - ndtr((x[0] - mean) / np.sqrt(variance))
    exact_integral = logsumexp(w + np.log(mass))
    reference_fine = log_reference(positions.ravel(), components).reshape(positions.shape)
    interpolated = interpolant(positions)
    fine_scale = max(np.max(reference_fine), np.max(interpolated))
    fine_exact = np.exp(reference_fine - fine_scale)
    fine_error = np.exp(interpolated - fine_scale) - fine_exact
    quadrature_weights = weights * np.diff(x)[:, None] / 2
    reconstructed_l2 = np.sqrt(np.sum(quadrature_weights * fine_error ** 2)
                               / np.sum(quadrature_weights * fine_exact ** 2))
    return dict(reconstructed_l2=float(reconstructed_l2), log_error_active=float(max(abs(delta[active]))),
                log_error_all=float(max(abs(delta))),
                relative_l1=float(l1), relative_l2=float(l2),
                relative_integral_error=float(abs(np.expm1(log_integrals[-1] - exact_integral))),
                quadrature_change_8_16=float(abs(np.expm1(log_integrals[0] - log_integrals[1]))),
                quadrature_change_16_32=float(abs(np.expm1(log_integrals[1] - log_integrals[2]))),
                mean_log_error_active=float(np.mean(delta[active])))


# Run a single branch or a two-child merge followed by another branch. Every branch
# receives its own analytic exterior data; child numerical errors are retained at the merge.
def experiment(name, dx, extent, order, atol, duration=DURATION):
    components, drift = CASES[name]
    components = np.array(components)
    x = np.linspace(-extent, extent, round(2 * extent / dx) + 1)
    initial = log_reference(x, components)
    output, stats = propagate(x, initial, components, duration, drift, DIFFUSION, order, atol)
    records = [("left" if name == "merge" else "branch", stats)]
    reference = advance(components, duration, drift, DIFFUSION)
    if name == "merge" and output is not None:
        right = np.array([(np.log(.3), -.5, .35 ** 2), (np.log(.7), .4, .25 ** 2)])
        right_output, right_stats = propagate(x, log_reference(x, right), right, duration,
                                              drift, DIFFUSION, order, atol)
        records.append(("right", right_stats))
        output = None if right_output is None else output + right_output + np.log(.15)
        reference = merge(reference, advance(right, duration, drift, DIFFUSION))
        if output is not None:
            output, parent_stats = propagate(x, output, reference, duration, drift, DIFFUSION, order, atol)
            records.append(("parent", parent_stats))
            reference = advance(reference, duration, drift, DIFFUSION)
    row = dict(case=name, dx=dx, extent=extent, order=order, atol=atol, points=len(x),
               drift=drift, diffusion=DIFFUSION, branch_duration=duration,
               status=next((s["status"] for _, s in records if s["status"] != "ok"), "ok"),
               seconds=sum(s["seconds"] for _, s in records),
               evaluations=sum(s["evaluations"] for _, s in records),
               steps=sum(max(0, len(s["times"]) - 1) for _, s in records))
    if output is not None:
        row.update(metrics(x, output, reference, order))
    return row, records


# Small analytic checks protect stencil signs/order and Gaussian-product amplitudes;
# end-to-end convergence alone can miss a drift sign or consistently wrong node scale.
def check_formulas():
    for order in (2, 4, 6):
        x = np.arange(-8, 9) * .125
        radius = order // 2
        interior = x[radius:-radius]
        for degree in range(order + 1):
            slope, curvature = derivatives(x ** degree, .125, order)
            exact1 = degree * interior ** (degree - 1) if degree else np.zeros_like(interior)
            exact2 = degree * (degree - 1) * interior ** (degree - 2) if degree > 1 else np.zeros_like(interior)
            np.testing.assert_allclose(slope, exact1, atol=1e-11)
            np.testing.assert_allclose(curvature, exact2, atol=1e-10)
    left = np.array(CASES["overlapping"][0])
    right = np.array(CASES["separated"][0])
    np.testing.assert_allclose(log_reference(x, merge(left, right)),
                               log_reference(x, left) + log_reference(x, right) + np.log(.15), atol=1e-12)


# Save each completed case immediately. The fixed matrix separates spatial refinement,
# time tolerance, and domain extent; extra tight runs probe the apparent temporal floor.
def measure():
    check_formulas()
    settings = [(name, dx, 4., order, atol, DURATION) for name in MAIN_CASES for dx in (.5, .25, .125, .0625)
                for order in (2, 4, 6) for atol in (1e-6, 1e-9)]
    settings += [(name, .125, extent, 6, 1e-9, DURATION) for name in ("narrow-offset", "separated") for extent in (2., 8.)]
    settings += [(name, .125, 4., 6, 1e-11, DURATION) for name in MAIN_CASES]
    settings += [(f"width-{width}", .125, 4., order, 1e-9, DURATION)
                 for width in (.16, .32) for order in (2, 4, 6)]
    settings += [(name, .125, 4., 6, atol, duration) for name in ("narrow-offset", "separated")
                 for duration in (.02, 1.) for atol in (1e-6, 1e-9)]
    # Finer baselines allow comparisons at matched error, rather than only equal grid size.
    settings += [(name, dx, 4., 2, 1e-9, DURATION) for name in ("overlapping", "separated", "merge")
                 for dx in (.03125, .015625, .0078125)]
    settings += [(name, .03125, 4., order, 1e-9, DURATION)
                 for name in ("overlapping", "separated", "merge") for order in (4, 6)]
    fieldnames = list(experiment("broad", .5, 4., 2, 1e-6)[0])
    with (HERE / "measurements.csv").open("w") as stream, gzip.open(HERE / "steps.csv.gz", "wt") as steps:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        steps_writer = csv.writer(steps)
        steps_writer.writerow(["case", "dx", "extent", "order", "atol", "duration", "branch", "time", "dt"])
        for i, setting in enumerate(settings):
            row, records = experiment(*setting)
            writer.writerow(row)
            stream.flush()
            for branch, stats in records:
                for t, dt in zip(stats["times"][1:], np.diff(stats["times"])):
                    steps_writer.writerow([*setting, branch, t, dt])
            steps.flush()
            print(f"{i + 1}/{len(settings)} {setting}: {row['status']}; "
                  f"{row['evaluations']} evaluations, {row['seconds']:.3f}s", flush=True)
    metadata = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                    platform=platform.platform(), rtol=1e-12, max_evaluations_per_branch=20000,
                    source_sha256={name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
                                   for name in ("solver.py", "run.py")})
    (HERE / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


# Plot only saved measurements. Failed cases are counted in the report and excluded
# from error/cost curves; successful integration does not imply a small spatial error.
def report():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FormatStrFormatter, LogLocator

    with (HERE / "measurements.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        for key in row.keys() - {"case", "status"}:
            row[key] = float(row[key]) if row[key] else float("nan")
    good = [r for r in rows if r["status"] == "ok"]
    base = [r for r in good if r["extent"] == 4 and r["atol"] == 1e-9 and r["branch_duration"] == DURATION]
    for filename, xkey, xlabel in (("spatial-error.png", "dx", "Grid spacing Δx"),
                                   ("cost-error.png", "seconds", "Seconds (single Python evaluation)")):
        fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout="constrained")
        for name, ax in zip(MAIN_CASES, axes.flat):
            for order in (2, 4, 6):
                group = sorted([r for r in base if r["case"] == name and r["order"] == order],
                               key=lambda r: r["dx"])
                ax.loglog([r[xkey] for r in group], [r["reconstructed_l2"] for r in group], "o-", label=f"Order {order}")
            ax.set(title=name, xlabel=xlabel, ylabel="Relative L2 error of log interpolant")
            ax.yaxis.set_major_locator(LogLocator(numticks=12))
            ax.yaxis.set_major_formatter(FormatStrFormatter("%.0e"))
            ax.grid(alpha=.25)
            ax.legend()
        fig.savefig(HERE / filename, dpi=140)
        plt.close(fig)
    with gzip.open(HERE / "steps.csv.gz", "rt") as stream:
        steps = list(csv.DictReader(stream))
    fig, axes = plt.subplots(1, 3, figsize=(15, 4), layout="constrained")
    for name in ("narrow-offset", "broad", "separated"):
        group = [s for s in steps if s["case"] == name and float(s["dx"]) == .125
                 and float(s["extent"]) == 4 and float(s["duration"]) == DURATION
                 and float(s["order"]) == 6 and float(s["atol"]) == 1e-9]
        axes[0].semilogy([float(s["time"]) for s in group], [float(s["dt"]) for s in group], label=name)
    axes[0].set(xlabel="Time", ylabel="Accepted Δt", title="Sixth order, Δx=0.125")
    axes[0].legend()
    for order in (2, 4, 6):
        group = sorted([r for r in base if r["case"] in ("narrow-offset", "width-0.16", "width-0.32", "broad")
                        and r["dx"] == .125],
                       key=lambda r: CASES[r["case"]][0][0][2])
        group = [r for r in group if r["order"] == order]
        axes[1].loglog([np.sqrt(CASES[r["case"]][0][0][2]) for r in group],
                      [r["evaluations"] for r in group], "o-", label=f"Order {order}")
        axes[2].loglog([np.sqrt(CASES[r["case"]][0][0][2]) for r in group],
                      [r["reconstructed_l2"] for r in group], "o-", label=f"Order {order}")
    axes[2].set(xlabel="Initial Gaussian SD", ylabel="Relative L2 error of log interpolant",
                title="Width and error, Δx=0.125")
    axes[2].legend()
    axes[1].set(xlabel="Initial Gaussian SD", ylabel="Derivative evaluations", title="Width and cost, Δx=0.125")
    axes[1].legend()
    fig.savefig(HERE / "timesteps.png", dpi=140)
    plt.close(fig)
    lines = ["# Measured log-PDE results", "", f"Successful runs: {len(good)}/{len(rows)}.", "",
             "Exact exterior values; fixed rtol=1e-12; single-run Python timings.", "",
             "## Spatial errors at atol=1e-9 and domain [−4,4]", "",
             "| Case | Order | Δx=0.5 | 0.25 | 0.125 | 0.0625 | Final observed order |",
             "|---|---:|---:|---:|---:|---:|---:|"]
    for name in MAIN_CASES:
        for order in (2, 4, 6):
            group = sorted([r for r in base if r["case"] == name and r["order"] == order], key=lambda r: -r["dx"])
            errors = {r["dx"]: r["relative_l2"] for r in group}
            rate = np.log2(errors[.125] / errors[.0625]) if .125 in errors and .0625 in errors else float("nan")
            cells = [f"{errors[dx]:.3g}" if dx in errors else "failed" for dx in (.5, .25, .125, .0625)]
            rate_text = f"{rate:.2f}" if name in ("overlapping", "separated", "merge") else "quadratic log; no spatial-order inference"
            lines.append(f"| {name} | {order} | " + " | ".join(cells) + f" | {rate_text} |")
    lines += ["", "## Domain and tight-tolerance checks", "",
              "| Case | Extent | Duration | atol | Grid L2 | Interpolated L2 | Integral error | Quadrature change 16→32 | Evaluations | Seconds |",
              "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in good:
        if r["dx"] == .125 and r["order"] == 6:
            lines.append(f"| {r['case']} | {r['extent']:g} | {r['branch_duration']:g} | {r['atol']:g} | "
                         f"{r['relative_l2']:.3g} | {r['reconstructed_l2']:.3g} | "
                         f"{r['relative_integral_error']:.3g} | {r['quadrature_change_16_32']:.3g} | "
                         f"{r['evaluations']:.0f} | {r['seconds']:.3f} |")
    lines += ["", "## Matched-accuracy comparisons", "",
              "Smallest observed propagation time among the measured grids at atol=1e-9, domain [−4,4],",
              "duration 0.2. Target is relative L2 error of the exponentiated log interpolant.",
              "Timings are single-run measurements, including analytic exterior evaluations; no native comparison.", "",
              "| Case | Target | Order | Δx | Points | Error | Evaluations | Seconds |",
              "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for name in ("overlapping", "separated", "merge"):
        for target in (1e-3, 1e-4, 1e-5):
            for order in (2, 4, 6):
                eligible = [r for r in base if r["case"] == name and r["order"] == order
                            and r["reconstructed_l2"] <= target]
                if not eligible:
                    lines.append(f"| {name} | {target:g} | {order} | — | — | not reached | — | — |")
                    continue
                r = min(eligible, key=lambda r: r["seconds"])
                lines.append(f"| {name} | {target:g} | {order} | {r['dx']:g} | {r['points']:.0f} | "
                             f"{r['reconstructed_l2']:.3g} | {r['evaluations']:.0f} | {r['seconds']:.3f} |")
    for r in rows:
        if r["status"] != "ok":
            lines.append(f"\nFailed: {r['case']}, dx={r['dx']}, order={r['order']}, "
                         f"extent={r['extent']}, atol={r['atol']}: {r['status']}.")
    (HERE / "results.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plot-only", action="store_true")
    args = parser.parse_args()
    if not args.plot_only:
        measure()
    report()
