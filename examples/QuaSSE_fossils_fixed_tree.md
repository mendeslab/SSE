# Fixed-tree fossil sampling in QuaSSE

[QuaSSE_fossils_fixed_tree.xml](QuaSSE_fossils_fixed_tree.xml) evaluates the joint likelihood of
continuous observations and a fixed dated tree, conditional on the root split having sampled
descendants on both sides. Samples may be living or fossil on either side. There is no stem
above the root. The 20-step MCMC illustrates inference of sampling parameters; it is not a
converged analysis or a demonstration of tree inference.

Build and run from a directory where BEAST may write its log and state files:

```sh
export BEAST2_HOME=/path/to/beast
cd /path/to/SSE
ant native
export LD_LIBRARY_PATH="$PWD/build/native${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
./beast-sse -seed 42 examples/QuaSSE_fossils_fixed_tree.xml
```

The example selects FFTW/native integration. To use Java only, set `fftBackend="sst"` and
`integrationBackend="java"`. Both produce the initial log likelihood approximately
**−11.0755288057359**; the posterior also includes the sampling-parameter priors.

## Sampling and tree representation

The optional scalar RealParameter inputs are `fossilSamplingRate` (ψ ≥ 0, default 0) and
`presentSamplingProbability` (0 ≤ ρ ≤ 1, default 1). Both must be finite. Fossil sampling is
non-removing: observing a lineage does not terminate it biologically. A terminal fossil has
an unobserved continuation, whereas a sampled ancestor has an observed continuing subtree.
Either explicit sampling input, or any fossil event, requires `dynDt=true`. Extant-only inputs
that omit both parameters keep their previous integration mode. Legacy `dynDt=false` drops
branch-length remainders; its numerical convention is not changed by this feature.

The Newick tree is `(A:2,(F:0.5,(S:0,B:0.75):0.25):1);`:

| Node/event | Age before present |
| --- | ---: |
| Living samples A and B | 0 |
| Terminal fossil F | 0.5 |
| Sampled ancestor S | 0.75 |
| Ordinary internal split | 1 |
| Root split | 2 |

BEAST represents an ancestor observation as a leaf on an exactly zero-length branch from a
binary sampling-event node. The other child is the continuing lineage. The API calls the parent
`isFake()` and the observation leaf `isDirectAncestor()`; the observation event is real.
Every sample label, including S, must have a continuous observation. Unary nodes, multifurcations,
ambiguous pairs of ancestor leaves, present-day ancestor events, and a sampling-event root are
not supported. No tree operators should be added to this fixed-tree example.

**At least one genuinely living sample is required.** BEAST measures heights from the youngest
sample. A fossil-only tree can therefore appear to have a zero-height tip; that does not make
it supported. This implementation has no separate present-age offset.

Use `adjustTipHeights="false"` to preserve fossil ages. The explicit parser `threshold="1e-5"`
sets smaller leaf ages to zero, addressing rounded extant Newick lengths. Units are the tree's
time units; samples younger than that threshold are treated as living. Because BEAST adjusts
leaves without their parents, ancestor observations younger than the threshold would lose their
zero-length connections. Such near-present ancestors are excluded from this example; use a
separately considered age treatment if they occur. `adjustTipHeights="true"` would move all
leaves to the present when date traits are absent and is unsuitable for this mixed tree.

## Likelihood and numerical conventions

With t increasing backward in time, E(x,t) is the probability of no sampled descendants, D(x,t)
is the subtree likelihood, and ℒ is the existing trait drift/diffusion operator:

```text
∂E/∂t = ℒE + μ − (λ+μ+ψ)E + λE²,       E(x,0) = 1−ρ
∂D/∂t = ℒD + [2λE − (λ+μ+ψ)]D
```

Observation density g(y|x) gives these event contributions:

| Event | D on the ancestral side |
| --- | --- |
| Living sample | ρg |
| Terminal fossil | ψgE |
| Sampled ancestor | ψgDcontinuing |
| Ordinary split | λDleftDright |

Sampling scalars are stored as log factors, and branch normalization factors are restored at the
root. E is never normalized. Each terminal fossil independently integrates E from the present
to its age, reusing one workspace per resolution; no shared trajectory or interpolation is used.
Work therefore grows with the sum of terminal-fossil ages divided by dtMax. Ancestor observations
use E from the continuing subtree. All integrations use the same resolution switch and kernel
guard; events exactly at tc use the coarse grid after transfer from the fine grid.

The root likelihood is L = ∫πDroot dx / ∫πλ(1−Eroot)² dx, with Droot including the ordinary
root split and restored scales. The example uses `Flat`; historical `Observed` weighting remains
available and is proportional to Droot, making it data-dependent. Equal-weight root sums use
bin centers. The validation grids keep the outer bin boundaries fixed at ±3.75; grid defaults
and parameter-dependent padding can otherwise change that domain. Flat does not create a fixed
biological prior domain across arbitrary parameter changes.

The local reaction has an exact analytic solution without trait transport; the combined Strang
reaction/transport calculation remains approximate. Derivations are in
[AnalyticSolution.tex](../../TODO/AnalyticSolution.tex) and
[AnalyticSolution2.tex](../../TODO/AnalyticSolution2.tex).
Invalid sampling proposals and identifiable zero-support states return −∞. Numerical failures
(including underresolved Gaussian kernels, nonfinite reaction results, or lost normalization)
still stop evaluation with a diagnostic; they are not silently converted to rejected proposals.
See [the validation report](../validation/FossilSampling.md) for numerical references and limits.
