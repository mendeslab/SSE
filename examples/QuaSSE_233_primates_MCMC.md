# QuaSSE: 233-primate execution example

`QuaSSE_233_primates_MCMC.xml` runs 20 MCMC iterations with FFTW and native T/X integration.
It is a self-contained fixed-tree example for BEAST 2.7, not a converged primate analysis.
No R installation is needed to run it. It does require the optional QuaSSE native library.

## Data and model

The tree and masses come from `diversitree/doc/data/Vos-2006.nex` and `Redding-2010.csv`,
used in the tutorial's "Primate analysis" section. All 233 taxa are retained. Traits are natural
logarithms of body mass in grams, range 3.881564–11.643954. Taxon names are matched explicitly;
tree lengths and log masses are exported with `ape::write.tree(phy, digits=17)` and
`sprintf("%.17g", log(d$mass))`. The data are not simulated.

The model has logistic speciation, constant extinction, zero drift, and inferred diffusion
(variance per tree-time unit). Observation SD is fixed at 0.02 in log-mass units, as in the
tutorial; that illustrative error assumption may be overconfident. The tree has no operators or
separate prior. Root treatment is `Observed`, with survival conditioning and complete sampling
assumed, matching the basic tutorial setup. Sampling completeness is not established for the data.

| Parameter | Start | Independent prior |
|---|---:|---|
| lambda0 | 0.15 | Exponential, mean 0.2 |
| lambda1 | 0.25 | Exponential, mean 0.2 |
| mu | 0.1 | Exponential, mean 0.1 |
| midpoint | 7.75 | Normal, mean 7.75, SD 2 |
| steepness | 1 | Lognormal, median 1, log-SD 1 |
| diffusion | 0.03 | Exponential, mean 0.05 |

These are illustrative testing priors, not priors calibrated for scientific inference. The midpoint
prior is centered near the observed log-mass range midpoint, so it is data-informed. Steepness is
per log-mass unit. The positive endpoints are unordered, allowing either increasing or decreasing
speciation. Starts are rounded plausible values, not the tutorial's fitted estimates. The midpoint
random-walk window is 0.5; the five scale operators retain scaleFactor=0.75 and equal weights.

Grid controls follow the same defaults as the smaller example: nX=1024, hiLoRatio=4,
flankWidthScaler=5, dX=5 times trait range/1024, xMid=trait-range midpoint, tc=tree height/10,
and dtMax=tree height/250 for Strang splitting. Here height is 65.0916859837, dX=0.0379022945773,
and dtMax=0.260366743935. Padding can change with diffusion; `dynDt=true` is not error-adaptive.
The current default grid can still fail the kernel-variance guard on short branch segments,
including at the initial parameters. Increasing `dtMax` does not resolve those kernels.

## Run

From the SSE checkout, with `BEAST2_HOME` configured:

```sh
make -C jni/quasse BUILD_DIR=../../build/gcc-16 CXXFLAGS='-O3 -g'
export LD_LIBRARY_PATH="$PWD/build/gcc-16${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
sse_dir="$PWD"
run_dir=$(mktemp -d "$PWD/build/primates.XXXXXX")
cd "$run_dir"
JAVA_TOOL_OPTIONS="${JAVA_TOOL_OPTIONS:+$JAVA_TOOL_OPTIONS }--enable-native-access=ALL-UNNAMED" \
    "$sse_dir/beast-sse" -seed 127 "$sse_dir/examples/QuaSSE_233_primates_MCMC.xml"
```

Use a fresh directory and another seed for independent runs. Outputs are `QuaSSE_233_primates.log`
and the XML's state checkpoint. Nothing is overwritten automatically. Copy the XML and change
`chainLength` for longer runs; increase logging/checkpoint intervals as appropriate. The default
logs every iteration, prints every 10, and checkpoints every 10. A 20-step run has 21 trace rows.

For Java integration with FFTW, set `integrationBackend="java"` in a copy. For the portable SST
backend also set `fftBackend="sst"`. Native-library availability never silently changes the backend.

### Prior-initialized pilot

Use `QuaSSE_233_primates_MCMC_2000.xml` in the command above for 2,000 iterations with all six
inferred parameters drawn from their existing priors. XML parameter values are placeholders;
tree, observations, observation SD and zero drift remain fixed. Different `-seed` values give
different starts; use separate run directories. The trace is `QuaSSE_233_primates_MCMC_2000.log`.
An invalid initial posterior stops the run rather than silently retrying another draw.
BEAST's `-resume` restores the checkpoint instead of drawing new starting values. This is a pilot,
not a guarantee of convergence; priors, operators and numerical settings match the short example.

## Validation

All embedded taxon/value pairs and the exported tree were checked against the source files.
Starting-point log prior is 2.386168944446946, independently reproduced in R.
Reference calculations use the isolated corrected diversitree described in
[QuaSSEReference.md](../validation/QuaSSEReference.md); its preflight checks passed.

| Starting-point calculation | Log likelihood |
|---|---:|
| Native T/X, FFTW | -829.9986330762332 |
| Java T/X, FFTW | -829.9986330762331 |
| Corrected diversitree fftC | -829.9986330762316 |
| Corrected diversitree fftR | -829.9986330762260 |
| Native, nX=2048, unchanged dX | -829.9986330762338 |
| Native, nX=2048, half dX | -829.9989777350315 |

A native likelihood took approximately 1.4 seconds on the development machine (release build).
Doubling domain width has negligible effect at the starting point; halving spacing changes log
likelihood by about 0.000345. These checks are not an error bound across the posterior. Small
diffusion can still under-resolve the Gaussian kernel. No solver changes or new numerical cutoffs
were made for this example, and it is not added to `ant test`.

Two sequential native runs completed with 21 finite trace rows each:

| Seed | Accepted | Rejected | BEAST calculation time | Process wall time |
|---:|---:|---:|---:|---:|
| 127 | 16 | 5 | 35.661 s | 37.61 s |
| 128 | 10 | 11 | 35.386 s | 37.32 s |

Both runs attempted all six operators. Diffusion did not move in seed 127; lambda1 did not move
in seed 128. Across the two runs every parameter moved. R independently reproduced all logged
prior contributions, prior sums and posterior sums within 1e-12. Corrected diversitree fftC
matched the two final log likelihoods within 2e-12. These are execution checks, not evidence of
mixing or convergence; operator suggestions from such short runs were not used to retune the XML.

Temporary checks and run outputs are under `build/PrimateReference.R`, `build/PrimateCheck.java`,
`build/PrimateTraceCheck.R`, `build/primate-*.log`, and `build/primate-runs/`.
