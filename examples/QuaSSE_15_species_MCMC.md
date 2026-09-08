# QuaSSE: short fixed-tree MCMC example

`QuaSSE_15_species_MCMC.xml` is a self-contained execution and timing example for BEAST 2.7.8.
It defaults to 20 iterations, not a converged posterior analysis. No R installation is needed to run
the XML. It uses the current Java solver without FFTW or numerical-method changes.

## Data and model

The tree and 15 named trait values come from `diversitree/inst/tests/phy.Rdata`, corresponding to the
small simulated QuaSSE example in `diversitree/doc/diversitree-tutorial.Rnw`. They were exported with
`ape::write.tree(phy, digits=17)` and `sprintf("%.17g", phy$tip.state)`, preserving taxon keys.
Using the saved data avoids changes caused by R's historical random-number-generator defaults.

The tree is fixed. Speciation is logistic, extinction is constant, drift is fixed at zero, and the
shared observation SD is fixed at 0.005. Diffusion is variance per unit time, not an SD. The likelihood
uses the existing `Observed` root treatment and survival conditioning, matching the R example.
There is no separate tree prior, and no tree operator.

| Parameter | Start | Independent prior |
|---|---:|---|
| lambda0 | 0.1 | Exponential, mean 0.2 |
| lambda1 | 0.2 | Exponential, mean 0.2 |
| mu | 0.03 | Exponential, mean 0.05 |
| midpoint | 0 | Normal, mean 0, SD 1 |
| steepness | 2.5 | Lognormal, median 2.5, log-SD 1 |
| diffusion | 0.025 | Exponential, mean 0.025 |

Starting values are the generating parameters, not estimates. These are illustrative priors informed
by the simulation; they are not claimed to be weakly informative. Neither endpoint is constrained
to exceed the other, so increasing and decreasing speciation are both possible. No numerical lower
cutoff has been imposed on diffusion to hide under-resolution of its kernel.

The grid controls reproduce diversitree's defaults for the saved data: nX=1024, dX=5 times the trait
range divided by 1024, xMid at the trait-range midpoint, hiLoRatio=4, flankWidthScaler=5, tc=tree
height/10, and dtMax=tree height/1000. Controls remain fixed, but padding can change with diffusion.
`dynDt=true` subdivides branch intervals; it does not select timesteps from an error estimate.

## Run without an IDE

Build the sibling BEAST repositories as described in the main README, then run these commands from
the SSE checkout. Set `sse_dir` to the absolute path to this checkout if it differs:

```sh
ant compile
sse_dir=/home/bredelings/Devel/LSU/SSE
mkdir -p "$sse_dir/build/quasse-runs/seed127"
cd "$sse_dir/build/quasse-runs/seed127"
java -cp "$sse_dir/build:$sse_dir/../beast2/build/dist/launcher.jar:$sse_dir/../beast2/build/dist/BEAST.base.jar:$sse_dir/../BeastFX/build/dist/BEAST.app.jar:$sse_dir/lib/*" \
  beastfx.app.beast.BeastMCMC -seed 127 "$sse_dir/examples/QuaSSE_15_species_MCMC.xml"
```

For an independent run, use a new directory, such as `seed128`, and `-seed 128`. The absolute build
path lets BEAST discover this checkout's package services. Do not include the historical `SSE.jar`.
Overwrite is deliberately not enabled. The working directory receives `QuaSSE_15_species.log` and
`QuaSSE_15_species_MCMC.xml.state`; keep separate directories to avoid output/state collisions.

The trace records the posterior, likelihood, total prior, six prior contributions, and six parameters
every iteration. Screen output is every 10 iterations; checkpointing is every 10 iterations.
The final operator table reports attempts and acceptances. A zero-attempt operator can have a NaN
acceptance fraction; that alone is not a NaN likelihood. BEAST includes sample 0 after a proposal, so
chainLength=20 produces samples 0 through 20. Its startup "Start likelihood" is the posterior.

To change the length, copy the XML into the run directory and edit `chainLength`, then pass that copy
to BEAST. Increase logging intervals for long runs. Changing the seed does not change the initial
parameter values. A short run checks execution and parameter movement, not convergence, mixing,
parameter recovery, or accuracy throughout the prior/posterior range.

For routine execution checks, start with 20 iterations for each of seeds 127 and 128. Check that
every operator has been attempted, every parameter has moved, and each run has both accepted and
rejected proposals overall. If coverage is incomplete, resume in the same output directory with
the same XML (chainLength=20) and a new explicit seed, for example:

```sh
java -cp "$sse_dir/build:$sse_dir/../beast2/build/dist/launcher.jar:$sse_dir/../beast2/build/dist/BEAST.base.jar:$sse_dir/../BeastFX/build/dist/BEAST.app.jar:$sse_dir/lib/*" \
  beastfx.app.beast.BeastMCMC -resume -seed 129 "$sse_dir/examples/QuaSSE_15_species_MCMC.xml"
```

This restores the parameter checkpoint and appends another batch to the trace. Record the seeds
for each segment; it is not an exact random-number continuation of an uninterrupted run. Stop
when the execution checks are satisfied. Keep this separate from the selected-point numerical
refinement checks below; routine MCMC checks use the original XML grid, not the expensive fine grid.

### Completed short-run checks

With the corrected solver, both runs exercised all six operators and moved all six parameters:

| Initial seed | Batches | Trace samples | Accepted proposals | Rejected proposals | BEAST calculation time |
|---:|---|---|---:|---:|---:|
| 127 | 20 iterations | 0–20 | 20 | 1 | 148.439 s |
| 128 | 20, then resume for 20 with seed 129 | 0–40 | 37 | 5 | 148.236 + 84.014 s |

The first seed-128 batch accepted all 21 proposals, so only that run was extended. Resuming restored
the final posterior exactly and appended sequential trace rows. BEAST also proposes at sample 0
on resumption but suppresses that duplicate log number; hence the second run has 42 proposals and
41 trace rows. All trace entries were finite, and R independently verified every prior contribution,
prior sum, and posterior sum within 1e-12. Both final checkpoints contain the six parameters and
operator statistics. Outputs from these checks are under `build/quasse-mcmc-check/seed127` and
`seed128` (untracked build artifacts).

The initial batches ran concurrently, so their elapsed times are not isolated benchmarks. The
MCMC timings include work beyond one likelihood per reported iteration, including consistency
checks and checkpointing. No tuning changes were made based on these very short runs.

## Validation status and numerical refinement

Compilation and XML initialization passed with BEAST 2.7.8. R independently reproduced all six
logged prior contributions and their sum within 1e-12. At the starting values the log prior is
4.0493197542526271. After correcting root normalization, the Java log likelihood is
-45.556924180196084; the corrected diversitree fftC value is -45.556923977099601.

Both implementations formerly multiplied D by dX after dividing by the root-conditioning sum.
They now divide by the integral (sum times dX). This removes an erroneous dX² likelihood factor.
The error entered diversitree in June 2011, commit `1c44f1c`, when two correct statements were
combined; Java's root calculation included it when added in November 2021. The factor cancelled
in this fixed-tree, fixed-spacing MCMC's ratios, but prevented meaningful raw spacing comparisons.
The correction does not change the priors, grid controls, or branch integration.

These measurements use the corrected Java calculation and an isolated build of diversitree
`6e991378c226e4a3b283236f15c78fe51496a2a4`; see [reference instructions](../validation/QuaSSEReference.md).
They are actual corrected likelihoods, with no offset subtracted afterwards. Only starting-point
likelihoods were evaluated, without MCMC proposals. Times below cover one likelihood evaluation,
excluding model construction, and are indicative measurements rather than a controlled benchmark.

| nX | dX relative to XML | Java log likelihood | R fftC log likelihood | Java seconds | R seconds |
|---:|---:|---:|---:|---:|---:|
| 1024 | 1 | -45.556924180196 | -45.556923977100 | 3.398 | 0.087 |
| 2048 | 1 | -45.556924180196 | -45.556923977099 | 7.233 | 0.175 |
| 2048 | 1/2 | -45.556929374604 | -45.556929171050 | 7.267 | 0.174 |
| 4096 | 1/4 | -45.556937195651 | -45.556936991192 | 15.331 | 0.379 |
| 8192 | 1/8 | -45.556943958240 | -45.556943752914 | 33.066 | 0.786 |

The second row doubles domain width without changing spacing; Java's change is about 8.5e-14.
Subsequent rows halve spacing while doubling bin count, preserving nominal domain width. Their
successive Java changes are about -5.19e-6, -7.82e-6, and -6.76e-6. The large -log(4) jumps are gone,
but these small differences alone do not establish an order of spatial convergence.

At 8192 coarse bins and dX equal to one-eighth of the XML spacing, timestep refinement gives:

| dtMax relative to XML | Java log likelihood | R fftC log likelihood | Java seconds | R seconds |
|---:|---:|---:|---:|---:|
| 1 | -45.556943958240 | -45.556943752914 | 33.066 | 0.786 |
| 1/2 | -45.557060219241 | -45.557060167591 | 65.955 | 1.552 |
| 1/4 | -45.557116927100 | -45.557116914060 | 130.442 | 3.053 |
| 1/8 | -45.557149926821 | -45.557149923472 | 260.801 | 6.070 |

Successive Java changes are approximately -1.16e-4, -5.67e-5, and -3.30e-5. Java/R disagreement
decreases from about 2.05e-7 to 3.35e-9. This is evidence of diminishing timestep sensitivity at
this parameter setting, not an error bound across the posterior. These expensive configurations
are for selected-point accuracy checks, not the default MCMC grid.

The short MCMC runs do not establish convergence or accuracy throughout the parameter space.
In particular, small-diffusion kernel resolution remains outside this correction's scope.
The example is not added to `ant test`.
