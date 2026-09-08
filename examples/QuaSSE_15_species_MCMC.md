# QuaSSE: short fixed-tree MCMC example

`QuaSSE_15_species_MCMC.xml` is a self-contained execution and timing example for BEAST 2.7.8.
It defaults to 100 iterations, not a converged posterior analysis. No R installation is needed to run
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
chainLength=100 produces samples 0 through 100. Its startup "Start likelihood" is the posterior.

To change the length, copy the XML into the run directory and edit `chainLength`, then pass that copy
to BEAST. Increase logging intervals for long runs. Changing the seed does not change the initial
parameter values. A short run checks execution and parameter movement, not convergence, mixing,
parameter recovery, or accuracy throughout the prior/posterior range.

## Validation status and numerical caveat

Compilation and XML initialization passed with BEAST 2.7.8. R independently reproduced all six
logged prior contributions and their sum within 1e-12. At the starting values the log prior is
4.0493197542526271; the Java log likelihood is -54.712828027798274 and the installed diversitree
0.10-1 fftC log likelihood is -54.712827824701591 (difference about 2.0e-7).

Starting-point sensitivity checks gave:

| Configuration | Java log likelihood | Change from baseline |
|---|---:|---:|
| XML controls | -54.712828027798274 | 0 |
| Double nX, unchanged dX | -54.712828027798189 | 8.5e-14 |
| Double nX, halve dX | -56.099127583326414 | -1.38629956 |
| Double nX, halve dX and dtMax | -56.099255315211295 | -1.38642729 |

The large spacing-dependent shift also occurs in diversitree. Both implementations' root survival
conditioning divides by a sum and then multiplies by dX, instead of dividing by the integral (sum
times dX). This contributes an extra dX² factor to the likelihood. Halving dX therefore contributes
-log(4). After removing that offset for comparison only, the latter two differences are approximately
-5.2e-6 and -1.33e-4. No correction has been applied to either implementation or this XML.

That factor is constant during this fixed-control, fixed-tree run, so by itself it cancels in MCMC
ratios. It does affect comparisons between grid spacings and absolute likelihoods. Agreement with
diversitree does not establish correctness because both implementations share the expression.

Validation paused at this finding, before the planned two 100-iteration runs. The completed checks
used temporary chainLength=0 copies, which each made one proposal. Thus full-run timing and behavior
are not yet validated. Broader numerical issues, including small-diffusion kernel resolution, remain
outside this example's scope; the example is not added to `ant test`.
