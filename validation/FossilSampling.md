# Fossil sampling implementation and validation

2026-09-30. Fossil tips, sampled ancestors, and present-day sampling are implemented for a fixed
binary tree with a living sample anchoring the present. The user chose equal-weight integration
with fixed useful-bin boundaries in the refinement study, and a parser threshold of 1e−5 for tiny
rounded tip ages. The original measurements and their diagnosis are preserved below.

## Implementation and tests

The [XML example](../examples/QuaSSE_fossils_fixed_tree.xml) and
[model/numerical guide](../examples/QuaSSE_fossils_fixed_tree.md) describe the input contract,
observation events, root conditioning, resolution switch, and numerical failure policy.
No interpolation or shared E trajectory is used. Each terminal fossil independently integrates
E from the present in reusable storage; a sampled ancestor uses its continuing child's E.

Production inventory relative to `utxu` (535452dbf9703ac6fec706769087579b6450a4b6): **77 net lines**
across four production files, versus the plan's rough estimate of 350–600. There are no new
production classes, structs, enums, or result/option types. Two public model inputs, two reusable
E-only workspace sets, and two optional one-row native owners supplement existing node arrays.
Native rate/step coefficients use four additional arrays net of the old exponent arrays. There is
no persistent trajectory, observation map, or alternative tree representation. The shared segment
helpers and private reaction helper replace duplicated code; parameter refresh, tree validation,
and sample initialization account for the other three new helpers.

The three affected test files grow by **295 net lines** (planned estimate 180–300), including one
test-only constant observation class. Reaction references check E and D's analytic multiplier
against a 70-digit matrix exponential, equal/unequal rates, zero rates/duration, and half-step
composition. Native tests cover one/two/three rows, changing ψ and rates, close/reinitialization,
and JNI array validation. Small-tree tests isolate sampling event factors; the resolved Gaussian
case checks the real XML, E transport, and BEAST's state lifecycle without adding a new framework.

`BEAST2_HOME=/home/bredelings/Applications/beast2/beast ant compile-test native` passes.
Final results:

| Run | Passed | Aborted | Failed |
| --- | ---: | ---: | ---: |
| `ant test` (Java/SST) | 74 | 2 | 16 |
| `ant test-native`, first launch (Java integration/FFTW plus FFT/native unit tests) | 30 | 0 | 16 |
| Explicit native distribution/reaction/segment tests | 42 | 0 | 16 |

All 16 failing names and exception types match the recorded pre-feature kernel-resolution
failures. No new failures remain and no kernel guard or test tolerance was weakened.
The first failing launch stops `ant test-native`, so the native-segment distribution tests were
also run explicitly with `-Xcheck:jni`. The two Java aborts are the existing optional native tests.

The finest Gaussian reference is −11.075528805735900 with 1e−8 tolerance. The permanent test parses
the actual XML, verifies living/fossil/ancestor ages, changes ψ, ρ, diffusion and a λ parameter
through BEAST State, compares with fresh distributions, rejects/restores each proposal, and checks
an out-of-bounds ψ proposal. Invalid tree structures and missing sample observations are rejected.
The example's actual 20-step MCMC runs on Java and native with seed 42 have 21 matching trace rows;
the maximum difference across their numeric entries is **8.17e−14**. Logs/state files remain under
`build/fossil-validation/example-java` and `example-native`.

Reaction references were generated with
`python3 build/fossil-validation/reaction-references.py` (mpmath 1.4.1, 70 decimal digits).
Its compact input/output table is retained in `src/test/PropagatesQuaSSETest.java`; the temporary
generator and output remain under build. Scalar tree references use the separate generator below.


## Independent scalar references

Generated using mpmath 1.4.1 with 70 decimal digits. For λ=0.3, μ=0.1, ψ=0.08, ρ=0.7, use
M = [[−(λ+μ+ψ), μ], [−λ, 0]], T = exp(hM),
E(h) = (T11 E(0) + T12)/(T21 E(0) + T22), and J = det(T)/(T21 E(0) + T22)².
Start E(0)=1−ρ; multiply J along each branch, ρ at living tips, ψE at terminal fossils,
ψ at ancestors, and λ at splits. Divide the root contribution by λ[1−E(2)]².
Observation density g=1; trees and grid are those in plan section 9.2.

| Tree | Reference log likelihood |
| --- | ---: |
| Two living tips | −1.30503697077521588852453732830365510993835118579730047 |
| Living tip and terminal fossil | −4.553502530096844596242260350955896024972035398044699016 |
| Living tips and sampled ancestor | −3.83076561508347132830882387329785382091405375997408066 |
| Mixed tree | −8.94450156751828936041078667050598304991597308940224921 |

All match within the planned 1e−8 tolerance on both backends, for both child orders and
for tc=0, 0.25, 0.5, 1.25. The parsed fossil and ancestor ages and ancestor flags are checked.

## Corrected Gaussian refinement: fixed bin boundaries

Keep the biological inputs from the initial Gaussian study below. Only validation grid settings
change; the production Java/C++ integration, root weights, and acceptance tolerances are unchanged.
The useful points represent bin centers. For N useful bins, dX = 7.5/N makes their bins cover
[−3.75, +3.75]; the first and last centers are −3.75+dX/2 and +3.75−dX/2.
Here N = nX − leftPadding − rightPadding − 1. Padding is still computed by production code.

| nX | useful bins N | dX | dtMax | support multiplier | padding per side | log likelihood (native) |
| ---: | ---: | --- | --- | --- | ---: | ---: |
| 2048 | 1919 | 7.5/1919 | 1/64 | sqrt(80) | 64 | −11.075528783765982 |
| 4096 | 3839 | 7.5/3839 | 1/128 | sqrt(160) | 128 | −11.075528796726125 |
| 8192 | 7679 | 7.5/7679 | 1/256 | sqrt(320) | 256 | −11.075528805735900 |

All 24 evaluations (both backends, both child orders, joint/temporal/domain settings) completed
with the ordinary kernel variance guard and `-Xcheck:jni`. The driver asserts expected padding
and useful-bin counts, finite likelihoods, and the outer bin boundaries for every fixed-domain
run. Boundary deviations from ±3.75 were at most 4.29e−13 (tolerance 1e−10).

- Finest joint refinement difference: **9.01e−9**, below the unchanged 1e−4 target.
- Finest temporal refinement difference on the 8192-bin grid: **9.04e−9**, below 1e−4.
- Spatial difference between 4096 and 8192 bins at the same dtMax=1/128: **3.46e−11**.
- Maximum Java/native difference: **1.19e−13**, below 1e−8.
- Maximum child-order difference: **3.55e−15**, below 1e−4.

At 8192 bins, temporal refinement gives −11.075528783718724 at dtMax=1/64,
−11.075528796691480 at 1/128, and −11.075528805735900 at 1/256. Support multipliers are
sqrt(80), sqrt(160), and sqrt(320), keeping the physical kernel support and grid geometry fixed.
The joint-refinement difference is consequently mostly temporal in this experiment.

The separate domain-sensitivity run doubles nX to 16384 at the finest dX, dtMax, and support.
It has outer bin boundaries approximately ±7.750520901159 and log likelihood
−11.803652070174213. This deliberately changes the Flat domain, so it is not a convergence test.

[Raw measurements](fossil-sampling-bin-boundaries.csv) retain all runs, root numerator and
conditioning integral, actual bin centers and boundaries, Fourier width, and elapsed time.
Times are single evaluations, not a warmed performance comparison. The numerator and denominator
include the existing Flat density 1/((nX−1)*dX); its grid-dependent scalar cancels in their ratio.

The local driver is `build/fossil-validation/FossilStudy.java`; its previous version and measurements
are preserved under `build/fossil-validation/fixed-center-endpoints/`. The corrected output is
`build/fossil-validation/continuous-study-bin-boundaries.csv`. Compile/run with:

```sh
export BEAST2_HOME=/home/bredelings/Applications/beast2/beast
ant compile-test native
study_cp="build/fossil-validation:build:build-test:lib/*"
study_cp="$study_cp:$BEAST2_HOME/lib/packages/BEAST.base.jar:$BEAST2_HOME/lib/launcher.jar"
study_cp="$study_cp:build/test-libs/junit-platform-console-standalone-1.8.2.jar"
javac -cp "$study_cp" -d build/fossil-validation build/fossil-validation/FossilStudy.java
java --enable-native-access=ALL-UNNAMED -Xcheck:jni -Djava.library.path=build/native \
    -cp "$study_cp" test.FossilStudy > build/fossil-validation/continuous-study-bin-boundaries.csv
python3 build/fossil-validation/check-bin-boundaries.py
```

The study driver and checks are local investigation scripts, not permanent test infrastructure.
The scalar-tree tests retain their original grid: they do not compare spatial resolutions, and
their constant-rate/constant-observation likelihood is independent of the trait-domain width.

## Initial Gaussian refinement: fixed outermost centers (superseded)

This initial experiment held the first/last useful centers fixed. With equal weights the covered
bin boundaries were ±(3.75+dX/2), changing the domain during refinement. Its apparent first-order
error therefore did not establish a fossil-integration regression. The corrected study above
supersedes its grid prescription and clears the spatial-validation blocker.

Use the mixed tree `(A:2,(F:0.5,(S:0,B:0.75):0.25):1);`, observations A=−0.2, B=0.1,
F=0.3, S=0, observation SD=0.1, λ(x)=0.2+0.1/(1+exp(−x)), μ=0.1, ψ=0.08,
ρ=0.7, diffusion=0.05, drift=0, xMid=0, hiLoRatio=1, tc=0, Flat root weights.
All kernel variance checks pass. Both useful endpoints stay at −3.75 and +3.75 within 1e−10.
The Fourier width is nX*dX and changes slightly with refinement.

| nX | dX | dtMax | support multiplier | padding per side | log numerator | conditioning integral | log likelihood |
| ---: | --- | --- | --- | ---: | ---: | ---: | ---: |
| 2048 | 3.75/959 | 1/64 | sqrt(80) | 64 | −13.2079385527473 | 0.118613296680971 | −11.076051867712062 |
| 4096 | 3.75/1919 | 1/128 | sqrt(160) | 128 | −13.2076615427743 | 0.118615120821201 | −11.075790236505647 |
| 8192 | 3.75/3839 | 1/256 | sqrt(320) | 256 | −13.2075231171417 | 0.118616032882509 | −11.075659500093442 |

These are native values; maximum Java/native log-likelihood difference across all measured
settings is 1.05e−13. Maximum child-order difference is 1.78e−15. Numerators include restored
likelihood scales and the existing Flat density 1/((nX−1)*dX); conditioning integrals use the
same density. They are reported separately to make weighting effects visible.

The last spatial/temporal refinement difference is **1.3073641220451293e−4**, exceeding the
planned 1e−4 threshold. At that stage no tolerance was changed and no permanent Gaussian reference
was added.
At fixed nX=8192 and dX=3.75/3839, temporal refinement gives:

| dtMax | support multiplier (same physical support) | log likelihood |
| --- | --- | ---: |
| 1/64 | sqrt(80) | −11.075659478084665 |
| 1/128 | sqrt(160) | −11.075659491053097 |
| 1/256 | sqrt(320) | −11.075659500093442 |

The last temporal difference is 9.04e−9. This points to spatial discretization rather than
reaction timestep error, but does not identify which spatial operation contributes most.
The roughly halving spatial differences are consistent with a first-order contribution;
root quadrature is a candidate to investigate, not an established diagnosis.

The separate doubled-domain measurement (nX=16384 at unchanged dX and dtMax) gives
endpoints ±7.751041938, log numerator −13.9007313384464, conditioning integral
0.122830643685753, and log likelihood −11.803782585461912. This changes the Flat domain
and is not an integration-convergence failure.

Local reproduction artifacts are under `build/fossil-validation/`: `tree-references.py`,
`continuous-study.csv`, `run-tests.sh`, and the test reports/logs. The initial `FossilStudy.java`
is archived in `fixed-center-endpoints/`; the current driver uses bin boundaries instead.
These investigation scripts are intentionally not tracked as permanent test infrastructure.

## Rounded extant tip ages: resolved at parsing

The rounded 15-species Newick has thirteen nominally living tips with inferred ages between
about 1e−10 and 4e−9; only sp8 and sp9 were exactly zero. The new exact-zero classification exposed
these input rounding errors, initially yielding −∞ at default ψ=0. The supplied QuaSSE XMLs and
that test now use TreeParser `threshold="1e-5"` with `adjustTipHeights="false"`, explicitly setting
small leaf heights to zero. This does not add a tolerance inside the likelihood. Near-present
sampled ancestors below the parser threshold are deferred, as discussed with the user.

The independent R generator `validation/r_scripts/QuaSSEStrangReference.R` retains the original
unadjusted result −52.116649459550381 and additionally reports the threshold-adjusted tree result
**−52.116649450126559**. Java/native agrees within the unchanged 1e−9 test tolerance. R extends
the terminal edges by the computed tip ages without moving internal nodes, matching the parser.
Historical `likelihood_reference2` measurements are unchanged.

The existing test for a zero-age root with two zero-length leaves now expects rejection: that
tree is ambiguous in BEAST's ancestor encoding and excluded by the planned validation.
Positive-height resolution-transition references are unchanged.

## Scope and remaining limitations

This completes the specified first implementation; fossil-only trees, tree proposals, stems,
discrete states, removal at sampling, and trait-dependent/time-dependent sampling remain outside
its scope. Fixed-step remainder handling is preserved for legacy extant-only inputs; new sampling
inputs require dynDt=true. Numerical failures still terminate evaluation rather than silently
altering model support. Performance is measured below; no speed-oriented rewrite is part of this
implementation.

## Warmed native timing comparison

Measured sequentially on an AMD Ryzen 9 7900X, GCC 16.2.0, `-O3`, default architecture flags,
with FFTW/native integration. The local `build/fossil-validation/timing/TimingStudy.java` driver
repeats BEAST posterior evaluation with fixed parameter values, excludes startup, and reports
the likelihood separately from priors. `timing/run.sh` runs each variant in a separate JVM.
[Per-evaluation measurements](fossil-sampling-timing.csv) retain all measured samples.

The primate case uses the existing 233-species Observed reference at 4096 bins, dtMax=H/8000,
padding sized for H/62.5, eight-SD support, and no resolution switch. All variants use the same
input with parser threshold=1e−5. Two warmups precede five measured likelihood evaluations.

| Variant | Median seconds | Range (seconds) |
| --- | ---: | --- |
| Archived pre-feature Java/native implementation | 8.729089 | 8.719501–8.740751 |
| Completed fossil implementation | 10.230472 | 10.223501–10.238248 |
| Temporary current build without per-bin reaction validity checks | 10.061454 | 10.049528–10.094627 |

The complete change is **17.2% slower** in this workload.
The checked version is **1.7% slower** than the otherwise matching
no-check experiment. Removing the checks therefore does not recover most of the baseline speed;
the remaining difference needs profiling before attributing it to particular arithmetic or memory
operations. All three likelihoods agree within 1.3e−11. This is a measured local comparison,
not a hardware-independent speed claim. Production checks were not removed.

The corrected finest fossil example uses five warmups and fifteen measured evaluations:
median **0.119223 seconds**. It cannot be compared with a pre-feature fossil likelihood because
the old implementation did not support fossils. A separate instrumented current build times the
E-only restart/integration/copy for the one terminal fossil: median **0.006423 seconds**, about
**5.1%** of its 0.126129-second total. This includes workspace clearing,
kernel preparation and the copied E boundary, and is specific to this tree's one fossil at age 0.5.
The sampled-ancestor observation requires no separate E trajectory. Instrumentation is confined to
an untracked build copy; the production Java and C++ sources contain no benchmark timing calls.

Inspection confirms that native exponentials, expm1, and square roots are evaluated once per
segment/bin in coefficient preparation, outside the timestep loop. No performance rewrite or
alternative ψ=0 numerical path has been added in response to these measurements.

Library SHA-256 values used in this comparison:

```text
baseline   613afb4cf711091c8cc9bf6a9e0e717e20d4e3ac183e043fb2850dbe1448074b
current    d0e6422311dd15002a38d3c19af28fc954440e6d8c56dd034a2db2dcc545839c
unchecked  a73f8c423b1a1ddc834f64b689e658dc43fb540f17c511e6761111bc6465494e
```
