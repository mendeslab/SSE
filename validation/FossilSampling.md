# Fossil sampling implementation: interim validation

2026-09-30. The corrected spatial study passes the numerical acceptance criteria in
`../TODO/FossilSampling2.md`, section 9.3. The user chose to retain equal-weight integration and
hold useful bin boundaries fixed during refinement. The initial study and its diagnosis are
preserved below. Fossil implementation remains incomplete: the permanent Gaussian/lifecycle
checks, XML example, legacy tip-age correction, and performance measurements remain to be done.

## Implementation state

The draft model change connects ψ and ρ to Java/native integration, restarts an E-only trajectory
for each terminal fossil, handles sampled ancestors, and tracks sampling factors in log space.
Its parent changes contain the generalized reaction and shared segment integration.

Production inventory for the model change: two new public inputs; no new production classes,
structs, or enums; three private helpers (sampling parameter refresh, tree validation, and sample
initialization). Two reusable E-only workspace sets and two optional native owners supplement the
existing node arrays. There is no stored E trajectory or separate observation map.
The model file grows by 152 lines net; its existing test file grows by 131 lines net, including
one test-only constant observation class and a model construction helper.

Build: `BEAST2_HOME=/home/bredelings/Applications/beast2/beast ant compile-test native` passes.
The new scalar-tree and sampling-support test groups pass under Java/SST and native/FFTW with
`-Xcheck:jni`. The broader native distribution suite has 21 passing tests and 17 failures:
16 unchanged kernel-resolution failures and the additional compatibility failure below.

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
planned 1e−4 threshold. No tolerance was changed and no new permanent Gaussian reference was added.
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
`FossilStudy.java`, `continuous-study.csv`, `run-tests.sh`, and the test reports/logs.
These investigation scripts are intentionally not tracked as permanent test infrastructure.

## Legacy rounded tree: compatibility assumption failed

`testPruneFifteenSpTree1024Bins` now returns −∞ instead of its existing −52.11664945955038.
Its Newick branch lengths are rounded and TreeParser is called with tip-height adjustment off.
Thirteen nominally living tips therefore have positive heights between about 1e−10 and 4e−9;
only sp8 and sp9 have exactly zero height. Under the specified positive-height classification,
the new code treats these as fossils. The omitted ψ input defaults to zero, giving zero support.

This is a conflict between automatic fossil detection and preserving legacy extant-only inputs.
Simply applying a time tolerance would change the supported age convention and could erase
real young fossils. No such tolerance or compatibility rule has been introduced.

Separately, the existing test for a zero-age root with two zero-length leaves now expects
rejection: that tree is ambiguous in BEAST's ancestor encoding and is explicitly excluded by
the planned tree validation. Positive-height resolution-transition references are unchanged.

## Remaining implementation work

The user resolved the spatial convention by choosing fixed bin boundaries and unchanged
production quadrature; the rerun passes all planned numerical acceptance criteria. Update the
legacy tree inputs with the agreed parser treatment of rounded extant tip ages separately.
Then finish lifecycle/tree-structure checks, rerun old and new tests, retain the accepted Gaussian
reference, add the documented XML example, and measure warmed performance.
