# Fossil sampling implementation: interim validation

2026-09-30. Implementation is paused at the numerical-design review required by
`../TODO/FossilSampling2.md`, section 9.3. No final Gaussian reference or XML example has been
established. Performance measurements and the remaining lifecycle/structure checks are deferred.

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

## Gaussian observation refinement: acceptance criterion failed

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

## Decisions needed before resuming

1. Resolve fossil opt-in versus legacy rounded tip ages. Explicit sampling-mode opt-in would
   preserve old input semantics; alternatively require old trees to declare contemporaneous
   tips when parsed. An age tolerance is possible but introduces a scale-dependent cutoff.
2. Extend the fixed-domain spatial study to a finer grid, or investigate/change spatial
   quadrature separately. A finer study preserves the numerical method but costs more;
   changing quadrature may improve convergence but affects existing reference values.

After those decisions: finish lifecycle/tree-structure checks, rerun old and new tests, establish
an accepted Gaussian reference, add the documented XML example, and measure warmed performance.
