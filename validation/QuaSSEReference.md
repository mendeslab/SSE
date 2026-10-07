# Corrected diversitree reference

The FFT references use diversitree 0.10-1 from the public
[bredelings/diversitree fork](https://github.com/bredelings/diversitree),
commit `7efd65f73068cb32fef3f5966a0e6b297fe2c893` on `strang`.
This is not an upstream release. The stacked `fixes` bookmark supplies:

- Kernel support matching the backward mean −drift × dt.
- Boundary restoration preserving nkr bins on the left and nkl on the right: convolution
  reads input[i − offset]. C restores only E; R restores E and D. The reference cases
  have negligible boundary D, and both backends must agree before a value is used.
- Padding covering every 0 < dt ≤ dtMax, including the opposing-drift side's interior maximum.
- Root conditioning using the integral `sum(root.p * lambda * (1-E)^2) * dx` as denominator.

The child `strang` bookmark changes both FFT integrators to half-reaction / diffusion /
half-reaction steps, matching SSE. The `fixes` version retains the original full-reaction /
diffusion sequence. Comparisons with SSE must use `strang`; matching results establish
implementation agreement at the selected discretization, not numerical convergence.

Independent direct-convolution and support tests are in
`diversitree/inst/tests/test-quasse-internal.R`. They fail against the unmodified installed
package and pass against the corrected package. These tests do not use Java as their oracle.

## Rebuild and regenerate

From the SSE repository, with R, FFTW and GSL development dependencies available:

```sh
git clone https://github.com/bredelings/diversitree.git build/diversitree-reference-repo
mkdir -p build/diversitree-strang-reference/library build/diversitree-strang-reference/source
git -C build/diversitree-reference-repo archive 7efd65f73068cb32fef3f5966a0e6b297fe2c893 | tar -x -C build/diversitree-strang-reference/source
R CMD INSTALL --library=build/diversitree-strang-reference/library build/diversitree-strang-reference/source
Rscript validation/r_scripts/QuaSSENonconstantMuReference.R build/diversitree-strang-reference/library
```

An existing sibling checkout may replace the clone as the archive source. The exact commit,
rather than the movable bookmark name, pins the source. Building the archive keeps generated
files out of the checkout; the isolated library does not replace the ordinary R installation.

The generator requires an explicit library path. Its preflight checks root normalization,
asymmetric/all-timestep padding, and Strang composition against independent pure-birth reaction
formulas and direct convolution. It retains the C/R boundary-policy distinction. A package
containing only `fixes`, without Strang splitting, must fail this preflight.

The printed expected source revision documents provenance; it does not detect the revision
of an arbitrary installed package. Use the pinned archive/build commands to obtain that revision.

The generator retains the original zero-drift check and adds diffusion 0.004 and drift ±0.1,
±1 (diffusion 0.001). Every case uses constant λ and logistic μ, checks C/R agreement within
10⁻¹², and checks stability within 10⁻¹⁰ on doubling nx from 128 to 256 at fixed dx=0.01.
It prints padding and coordinate minima as well as likelihoods. Java checks both fresh and
reused objects against these values within 10⁻⁹, and checks the grid coordinates and transfer map.

## Strang reference

Regenerate the 15-species test with:

```sh
Rscript validation/r_scripts/QuaSSEStrangReference.R build/diversitree-strang-reference/library
```

The generator uses both installed FFT implementations directly, without replacing R functions.
It uses the Java test's rounded tree/traits and grid. TreeParser's `threshold="1e-5"`
extends terminal edges by their tiny implied ages, leaving internal ages fixed; the generator
applies the same adjustment. The adjusted R value is **−52.116649450126559**, used by Java
with a 1e−9 tolerance. The script also prints the unadjusted value to distinguish input-tree
rounding from integration differences.

## R tests

Install missing test dependencies into the isolated library, without replacing normal packages:

```sh
Rscript -e 'install.packages(c("testthat", "expm", "caper", "lubridate", "minqa"), lib="build/diversitree-strang-reference/library", repos="https://cloud.r-project.org")'
Rscript -e '.libPaths(c("build/diversitree-strang-reference/library", .libPaths())); testthat::test_dir("../diversitree/inst/tests", filter="quasse", reporter="summary")'
```

Omit `filter="quasse"` to run the full package suite. On the September 2026 R environment,
the corrected QuaSSE tests pass their new assertions but the old suite cannot pass completely:
method-of-lines tests cannot resolve `initmod_quasse_mol`, and split-model validation uses a
vector in a scalar `||` condition. Both errors also occur with the unmodified installed package.
The full package suite has additional failures outside QuaSSE. Before the root correction, running
the same tests against
the ordinary installed package produced 296 failed expectations and 15 errors, versus 287 failed
expectations and the same 15 errors with the correction: the nine changed outcomes are the new
padding/boundary checks. The four new test groups contain 24 passing assertions in total.
These older failures are not repaired or suppressed by this reference correction.
The C/R FFT generator above is independently runnable
and does not rely on method-of-lines or the package's test helper dependencies.

After the root correction, the QuaSSE suite passes 40 assertions and encounters the same three
errors (two method-of-lines loading errors and one split-model scalar-condition error).
Nine of those assertions are the new analytical root checks; they use synthetic arrays, not PDE
output. The FFT root-mode expectations after the blocked method-of-lines call were also checked
separately, including the unchanged unconditioned likelihood. The generator passes all preflight,
backend, and domain-width checks; the previous reference package fails its normalization preflight.
Java's `ant test` passes 74 tests, with two expected aborts for absent validation data and no failures.

## Fossil sampling references

[FossilSampling.md](FossilSampling.md) records independent matrix-exponential reaction and scalar-tree
references, the fixed-bin-boundary Gaussian refinement study, and Java/native comparisons.
The [fixed-tree fossil guide](../examples/QuaSSE_fossils_fixed_tree.md) gives the equations,
non-removing observation events, sampled-ancestor encoding, present-anchor restriction, and root
conditioning convention. Its Flat example does not replace historical Observed references here
or in `likelihood_reference2`; changing root weights changes the likelihood convention.
