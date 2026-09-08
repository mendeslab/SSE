# Corrected diversitree reference

The nonconstant-extinction references use diversitree 0.10-1 from the sibling checkout
`LSU/diversitree`, commit `6e991378c226e4a3b283236f15c78fe51496a2a4`. This is a local
correction, not an upstream diversitree release. It changes:

- Kernel support to match the backward mean −drift × dt.
- Boundary restoration to preserve nkr bins on the left and nkl on the right: convolution
  reads input[i − offset]. C still restores only E; R restores E and D. The reference cases
  have negligible boundary D, and both backends must agree before a value is used.
- Padding to cover every 0 < dt ≤ dtMax, including the opposing-drift side's interior maximum.
  Split-model extents still maximise over all parameter regimes.
- Root conditioning divides D by the integral `sum(root.p * lambda * (1-E)^2) * dx`.
  The old expression multiplied D by dx after division and introduced an extra dx² factor.
  Conditioned log likelihoods increase by −2 log(dx); unconditioned values are unchanged.

Independent direct-convolution and support tests are in
`diversitree/inst/tests/test-quasse-internal.R`. They fail against the unmodified installed
package and pass against the corrected package. These tests do not use Java as their oracle.

## Rebuild and regenerate

From the SSE repository, with R, FFTW and GSL development dependencies available:

```sh
mkdir -p build/diversitree-root-reference/library build/diversitree-root-reference/source
git -C ../diversitree archive 6e991378c226e4a3b283236f15c78fe51496a2a4 | tar -x -C build/diversitree-root-reference/source
R CMD INSTALL --library=build/diversitree-root-reference/library build/diversitree-root-reference/source
Rscript validation/r_scripts/QuaSSENonconstantMuReference.R build/diversitree-root-reference/library
```

The Git archive command only reads a jj-created commit; all change management remains in jj.
Building from the archive keeps generated native files out of the source checkout. Installation
does not replace the ordinary R package. The generator requires an explicit library path and
prints the loaded package path, version, and a separately labelled expected source revision.
Before calculating reference likelihoods, it checks that conditioning with constant lambda and
zero extinction changes log likelihood by −log(lambda), independently of dx. This rejects the
previous locally corrected package as well as ordinary installations with the normalization bug.
It also checks asymmetric/all-timestep padding and compares
both FFT backends with direct convolution at asymmetric boundaries. These checks reject an ordinary
uncorrected installation even when its end-to-end likelihood comparisons would pass. The checks use
the generator's existing dependencies, not the broader R test suite.

Passing the preflight qualifies the checked behavior; it does not verify the installed package's
exact source revision. Use the archive/build commands to obtain the expected revision. The printed
commit is documented provenance for the reference build, not a detected property of an arbitrary
installed package. The preflight preserves the C/R boundary-policy distinction described above.

The generator retains the original zero-drift check and adds diffusion 0.004 and drift ±0.1,
±1 (diffusion 0.001). Every case uses constant λ and logistic μ, checks C/R agreement within
10⁻¹², and checks stability within 10⁻¹⁰ on doubling nx from 128 to 256 at fixed dx=0.01.
It prints padding and coordinate minima as well as likelihoods. Java checks both fresh and
reused objects against these values within 10⁻⁹, and checks the grid coordinates and transfer map.

## R tests

Install missing test dependencies into the isolated library, without replacing normal packages:

```sh
Rscript -e 'install.packages(c("testthat", "expm", "caper", "lubridate", "minqa"), lib="build/diversitree-root-reference/library", repos="https://cloud.r-project.org")'
Rscript -e '.libPaths(c("build/diversitree-root-reference/library", .libPaths())); testthat::test_dir("../diversitree/inst/tests", filter="quasse", reporter="summary")'
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
