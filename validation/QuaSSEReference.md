# Corrected diversitree reference

The nonconstant-extinction references use diversitree 0.10-1 from the sibling checkout
`LSU/diversitree`, commit `5af1bf7bef1742512f3bc84beaf0e6115fc08e42`. This is a local
correction, not an upstream diversitree release. It changes:

- Kernel support to match the backward mean −drift × dt.
- Boundary restoration to preserve nkr bins on the left and nkl on the right: convolution
  reads input[i − offset]. C still restores only E; R restores E and D. The reference cases
  have negligible boundary D, and both backends must agree before a value is used.
- Padding to cover every 0 < dt ≤ dtMax, including the opposing-drift side's interior maximum.
  Split-model extents still maximise over all parameter regimes.

Independent direct-convolution and support tests are in
`diversitree/inst/tests/test-quasse-internal.R`. They fail against the unmodified installed
package and pass against the corrected package. These tests do not use Java as their oracle.

## Rebuild and regenerate

From the SSE repository, with R, FFTW and GSL development dependencies available:

```sh
mkdir -p build/diversitree-reference/library build/diversitree-reference/source-5af1bf7b
git -C ../diversitree archive 5af1bf7bef1742512f3bc84beaf0e6115fc08e42 | tar -x -C build/diversitree-reference/source-5af1bf7b
R CMD INSTALL --library=build/diversitree-reference/library build/diversitree-reference/source-5af1bf7b
Rscript validation/r_scripts/QuaSSENonconstantMuReference.R build/diversitree-reference/library
```

The Git archive command only reads a jj-created commit; all change management remains in jj.
Building from the archive keeps generated native files out of the source checkout. Installation
does not replace the ordinary R package. The generator requires an explicit library path and
prints the loaded package path, version, and a separately labelled expected source revision.
Before calculating reference likelihoods, it checks asymmetric/all-timestep padding and compares
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
Rscript -e 'install.packages(c("testthat", "expm", "caper", "lubridate", "minqa"), lib="build/diversitree-reference/library", repos="https://cloud.r-project.org")'
Rscript -e '.libPaths(c("build/diversitree-reference/library", .libPaths())); testthat::test_dir("../diversitree/inst/tests", filter="quasse", reporter="summary")'
```

Omit `filter="quasse"` to run the full package suite. On the September 2026 R environment,
the corrected QuaSSE tests pass their new assertions but the old suite cannot pass completely:
method-of-lines tests cannot resolve `initmod_quasse_mol`, and split-model validation uses a
vector in a scalar `||` condition. Both errors also occur with the unmodified installed package.
The full package suite has additional failures outside QuaSSE. Running the same tests against
the ordinary installed package produced 296 failed expectations and 15 errors, versus 287 failed
expectations and the same 15 errors with the correction: the nine changed outcomes are the new
padding/boundary checks. The four new test groups contain 24 passing assertions in total.
These older failures are not repaired or suppressed by this reference correction.
The C/R FFT generator above is independently runnable
and does not rely on method-of-lines or the package's test helper dependencies.
