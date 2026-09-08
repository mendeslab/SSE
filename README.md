# SSE


SSE is a [BEAST2](http://beast2.org) package for analyses employing SSE-type models (State-dependent Speciation and Extinction models) in a Bayesian framework.
The SSE package currently implements a well calibrated ClaSSE model (Goldberg and Igic, 2012), with support for its submodels BiSSE, MuSSE and GeoSSE.
We are finalizing the implementation of a stochastic character mapper as in Freyman and Hohna (2018).

## Command-line build with BEAST 2.7

Use JDK 17 or newer and Ant. Point `BEAST2_HOME` at a BEAST 2.7.8 installation
containing `bin/` and `lib/`, then run directly from source (no SSE jar needed):

```sh
export BEAST2_HOME=/path/to/beast
./beast-sse -seed 127 examples/BiSSE_fixed_tree_SDSEP.xml
```

`beast-sse` compiles incrementally and invokes BEAST's standard launcher with normal
package discovery. It preserves your working directory and passes BEAST options unchanged.
Use `JAVA_HOME` to select Java and `JAVA_TOOL_OPTIONS` for JVM memory/profiling options.

`BEAST2_HOME` is required; `ant show-config` prints resolved dependency paths.

To compile without running, use `ant compile`. To run the portable tests:

```sh
ant test
```

The two CLaSSE sampling studies are skipped unless their untracked validation
data are present under `data/`. MoSSE uses a separate native test target; see
`NOTES.md` for its FFTW build requirements and current limitations.

First run `ant fetch-test-deps` to obtain JUnit, or set `junit.jar`.

Do not add the historical `SSE.jar` to this class path: it bundles BEAST 2.6
classes that conflict with BEAST 2.7.
If an installed SSE conflicts with development classes, uninstall it or update its package files.

## Install a package

`ant build` creates `dist/biogeo.v0.0.1.zip` without bundling BEAST. For a local Linux install,
unzip it into `~/.beast/2.7/biogeo/` (remove an older installation first), then use BEAUti's
`File > Clear class path`. Run analyses with `"$BEAST2_HOME/bin/beast" analysis.xml`.
See BEAST's [manual installation instructions](https://www.beast2.org/managing-packages/index.html).

## IDE development

Use `beast.pkgmgmt.launcher.BeastLauncher` as the run configuration's main class.
`ant show-ide` prints the classpath and `-version_file` program arguments without compiling;
replace the SSE `build` entry with your IDE's compiled output, and append the XML path.
Set the working directory to your run directory. Let IntelliJ/Eclipse compile SSE for debugging.
Use dependencies from the configured BEAST installation; sources can be attached in the IDE.

For a short fixed-tree QuaSSE MCMC example, see
[the 15-species example](examples/QuaSSE_15_species_MCMC.md), including its priors,
run commands, validation status, and numerical caveats.
