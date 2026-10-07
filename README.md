# SSE


SSE is a [BEAST2](http://beast2.org) package for analyses employing SSE-type models (State-dependent Speciation and Extinction models) in a Bayesian framework.
The SSE package currently implements a well calibrated ClaSSE model (Goldberg and Igic, 2012), with support for its submodels BiSSE, MuSSE and GeoSSE.
We are finalizing the implementation of a stochastic character mapper as in Freyman and Hohna (2018).

## Mathematical documentation

See [the mathematical documentation](docs/README.md) for the QuaSSE guide, derivations,
and ongoing numerical explorations, including build instructions and contribution conventions.

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
The [233-primate example](examples/QuaSSE_233_primates_MCMC.md) uses FFTW/native integration.
The [fixed-tree fossil example](examples/QuaSSE_fossils_fixed_tree.md) includes terminal fossils,
sampled ancestors, and inference of the scalar sampling parameters ψ and ρ. It requires a living
sample to anchor the present and `dynDt=true`; it uses non-removing sampling and Flat root weights.
The linked guide explains the event equations, BEAST encoding, root conditioning, numerical domain,
and the repeated E-integration cost. Supplied QuaSSE XMLs use a parser tip-age threshold of 1e−5
for rounded living-tip lengths; mixed fossil trees must retain `adjustTipHeights="false"`.

QuaSSE can omit `grid` when `q2d` is `NormalCenteredAtObservedLinkFn`: supply likelihood-level
`drift` and `diffusion`, and defaults use the initial tree and observed traits. Reinitialize the
likelihood to resolve these defaults again. Other links or numerical overrides require an explicit
`QuaSSEGrid`, taking `tree`, `traits`, `drift`, and `diffusion`. Optional scalar
attributes `nX`, `dX`, `xMid`, `hiLoRatio`, `flankWidthScaler`, `rangeMultiplier`, `dtMax`, and `tc`
override defaults resolved at initialization. With Strang splitting, `dtMax` defaults to the initial
tree height / 250; an explicit value overrides it. Reinitialize the grid and likelihood
after editing these configuration values. Parameter proposals use BEAST's invalidation lifecycle;
standalone callers should likewise notify through `State` before evaluating changed inputs.
If drift/diffusion are supplied both on the likelihood and grid, they must reference the same objects.

QuaSSE requires each sampled Gaussian kernel to retain at least 95% of its intended variance.
Set `minimumKernelVarianceRatio` on the distribution to change this threshold (finite, in (0, 1]).
Insufficient variance or invalid normalization throws `QuaSSEKernelException`, stopping evaluation
rather than rejecting an MCMC proposal. Its message identifies the timestep and grid; finer spacing
may help. Short branch segments can fail even with a large `dtMax`. This checks variance deficiency,
not drift accuracy or total likelihood error. There is no automatic grid refinement or retry.

## Optional QuaSSE FFTW transforms (Linux)

Requires a C++20 compiler, FFTW development files, `pkg-config`, and JDK headers
(`JAVA_HOME` if not using the `javac` on PATH). MoSSE's native library is separate.
Compilation and the FFT/native-integrator unit tests have been checked with GCC 13 and GCC 16;
other C++20 compilers have not been verified.

```sh
ant test-native
export LD_LIBRARY_PATH="$PWD/build/native${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
```

Set `fftBackend="fftw"` on the XML's `QuaSSEDistribution`, then run with `beast-sse` as usual.
The default is `sst`; selecting FFTW requires the library and never falls back silently.
For an IDE run, set `-Djava.library.path=/path/to/SSE/build/native` in JVM options.
On newer JDKs, add `--enable-native-access=ALL-UNNAMED` to JVM options (or `JAVA_TOOL_OPTIONS`).
By default only transforms move into C++; the numerical method is unchanged.

Both integration backends use Strang splitting: half a reaction step, a full drift/diffusion
step, then half a reaction step. This does not remove spatial or boundary error or
automatically select a timestep.

To run repeated steps in C++, set `integrationBackend="native"` alongside
`fftBackend="fftw"` on the QuaSSE distribution. The default is `integrationBackend="java"`;
Java still selects time steps and handles grids, tree operations, and normalization.
Rebuild with `ant native` after updating the JNI interface.

Native builds default to release mode. Each mode selects its own directory and default flags:

| Mode | Directory | Default compiler flags |
|---|---|---|
| `release` | `build/native` | `-O3` |
| `debug-O` | `build/native-debug-O` | `-O -g` |
| `debug` | `build/native-debug` | `-g` |

Select a mode with `ant native -Dnative.mode=debug-O` or
`ant test-native -Dnative.mode=debug-O`. Native tests build and load that same configuration.
To run an analysis with it, use its directory in `LD_LIBRARY_PATH` or `java.library.path`.

Make uses the default C++ compiler (`g++`), or `CXX` from the environment or command line.
For example, `CXX=g++-13 ant native` selects GCC 13. `CXXFLAGS` replaces the mode's default flags;
`CPPFLAGS`, `LDFLAGS`, and `LDLIBS` supply additional preprocessor, linker, and library options.
Required C++20, shared-library, JNI, and FFTW options remain part of the build command.
No architecture flags are added automatically: use, for example,
`CXXFLAGS='-O3 -march=ivybridge' ant native` to select one explicitly. These flags affect SSE's
C++ code, not the installed FFTW library. `-march=native` targets the build machine's CPU,
which may differ from the compute nodes where the library will run.

Direct Make builds accept the same modes, for example `make -C jni/quasse MODE=debug`.
For a custom output directory, Make resolves relative paths from `jni/quasse`, while Ant resolves
them from the project root:

```sh
make -C jni/quasse MODE=debug CXX=g++ BUILD_DIR=../../build/custom
CXX=g++ ant test-native -Dnative.mode=debug -Dnative.build=build/custom
```

Make tracks source timestamps, not compiler selection, Java paths, or flags. Clean and rebuild
before changing those settings within an existing directory. Switching between modes requires
no cleaning because their default directories are separate. Repeat the same overrides when
running tests so any rebuild uses the intended settings.

```sh
ant clean-native -Dnative.mode=debug-O
CXX=g++-13 ant native -Dnative.mode=debug-O
# Equivalent direct Make commands:
make -C jni/quasse MODE=debug-O clean
make -C jni/quasse MODE=debug-O CXX=g++-13
```

Cleaning removes only the selected directory's library, generated JNI headers/classes, and any
obsolete configuration records. It preserves unrelated files and other builds. Make's `clean`
requires no Java, compiler, or FFTW installation; `ant clean-native` needs only Ant's Java runtime.
Use the same directory override for cleaning a custom build.
Checkout paths may contain spaces when using the relative output paths above; arbitrary output
directory names containing spaces are not supported by the Makefile.
