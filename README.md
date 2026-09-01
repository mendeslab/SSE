# SSE


SSE is a [BEAST2](http://beast2.org) package for analyses employing SSE-type models (State-dependent Speciation and Extinction models) in a Bayesian framework.
The SSE package currently implements a well calibrated ClaSSE model (Goldberg and Igic, 2012), with support for its submodels BiSSE, MuSSE and GeoSSE.
We are finalizing the implementation of a stochastic character mapper as in Freyman and Hohna (2018).

## Command-line build with BEAST 2.7

Build the sibling `beast2` and `BeastFX` repositories first. Then compile SSE
classes directly into `build` without creating an SSE jar:

```sh
ant compile
```

If the BEAST source tree is elsewhere, pass it explicitly:

```sh
ant -Dbeast.dir=/path/to/beast2 compile
```

Run an example with the compiled class directory on the Java class path:

```sh
java -cp 'build:../beast2/build/dist/launcher.jar:'\
'../beast2/build/dist/BEAST.base.jar:'\
'../BeastFX/build/dist/BEAST.app.jar:lib/*' \
  beastfx.app.beast.BeastMCMC -overwrite examples/BiSSE_fixed_tree_SDSEP.xml
```

Do not add the historical `SSE.jar` to this class path: it bundles BEAST 2.6
classes that conflict with BEAST 2.7. The separate `src/mosse` source tree is
retained but is not part of this compile target.
