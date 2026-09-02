# Code base notes

This document is meant to guide developers through the code base.
It will:

* Explain the purposes of the different classes, giving some explanation on
relevant implementation details
* Go over the unit tests in a logical order

## Configuring JNI for MoSSE C libraries
Building requires FFTW3, a C compiler, and a JDK compatible with the Java used
to run BEAST. SSE targets Java 17 and also builds with newer JDKs.

Install fftw3 on ubuntu  
```
sudo apt-get install -y libfftw3-dev
```

Install Java 17
```
sudo apt-get install openjdk-17-jdk
```

Build the JNI methods from the repository root. The Makefile uses `JAVA_HOME`
when it is set and otherwise derives the JDK location from `javac`:

```
make -C jni
```

Run the portable tests and the native MoSSE tests separately:

```
ant test
ant test-mosse
```

After building the JNI library, `ant test-all` runs both suites. `test-mosse`
and `test-all` supply the JNI library path. The native tests cover
`MosseTreeLikelihood` initialization, but its likelihood calculation remains
disabled because that unfinished path can abort inside the native integration
code.

The native build currently defaults to Linux. Set `JNI_PLATFORM` to the name of
the appropriate JDK include subdirectory when porting it to another platform.

## QuaSSE classes

QuaSSE currently performs its active FFTs in Java using `JavaFftService` from
the bundled shared library. A JTransforms implementation is retained and
covered by lower-level tests, but is not selected by the likelihood. The JNI
setup above applies only to the separate MoSSE implementation. The retained
QuaSSE C wrapper is not connected to the Java implementation.

The QuaSSE computational grid is fixed when the likelihood is initialized.
Its padding and trait rulers are derived from the initial constant drift,
diffusion, and `dtMax`. Fixed nonzero drift is supported, but changing drift or
diffusion during an MCMC run can require different padding and is not yet
supported safely. Dynamic grid reconstruction and its cache invalidation must
be implemented before those parameters are inferred.

Use `dynDt=true`. This mode divides each branch into an integer number of equal
steps no longer than `dtMax`. The retained `dynDt=false` path does not integrate
a branch-length remainder and is unfinished. FFTW integration, incomplete
sampling, and other unfinished QuaSSE extensions are also deferred.

Each `calculateLogP()` currently performs a full pruning calculation from
freshly populated tip values. Incremental BEAST node-partial caching is a later
performance feature rather than part of the current correctness model.

### QuaSSEProcess

TODO

### QuaSSEDistribution

TODO

### MoSSEDistribution

TODO

## QuaSSE unit tests

### QuaSSEFunctionsTest

See `validation/r_scripts/QuaSSEFunctionsTest_JUnitTest.R` for where
expectations come from.

* `testLogistic`: tests the LogisticFunction class, which converts an
  x (continuous trait) value into a y (macroevolutionary parameter, e.g.,
  the birth rate) value

* `testConstant`: tests the ConstantLinkFn class, which converts an
  x (continuous trait) value and maps to a constant scalar value for y
  (a macroevolutionary parameter, e.g., the death rate)

* `testQu2MacroevolFailLogistic`: tests LogisticFunction throws an
  exception if the number of bins for x and y differ

* `testQu2MacroevolFailConstant`: tests ConstantLinkFn throws an
  exception if the number of bins for x and y differ

### PropagatesQuaSSETest

See `validation/r_script/PropagatesQuaSSETest_JUnitTest.R` for where
expectations come from

None of the tests here depend on, nor test, the initialization of
likelihood classes (i.e., we do not test if arrays are declared with
the right dimensions and correctly populated with initial values).

* `testPropagateTimeOneChQuaSSETest`: tests
  SSEUtils.propagateEandDinTQuaSSEinPlace.
  
  It uses a manually creates array of random values as input, with a
  total of 32 bins. It checks that both E's and D's are propagated in
  time correctly.

* `testPropagateTimeOneChQuaSSETestSSTJavaFftService`: tests
  SSEUtils.propagateEandDinTQuasseInPlaceSSTJavaFftService.
  
  It uses a manually created array of random values, spread over an
  array while skipping even indices, as required by SST's
  JavaFftService (which is tested elsewhere!).
  
  (Note that we do not use SST in this test, this is just to make
  sure this method skips indices correctly.)

* `testMakeNormalKernInPlaceAndFftAndIfft`: tests both
  SSEUtils.makeNormalKernelInPlace and FFTing with JTransforms.
  
  The Normal kernel assigns a different probability to each of
  'nXbins' bins. Note that x here is not a value a continuous trait
  can have, but instead, a certain delta-trait-value. So most of the
  probability mass of this kernel is centered around zero (we are
  assuming the continuous trait evolves according to a Brownian motion
  model).
  
  Now for some technical reason, we put the bins that would have the
  highest probabilities in the first 'nRightFlankBins' to the left,
  and the last 'nLeftFlankBins' to the right.
  The middle x bins are assigned probabilities of zero.
  
  This test assumes 32 bins and 2 flanking bins to the left, and to
  the right.

* `testMakeNormalKernInPlaceAndFftAndIfftSSTModalFftService`: tests
  both SSEUtils.makeNormalKernelInPlace and FFTing with the SST
  library, using RealArray methods (which invoke ModalFfftService).
  
  (See test above for more details.)

* `testConvolve`:

* `testConvolveSSTModalFftService`:

* `testConvolveSSTJavaFftService`:

* `testPropagateChOneCh32QuaSSE`:

* `testPropagateChOneCh1024QuaSSE`:

* `testPropagateChOneCh4096QuaSSE`:

* `testPropagateChOneCh4096QuaSSEModalFftService`:

* `testPropagateChOneCh4096QuaSSEJavaFftService`:

### QuaSSEDistributionTest

TODO
