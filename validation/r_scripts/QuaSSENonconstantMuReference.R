## Generate the diversitree reference value used by
## QuaSSEDistributionTest.testTraitDependentExtinctionLikelihoodAgainstDiversitree.
## Nonzero-drift references require the locally corrected diversitree, not the ordinary installation.
## Expected source: ../diversitree, commit 5af1bf7bef1742512f3bc84beaf0e6115fc08e42, version 0.10-1.
## See validation/QuaSSEReference.md for the isolated build and exact invocation.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L)
  stop("Usage: Rscript validation/r_scripts/QuaSSENonconstantMuReference.R ISOLATED_LIBRARY")
reference.library <- normalizePath(args[1], mustWork=TRUE)
.libPaths(c(reference.library, .libPaths()))
library(diversitree, lib.loc=reference.library)

## Sum input[i-offset] for every circular kernel offset, independently of either FFT backend.
## Entry j represents offset j-1 modulo n; wrapping here fixes the convolution orientation.
quasse.direct.convolution <- function(x, kernel) {
  n <- length(x)
  vapply(seq_len(n), function(i)
    sum(kernel * x[((i - seq_len(n)) %% n) + 1L]), numeric(1))
}

## Reject uncorrected references before generating values: backend agreement with negligible
## boundary D misses these defects. Adapted from test-quasse-internal.R; retain while references
## require the local padding/boundary corrections rather than an upstream implementation.
support.control <- list(nx=1024, dx=.1, dt.max=1, xmid=0, r=4L, w=3, method="fftC")
for (v in c(2, -2)) {
  ## For v=2 the left extent peaks at t=1: 2+3=5; the right peaks at t=(3/4)^2:
  ## -2*t+3*sqrt(t)=1.125. Divide by dx and round up; negative drift swaps sides.
  expected.padding <- if (v > 0) c(50, 12) else c(12, 50)
  extent <- diversitree:::quasse.extent(support.control, v, 1)
  if (!isTRUE(all.equal(unname(extent$padding[2,]), expected.padding)) ||
      !isTRUE(all.equal(unname(extent$padding[1,]), 4 * expected.padding)))
    stop("Reference preflight: diversitree lacks corrected asymmetric/all-timestep padding.")
}

nx <- 32L
dx <- .1
ptr <- .Call(diversitree:::r_make_quasse_fft, nx, dx, 2L, -1L)
for (padding in list(c(5L, 2L), c(2L, 5L))) {
  nkl <- padding[1]
  nkr <- padding[2]
  ndat <- nx - nkl - nkr - 1L
  drift <- if (nkl > nkr) .2 else -.2
  offsets <- c(0:nkr, rep(NA, nx-nkl-nkr-1L), -nkl:-1)
  kernel <- dnorm(offsets * dx, -drift * .1, sqrt(.01 * .1))
  kernel[is.na(kernel)] <- 0
  kernel <- kernel / sum(kernel)
  vars <- cbind(c(seq_len(ndat) / ndat, rep(0, nx-ndat)),
                c(rev(seq_len(ndat)) / ndat, rep(0, nx-ndat)))

  ## The C entry first propagates time. For pure birth, E'=E/h and D'=z*D/h^2,
  ## with z=exp(lambda*dt), h=z*(1-E)+E. Convolve that result, zero unused bins,
  ## then restore the dependency strips: nkr on the left, nkl on the right.
  z <- exp(.1 * .1)
  h <- z * (1-vars[,1]) + vars[,1]
  vars.t <- cbind(vars[,1] / h, z * vars[,2] / h^2)
  expected <- apply(vars.t, 2, quasse.direct.convolution, kernel=kernel)
  boundary <- c(seq_len(nkr), ndat-seq_len(nkl)+1L)
  expected[(ndat+1L):nx,] <- 0
  expected[boundary,1] <- vars.t[boundary,1]
  actual.c <- .Call(diversitree:::r_do_integrate, ptr, vars, rep(.1, ndat),
                    rep(0, ndat), drift, .01, 1L, .1, padding)
  if (!all(is.finite(actual.c)) || max(abs(actual.c-expected)) > 1e-12)
    stop("Reference preflight: fftC lacks corrected asymmetric E boundary restoration.")

  ## Preserve the existing backend distinction: fftR restores both E and D, fftC only E.
  expected[boundary,2] <- vars.t[boundary,2]
  actual.r <- diversitree:::fftR.propagate.x(vars.t, nx, fft(kernel), nkl, nkr)
  if (!all(is.finite(actual.r)) || max(abs(actual.r-expected)) > 1e-12)
    stop("Reference preflight: fftR lacks corrected asymmetric E/D boundary restoration.")
}

tree <- ape::read.tree(text="(sp1:0.1,sp2:0.1);")
states <- c(sp1=0.0, sp2=0.1)
states.sd <- 0.05

control <- list(tc=0.005,
                dt.max=0.005,
                nx=128,
                dx=0.01,
                r=4L,
                xmid=0.0,
                w=10)

## Parameters are constant lambda, logistic-mu base, maximum, midpoint, and growth rate,
## followed by drift and diffusion.
parameters <- c(0.15, 0.01, 0.08, 0.05, 20.0, 0.0, 0.001)

control.fftR <- c(control, list(method="fftR"))
likelihood.fftR <- make.quasse(tree, states, states.sd, constant.x, sigmoid.x, control.fftR)
log.likelihood.fftR <- likelihood.fftR(parameters)

control.fftC <- c(control, list(method="fftC"))
likelihood.fftC <- make.quasse(tree, states, states.sd, constant.x, sigmoid.x, control.fftC)
log.likelihood.fftC <- likelihood.fftC(parameters)

## Agreement between independent diversitree FFT paths makes the scalar a stronger reference.
stopifnot(abs(log.likelihood.fftR - log.likelihood.fftC) < 1e-12)

## The former Java bug evaluated every mu bin at the initial zero-filled rate array. Confirm that
## this test would distinguish that constant mu(0) behavior from trait-dependent extinction.
mu.at.zero <- sigmoid.x(0.0, 0.01, 0.08, 0.05, 20.0)
constant.mu.likelihood <- make.quasse(tree, states, states.sd, constant.x, constant.x, control.fftR)
constant.mu.parameters <- c(0.15, mu.at.zero, 0.0, 0.001)
log.likelihood.constant.mu <- constant.mu.likelihood(constant.mu.parameters)
stopifnot(abs(log.likelihood.fftR - log.likelihood.constant.mu) > 1e-6)

cat("diversitree version:", as.character(packageVersion("diversitree")), "\n")
cat("parameters:", paste(format(parameters, digits=17), collapse=", "), "\n")
cat("fftR log likelihood:", format(log.likelihood.fftR, digits=17), "\n")
cat("fftC log likelihood:", format(log.likelihood.fftC, digits=17), "\n")
cat("constant mu(0) log likelihood:", format(log.likelihood.constant.mu, digits=17), "\n")

## Both backends must agree, and widening the domain without changing dx must not change the answer.
## This guards against legitimising a shared boundary artefact by comparing only Java with one backend.
cases <- rbind(c(0, .001), c(0, .004), c(-.1, .001), c(.1, .001), c(-1, .001), c(1, .001))
cat("reference library:", find.package("diversitree"), "\n")
cat("expected reference source commit: 5af1bf7bef1742512f3bc84beaf0e6115fc08e42\n")
cat("reference preflight: passed; installed source revision is not verified\n")
cat("drift diffusion left right xLoMin xHiMin logP\n")
for (i in seq_len(nrow(cases))) {
  parameters[6:7] <- cases[i,]
  c.value <- likelihood.fftC(parameters)
  r.value <- likelihood.fftR(parameters)
  stopifnot(abs(c.value-r.value) < 1e-12)
  wider <- control.fftC
  wider$nx <- 256L
  wide.c <- make.quasse(tree, states, states.sd, constant.x, sigmoid.x, wider)(parameters)
  wider$method <- "fftR"
  wide.r <- make.quasse(tree, states, states.sd, constant.x, sigmoid.x, wider)(parameters)
  stopifnot(abs(wide.c-wide.r) < 1e-12, abs(c.value-wide.c) < 1e-10)
  extent <- diversitree:::quasse.extent(control.fftC, parameters[6], parameters[7])
  stopifnot(isTRUE(all.equal(extent$x[[1]][extent$tr], extent$x[[2]])))
  cat(sprintf("%.17g %.17g %d %d %.17g %.17g %.17g\n",
              parameters[6], parameters[7], extent$padding[2,1], extent$padding[2,2],
              extent$x[[2]][1], extent$x[[1]][1], c.value))
}
