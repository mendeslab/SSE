## Generate the diversitree reference value used by
## QuaSSEDistributionTest.testTraitDependentExtinctionLikelihoodAgainstDiversitree.
## Nonzero-drift references require the locally corrected diversitree, not the ordinary installation.
## Reference source: ../diversitree, commit 5af1bf7bef1742512f3bc84beaf0e6115fc08e42, version 0.10-1.
## See validation/QuaSSEReference.md for the isolated build and exact invocation.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L)
  stop("Usage: Rscript validation/r_scripts/QuaSSENonconstantMuReference.R ISOLATED_LIBRARY")
reference.library <- normalizePath(args[1], mustWork=TRUE)
.libPaths(c(reference.library, .libPaths()))
library(diversitree, lib.loc=reference.library)

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
cat("reference source commit: 5af1bf7bef1742512f3bc84beaf0e6115fc08e42\n")
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
