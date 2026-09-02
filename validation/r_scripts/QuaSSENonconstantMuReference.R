## Generate the diversitree reference value used by
## QuaSSEDistributionTest.testTraitDependentExtinctionLikelihoodAgainstDiversitree.
## Drift is zero so that the suspected asymmetric-padding error in diversitree is irrelevant.

library(diversitree)

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
