## Reproduce the rounded 15-species Java test with the fixes + Strang diversitree fork.
## https://github.com/bredelings/diversitree, commit 7efd65f73068cb32fef3f5966a0e6b297fe2c893.
## Usage: Rscript validation/r_scripts/QuaSSEStrangReference.R ISOLATED_LIBRARY
## Build the corrected package as documented in validation/QuaSSEReference.md.
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 1L) stop("Supply the isolated corrected diversitree library")
library(diversitree, lib.loc=normalizePath(args[1], mustWork=TRUE))
tr <- ape::read.tree(text=paste0(
  "(sp2:13.77320255,(sp1:12.76688384,((((sp12:1.170387028,sp13:1.170387028)nd16:0.9837720325,",
  "sp9:2.154159061)nd11:5.451401092,((sp5:4.311645343,(sp14:0.8910055279,sp15:0.8910055279)nd",
  "14:3.420639815)nd9:2.536663776,((sp16:0.3011866125,sp17:0.3011866125)nd12:4.264383667,(sp6",
  ":3.95083843,sp7:3.95083843)nd13:0.6147318498)nd10:2.282738839)nd8:0.7572510339)nd5:2.55473",
  "9141,((sp10:2.059478202,sp11:2.059478202)nd15:0.4198789018,sp8:2.479357104)nd6:7.68094219)",
  "nd4:2.60658455)nd3:1.006318707)nd1;"))
traits <- setNames(c(-.05384594,-.37091896,.59169195,.14947513,.46156791,
                    .27345680,1.73358959,.38883347,.42233625,1.55011787,
                    1.17169681,.72422971,.84092251,.19645523,.48495092),
                  c("sp1","sp2","sp5","sp6","sp7","sp8","sp9","sp10","sp11",
                    "sp12","sp13","sp14","sp15","sp16","sp17"))
pars <- c(.1,.2,0,2.5,.03,0,.01)
for (method in c("fftR", "fftC")) {
  lik <- make.quasse(tr,traits,.02,sigmoid.x,constant.x,
                    list(nx=1024L,dx=.01027592,xmid=.6813353,tc=1.37732,
                         dt.max=.005,r=4L,w=5,method=method))
  cat(method, "Strang log likelihood:",format(lik(pars),digits=17),"\n")
}

## Match TreeParser threshold=1e-5: extend terminal edges to the latest tip without moving
## internal nodes. Keep the unadjusted result above as the historical rounded-tree reference.
depths <- ape::node.depth.edgelength(tr)
ages <- max(depths[seq_len(length(tr$tip.label))]) - depths[seq_len(length(tr$tip.label))]
stopifnot(all(ages >= 0), all(ages < 1e-5))
terminal <- which(tr$edge[,2] <= length(tr$tip.label))
tr$edge.length[terminal] <- tr$edge.length[terminal] + ages[tr$edge[terminal,2]]
for (method in c("fftR", "fftC")) {
  lik <- make.quasse(tr,traits,.02,sigmoid.x,constant.x,
                    list(nx=1024L,dx=.01027592,xmid=.6813353,tc=1.37732,
                         dt.max=.005,r=4L,w=5,method=method))
  cat(method, "Strang log likelihood, threshold-adjusted tips:",format(lik(pars),digits=17),"\n")
}
