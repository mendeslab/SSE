# Generating input, running and parsing example .xml's    

This will guide you through the preparation of input for the .xml's in the examples/ folder.
You are going to need the following R packages:

* ape
* diversitree
* hisse
* ggplot2
* gridExtra
* RColorBrewer

As you complete step (1), you will produce the .xml's named in the headers and should run BEAST (the SSE package) on them.
In the end, you will have the following files:

* examples/BiSSE_fixed_tree_SDSEP.log
* examples/ClaSSE_fixed_tree_SDSEP.log
* examples/BiSSE_fixed_tree_on_HiSSE_HSDSEP.log
* examples/HiSSE_fixed_tree_on_HiSSE_HSDSEP.log
* examples/ModelAveraging_fixed_tree_on_HiSSE_BSSVSSDSEP.log
* validation/BiSSE_fixed_tree_SDSEP_SCM_parsed.txt

You can skip step (1) and go straight to plotting the posteriors and the ancestral state reconstructions in step (2).
In this case, you will need to use the provided files (inside the examples/ and validation/ folders).

## (1) Running simulations in R (see examples_xml_inputs.R)    
### (1.1) BiSSE_fixed_tree_SDSEP.xml    
Simulating and computing MLE under BiSSE.    

```
library(ape)
library(diversitree)
library(hisse)

pars <- c(.1, .5, .05, .05, .1, .1) # lambdas, mus, qs

set.seed(12345)
phy <- tree.bisse(pars, max.taxa=120, include.extinct=FALSE, x0=NA)

# tree (for .xml)
write.tree(phy) # see tree

# taxa (for .xml)
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

# tip data (for .xml)
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

# find MLE
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)
```

### (1.2) ClaSSE_fixed_tree_SDSEP.xml
Simulating under ClaSSE (computing MLE with diversitree bombs).    

```
argnames(make.classe(tree=phy, states=phy$tip.state+1, k=3, strict=FALSE)) # just to show order of parameters in classe object

pars <- c(.1, 0.0, 0.0, 0.0, 0.0, 0.0,
          0.0, 0.0, 0.0, .1, 0.0, 0.0,
          0.0, .5, .3, 0.0, .3, .1,
          .05, .05, .05, # mu1, mu2, mu3
          .1, .1, .1, # q12, q13, q21
          .1, .1, .1) # q23, q31, q32

set.seed(12345)
phy <- tree.classe(pars, max.taxa=120, x0=NA, include.extinct=FALSE)

# tree
write.tree(phy)

# taxa
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

# tip data
paste(paste(phy$tip.label, (phy$tip.state), sep="="), collapse=",")

# finding mle
lik <- make.classe(phy, phy$tip.state, k=3)
fit <- find.mle(lik, pars) # bombs!
```

### (1.3) HiSSE_fixed_tree_on_HiSSE_HSDSEP.xml and BiSSE_fixed_tree_on_HiSSE_HSDSEP.xml    
Simulating with hidden trait (one hidden state; 0 ⇔ 1 ⇔ 1H), and computing MLEs.    

```
pars <- c(.1,  .1,  .5,  # lambda 0, 1, 2
.05, .05, .05, # mu 0, 1, 2
.1, 0.0, # q01, q02
.1, .1, # q10, q12
0.0, .1 # q20, q21
) # pars above are equivalent to Fig. 1 in HiSSE paper

set.seed(10000)
phy <- tree.musse(pars, max.taxa=120, include.extinct=FALSE, x0=1)
phy$tip.state[phy$tip.state==3] <- 2 # hiding states
phy$tip.state <- phy$tip.state - 1

# tree (for .xml)
write.tree(phy) # see tree

# taxa (for .xml)
cat(paste0("<taxon id=\"", phy$tip.label, "\" spec=\"Taxon\"/>"), sep="\n")

# tip data (for .xml)
paste(paste(phy$tip.label, (phy$tip.state + 1), sep="="), collapse=",")

# finding MLE under BiSSE
pars.to.estimate <- pars <- c(.1, .3, .05, .05, .1, .1) # lambdas, mus, qs
lik <- make.bisse(phy, phy$tip.state)
fit <- find.mle(lik, pars)

# finding MLE under HiSSE
turnover.anc <- c(1,2,0,3)
eps.anc <- c(1,2,0,3)

trans.rates <- TransMatMaker(hidden.states=TRUE)
trans.rates.nodual.no0B <- ParDrop(trans.rates,c(2,3,5,7,8,9,10,12))
trans.rates.nodual.no0B

sim.dat <- data.frame(names(phy$tip.state), phy$tip.state)

pp <- hisse(phy, sim.dat, f=c(1,1), hidden.states=TRUE, turnover.anc=turnover.anc, eps.anc=eps.anc, trans.rate=trans.rates.nodual.no0B, output.type="raw", root.type="equal", condition.on.survival=FALSE, root.p=NULL)
```

### (1.4) ModelAveraging_fixed_tree_on_HiSSE_BSSVSSDSEP.xml
Same simulation of (1.2).
After running .xml, remove header lines (starting with '#') from .log file, and save new file as 'ModelAveraging_fixed_tree_on_HiSSE_BSSVSSDSEP_noheader.log'.

### (1.5) BiSSE_fixed_tree_SDSEP_SCM.trees (stochastic character mapping on 60-sp tree under BiSSE)

The historical `scripts/parse_asm_treesfile.py` helper is not present in this
repository or its available history. This post-processing step cannot currently
be reproduced without obtaining or replacing that parser.

## (2) Plotting all graphs in R (see examples_xml_plots.R)    

You can find all the code to plot the posterior distributions and other graphs in examples_xml_plots.R.
Below, "/path/to/SSE/validation" should be replaced with your local path to the validation folder that comes with this git repository.
We will just call the script and let it do all the work.    

```
mkdir /path/to/SSE/validation/plots
Rscript r_scripts/example_xml_plots.R /path/to/SSE/validation
```
