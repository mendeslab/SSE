### author: Fabio K. Mendes
## Historical inline likelihood output below predates the root-conditioning correction.
## Conditioned log likelihoods increase by -2*log(root-grid dx); pre-root arrays are unchanged.
## Use the isolated corrected package documented in ../QuaSSEReference.md for current references.
##
## This R script gives us the expected values for JUnit tests
## inside QuaSSEDistributionTest.java
##
## (1) testDimensions
## (2) testInitializationOfTips
## (3) testIntegrateOneBranchLoRes32BinsOutsideClassJustT
## (4) testIntegrateOneBranchLoRes32BinsOutsideClassJustX
## (5) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX
## (6) testIntegrateOneBranchHiRes4096BinsOutsideClassJustX (note that I'm using this in 4096 low res, but I'll co-opt 4096 high res in Java)
## (7) testBothXandTPropagateMethodsInsideClassOneBranchLoRes32Bins
## (8) testBothXandTPropagateMethodsInsideClassOneBranchLoRes32BinsDt002
## (9) testIntegrateOneBranchHiRes32BinsInsideClassBothXandTDt002
## (10) testPriorProbAtRootObserved
## (11) testPriorProbAtRootFlat
## (12) testRootCalcProcedure
## (13) testPruneBifTree32Bins
## (14) testPruneBifTree0025Height32Bins
## (15) testPruneTrifTree32Bins

## (20) testPruneFifteenSpTree1024Bins

## IMPORTANT: expectations (outputs from paste(xyz, collapse=", ") will vary depending on machine architecture; when dnorm() gets very small inputs (e.g., dx * nx are both large), the initial values of D will be tiny, and then the result of FFT on different machines will differ
##
## This is the reason for 48 bins having a dx = 0.01, and 1024/4096 bins having a dx = 0.0005


library(diversitree)

## First, diversitree's functions
source("rev_eng_diversitree.R")

## Imports that are generally hidden.
## make.pde.quasse.fftC <- diversitree:::make.pde.quasse.fftC
## make.pde.quasse.fftR <- diversitree:::make.pde.quasse.fftR
## make.pde.quasse.mol <- diversitree:::make.pde.quasse.mol

## make.branches.quasse.fftC <- diversitree:::make.branches.quasse.fftC
## make.branches.quasse.mol <- diversitree:::make.branches.quasse.mol

## make.cache.quasse <- diversitree:::make.cache.quasse

## Test expectations
##
## (1) testDimensions

dt <- 0.01
dx <- 0.0005
nx <- 1024
hi.lo.ratio <- 4
w <- 10
drift <- 0.0
diffusion <- 0.001
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y

mean4test <- drift * dt
sd4test <- sqrt(diffusion * dt);

nkleft <- max(ceiling(-(mean4test - w * sd4test)/dx)) * c(hi.lo.ratio, 1)
nkright <- max(ceiling((mean4test + w * sd4test)/dx)) * c(hi.lo.ratio, 1)

ndat <- nx * c(hi.lo.ratio, 1) - (nkleft + 1 + nkright)
ndat.lo <- ndat[2]; ndat.lo # 895
ndat.hi <- ndat[1]; ndat.hi # 3583

xmin.lo <- xmid - dx * ceiling((ndat.lo - 1)/2) # x.02
xmin.hi <- xmin.lo - dx * (1 - 1/hi.lo.ratio) # x.01
x.lo <- seq(xmin.lo, length.out=ndat.lo, by = dx) # same as ext.fft$x[[2]]
x.hi <- seq(xmin.hi, length.out=ndat.hi, by = dx/hi.lo.ratio) # same as ext.fft$x[[1]]

paste(x.lo[1:10], collapse=", ")
## expectedXLoFirst10 = -0.2235, -0.223, -0.2225, -0.222, -0.2215, -0.221, -0.2205, -0.22, -0.2195, -0.219

paste(x.lo[886:895], collapse=", ")
## expectedXLoLast10 = 0.219, 0.2195, 0.22, 0.2205, 0.221, 0.2215, 0.222, 0.2225, 0.223, 0.2235

paste(x.hi[1:10], collapse=", ")
## expectedXHiFirst10 = -0.223875, -0.22375, -0.223625, -0.2235, -0.223375, -0.22325, -0.223125, -0.223, -0.222875, -0.22275

paste(x.hi[3574:3583], collapse=", ")
## expectedXHiFirst10 = 0.22275, 0.222875, 0.223, 0.223125, 0.22325, 0.223375, 0.2235, 0.223625, 0.22375, 0.223875

lambda <- sigmoid.x # lambda is now a logistic function
mu <- constant.x # mu is now a constant function

# lambdas
paste(do.call(lambda, c(list(x.lo), c(y0, y1, xmid, r)))[1:10], collapse=", ")
## expectedLambdaLoFirt10 = 0.13638367349015, 0.136412610857295, 0.13644155805569, 0.136470515067719, 0.136499481875741, 0.136528458462084, 0.136557444809052, 0.13658644089892, 0.136615446713936, 0.136644462236322

paste(do.call(lambda, c(list(x.lo), c(y0, y1, xmid, r)))[886:895], collapse=", ")
## expectedLambdaLoLast10 = 0.163355537763678, 0.163384553286064, 0.16341355910108, 0.163442555190948, 0.163471541537916, 0.163500518124259, 0.163529484932281, 0.16355844194431, 0.163587389142705, 0.16361632650985

paste(do.call(lambda, c(list(x.hi), c(y0, y1, xmid, r)))[1:10], collapse=", ")
## expectedLambdaHiFirt10 = 0.136361976927131, 0.136369208498794, 0.136376440686558, 0.13638367349015, 0.136390906909294, 0.136398140943716, 0.136405375593141, 0.136412610857295, 0.136419846735901, 0.136427083228686

paste(do.call(lambda, c(list(x.hi), c(y0, y1, xmid, r)))[3574:3583], collapse=", ")
## expectedLambdaHiLast10 = 0.163572916771314, 0.163580153264099, 0.163587389142705, 0.163594624406859, 0.163601859056284, 0.163609093090706, 0.16361632650985, 0.163623559313442, 0.163630791501206, 0.163638023072869

# mus
paste(do.call(mu, c(list(x.lo), death))[1:10], collapse=", ") # expectedMuLoFirst10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(do.call(mu, c(list(x.lo), death))[886:895], collapse=", ") # expectedMuLoLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03

paste(do.call(mu, c(list(x.hi), death))[1:10], collapse=", ") # expectedMuHiFirst10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(do.call(mu, c(list(x.hi), death))[3574:3583], collapse=", ") # expectedMuHiLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03



# (2) testInitializationOfTips

dx <- 0.0005
hi.lo.ratio <- 4
w <- 10
drift <- 0.0
diffusion <- 0.001
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y
lambda <- function(x) sigmoid.x(x, y0, y1, xmid, r)
mu <- function(x) constant.x(x, death)

char <- make.brownian.with.drift(drift, diffusion)

set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=2, x0=0, single.lineage=FALSE, verbose=FALSE)
pars <- c(y0, y1, xmid, r, death, drift, diffusion) # 6th and 7th elements are drift and diffusion
sd <- 1/20
control.C.1 <- list(dt.max=0.01) # dt = 0.01
phy$tip.state[1] <- 0.0
phy$tip.state[2] <- 0.1

cache <- diversitree:::make.cache.quasse(phy, phy$tip.state, sd, lambda, mu, control.C.1, NULL)
cache$height[1] <- cache$height[2] <- 0.01 # setting branch length (node height) manually
cache$depth[3] <- 0.01 # setting depth of ancestral node "nd1"
cache$edge.length <- c(0.01, 0.01)
cache$len <- c(0.01, 0.01, NA)

## these four lines are to make this cache match the previous unit tests
cache$args$drift <- 6
cache$args$diffusion <- 7
cache$control$dx <- dx ## to match unit test
cache$control$xmid <- xmid ## to match unit test
cache$control$w <- w

## From now on: taking following lines from R/model-quasse.R, make.quasse()

## all.branches <- make.all.branches.quasse(cache, cache$control) # all.branches is a function
f.pars <- diversitree:::make.pars.quasse(cache) # f.pars is a function
pars2 <- f.pars(pars)

# understanding tip initialization
sampling.f <- 1
e0 <- 1 - sampling.f
nx <- cache$control$nx * cache$control$r # 1024 * 4 = 4096
npad <- nx - length(pars2[[1]]$x)

# this is how the tip y's are initialized
sp1.y <- c(rep(e0, nx), dnorm(pars2[[1]]$x, cache$states[1], cache$states.sd[1]), rep(0, npad)) # E's and D's (8192 elements)
sp2.y <- c(rep(e0, nx), dnorm(pars2[[1]]$x, cache$states[2], cache$states.sd[2]), rep(0, npad))

# we can compare it to running the function as is
cache$y <- diversitree:::initial.tip.quasse(cache, cache$control, pars2[[1]]$x) # this is a list

# cache$y$t # branch lengths
# cache$y$target # tip indices (I think)
identical(cache$y$y$sp1, sp1.y) # TRUE
identical(cache$y$y$sp2, sp2.y) # TRUE
identical(cache$y$y$sp1, sp2.y) # FALSE (just as a control)

sp1.y.ds <- sp1.y[nx+1:(length(sp1.y)-nx)] # skipping the E's, grabbing second half (D's)
sp2.y.ds <- sp2.y[nx+1:(length(sp2.y)-nx)] # skipping the E's, grabbing second half (D's)

sp1.y.ds.exp <- sp1.y.ds[2001:2048] # 48 arbitrarily consecutive elements for checking against (if using JTransforms, which takes consecutive doubles at the head of array)
sp2.y.ds.exp <- sp2.y.ds[2001:2048] # 48 arbitrarily consecutive elements for checking against  (if using JTransforms, which takes consecutive doubles at the head of array)

paste(sp1.y.ds.exp, collapse=", ")
## (if using JTransforms) expectedSp1Ds = 6.96077358436849, 6.95166528584696, 6.94252551477923, 6.93335442671583, 6.92415217758908, 6.91491892370871, 6.90565482175746, 6.89636002878667, 6.88703470221182, 6.87767899980818, 6.86829307970631, 6.85887710038768, 6.84943122068018, 6.83995559975374, 6.83045039711584, 6.82091577260705, 6.81135188639661, 6.80175889897795, 6.79213697116422, 6.78248626408384, 6.772806939176, 6.76309915818623, 6.75336308316186, 6.74359887644761, 6.73380670068105, 6.72398671878815, 6.71413909397875, 6.70426398974212, 6.69436156984244, 6.68443199831428, 6.67447543945815, 6.66449205783599, 6.65448201826663, 6.64444548582134, 6.63438262581928, 6.62429360382306, 6.61417858563415, 6.60403773728847, 6.59387122505179, 6.58367921541529, 6.57346187509105, 6.5632193710075, 6.55295187030495, 6.54265954033109, 6.53234254863645, 6.52200106296994, 6.51163525127429, 6.50124528168164

paste(sp2.y.ds.exp, collapse=", ")
## (if using JTransforms) expectedSp2Ds = 2.67859021074856, 2.68849414736158, 2.69841783805887, 2.70836123148143, 2.71832427571076, 2.72830691826807, 2.73830910611356, 2.74833078564564, 2.75837190270027, 2.76843240255022, 2.77851222990441, 2.78861132890721, 2.79872964313779, 2.80886711560949, 2.81902368876922, 2.82919930449678, 2.83939390410431, 2.84960742833573, 2.85983981736613, 2.87009101080125, 2.88036094767694, 2.89064956645866, 2.90095680504094, 2.91128260074695, 2.92162689032799, 2.93198960996305, 2.9423706952584, 2.95277008124712, 2.96318770238874, 2.97362349256884, 2.9840773850987, 2.9945493127149, 3.00503920757903, 3.01554700127736, 3.02607262482052, 3.03661600864323, 3.04717708260405, 3.05775577598507, 3.06835201749176, 3.07896573525268, 3.0895968568193, 3.10024530916586, 3.11091101868917, 3.12159391120842, 3.13229391196515, 3.14301094562307, 3.15374493626797, 3.16449580740766



# (3) testIntegrateOneBranchLoRes32BinsOutsideClassJustT

dt <- 0.01
dx <- 0.01
w <- 10
drift <- 0.0
diffusion <- 0.001
sd <- 0.05
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y
lambda <- sigmoid.x # lambda is now a function
mu <- constant.x # mu is now a function

control.fft.32 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=dt, # dt
                    nx=32, # number of X bins
                    dx=dx, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=w, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

args.fft.just.t <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.t <- c(y0, y1, xmid, r, death, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.t <- quasse.extent.debug(control.fft.32, drift, diffusion) # prepares X axis stuff
## ext.fft.just.t$padding

pars.fft.just.t <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.just.t, ext.fft.just.t, pars.just.t) # adds lambda and mu vectors to ext.fft.just.x

## initialization
vars.fft.just.t <- matrix(0, control.fft.32$nx, 2) # low resolution

## populating states in all bins, at tips, being drawn from a normal distribution centered at a species' observed value
## which here is assumed to be 0.0
##
## also note that the number of left and right flanking bins is 4
##
## pars.fft.just.t$lo$padding
## nkl nkr
##   4   4
##
## and if you look at the last (4 + 4) rows
## in vars.fft.just.r, should be 0.0
vars.fft.just.t[seq_len(ext.fft.just.t$ndat[2]),2] <- dnorm(ext.fft.just.t$x[[2]], 0.0, sd)
## states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x$lo$padding should be 0.0

## alternating elements with 0's to match SST's JavaFftService requirements
exp.sp1.initial.ds <- rep(0, length(vars.fft.just.t[,2])*2)
idxs <- seq(1, length(vars.fft.just.t[,2])*2, by=2)
j <- 1
for (i in vars.fft.just.t[,2]) {
    exp.sp1.initial.ds[idxs[j]] = i
    j = j+1
}

paste(exp.sp1.initial.ds, collapse=", ")
## expectedInitialDs = 0.709491856924629, 0, 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.9945493127149, 0, 3.88372109966426, 0, 4.83941449038287, 0, 5.79383105522966, 0, 6.66449205783599, 0, 7.36540280606647, 0, 7.82085387950912, 0, 7.97884560802865, 0, 7.82085387950912, 0, 7.36540280606647, 0, 6.66449205783599, 0, 5.79383105522966, 0, 4.83941449038287, 0, 3.88372109966426, 0, 2.9945493127149, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## paste(vars.fft.just.t[,2], collapse=", ")
## (if using JTransforms) expectedInitialDs = 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0

vars.fft.just.t.tmp <- fftR.propagate.t.debug(vars.fft.just.t, pars.fft.just.t$lo$lambda, pars.fft.just.t$lo$mu, control.fft.32$dt.max, ext.fft.just.t$ndat[2])

## alternating elements with 0's to match SST's JavaFftService requirements
exp.sp1.es.after.prop.t <- rep(0, length(vars.fft.just.t.tmp[,1])*2)
idxs <- seq(1, length(vars.fft.just.t.tmp[,1])*2, by=2)
j <- 1
for (i in vars.fft.just.t.tmp[,1]) {
    exp.sp1.es.after.prop.t[idxs[j]] = i
    j = j+1
}

paste(exp.sp1.es.after.prop.t, collapse=", ")
## (if using SST) expectedSp1EsAfterPropT = 0.000299740440744498, 0, 0.000299739520470888, 0, 0.000299738597335782, 0, 0.00029973767161558, 0, 0.000299736743589792, 0, 0.000299735813540893, 0, 0.000299734881753716, 0, 0.000299733948515344, 0, 0.00029973301411462, 0, 0.00029973207884177, 0, 0.000299731142988294, 0, 0.000299730206846208, 0, 0.00029972927070804, 0, 0.000299728334866242, 0, 0.000299727399612858, 0, 0.000299726465239364, 0, 0.000299725532035927, 0, 0.000299724600291355, 0, 0.000299723670292635, 0, 0.000299722742324704, 0, 0.000299721816669763, 0, 0.000299720893607378, 0, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## paste(vars.fft.just.t.tmp[,1], collapse=", ")
## (if using JTransforms) expectedSp1EsAfterPropT = 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0

## alternating elements with 0's to match SST's JavaFftService requirements
exp.sp1.ds.after.prop.t <- rep(0, length(vars.fft.just.t.tmp[,2])*2)
idxs <- seq(1, length(vars.fft.just.t.tmp[,2])*2, by=2)
j <- 1
for (i in vars.fft.just.t.tmp[,2]) {
    exp.sp1.ds.after.prop.t[idxs[j]] = i
    j = j+1
}

paste(exp.sp1.ds.after.prop.t, collapse=", ")
## (if using SST) expectedSp1DsAfterPropT = 0.708264611249434, 0, 1.07794488915543, 0, 1.57625248872262, 0, 2.21453845459856, 0, 2.98929572353432, 0, 3.87688349797518, 0, 4.83086427268124, 0, 5.78355856417974, 0, 6.65263439389469, 0, 7.35225216820117, 0, 7.80684129110362, 0, 7.96450018578724, 0, 7.80674374052118, 0, 7.35206845754541, 0, 6.65238511637526, 0, 5.78326972079126, 0, 4.83056283595408, 0, 3.87659337350287, 0, 2.98903491521347, 0, 2.21431781312691, 0, 1.57607596745573, 0, 1.07781089197137, 0, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## paste(vars.fft.just.t.tmp[,2], collapse=", ")
## (if using JTransforms) expectedSp1DsAfterPropT = 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052118, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0

pde.fftR.just.t <- with(control.fft.32, make.pde.quasse.fftR.debug(nx, dx, dt.max, 2L))
ans.fftR.just.t <- pde.fftR.just.t(vars.fft.just.t, control.fft.32$dt.max, pars.fft.just.t$lo, 0)

## same as above (just comparing)
paste(ans.fftR.just.t[[2]][,1], collapse=", ")

## (again, if using JTransforms, same expectedSp1EsAfterPropT =) 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0

paste(ans.fftR.just.t[[2]][,2], collapse=", ")
## (again, if using JTransforms, same expectedSp1DsAfterPropT =) 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052118, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0



# (4) testIntegrateOneBranchLoRes32BinsOutsideClassJustX

dt <- 0.01
dx <- 0.01
w <- 10
drift <- 0.0
diffusion <- 0.001
sd <- 0.05
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y
lambda <- sigmoid.x # lambda is now a function
mu <- constant.x # mu is now a function

control.fft.32 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=dt, # dt
                    nx=32, # number of X bins
                    dx=dx, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=w, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

args.fft.just.x <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x <- c(y0, y1, xmid, r, death, drift, diffusion) # specifies parameter values

## args.fft.just.x <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
## pars.just.x <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
## # note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x <- quasse.extent.debug(control.fft.32, drift, diffusion) # prepares X axis stuff
ext.fft.just.x$padding # 16/16, 4/4

pars.fft.just.x <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.just.x, ext.fft.just.x, pars.just.x) # adds lambda and mu vectors to ext.fft.just.x

## initialization
vars.fft.just.x <- matrix(0, control.fft.32$nx, 2) # low resolution

## populating (see previous test for explanation if necessary at this point)
vars.fft.just.x[seq_len(ext.fft.just.x$ndat[2]),2] <- dnorm(ext.fft.just.x$x[[2]], 0.0, sd)

## alternating elements with 0's to match SST's JavaFftService requirements
exp.sp1.initial.ds <- rep(0, length(vars.fft.just.x[,2])*2)
idxs <- seq(1, length(vars.fft.just.x[,2])*2, by=2)
j <- 1
for (i in vars.fft.just.x[,2]) {
    exp.sp1.initial.ds[idxs[j]] = i
    j = j+1
}

paste(exp.sp1.initial.ds, collapse=", ")
## expectedInitialDs = 0.709491856924629, 0, 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.9945493127149, 0, 3.88372109966426, 0, 4.83941449038287, 0, 5.79383105522966, 0, 6.66449205783599, 0, 7.36540280606647, 0, 7.82085387950912, 0, 7.97884560802865, 0, 7.82085387950912, 0, 7.36540280606647, 0, 6.66449205783599, 0, 5.79383105522966, 0, 4.83941449038287, 0, 3.88372109966426, 0, 2.9945493127149, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## paste(vars.fft.just.x[,2], collapse=", ")
## (if using JTransforms) expectedInitialDs = 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0

# checking kernel matches (see PropagatesQuaSSETest -> testMakeNormalKernInPlaceAndFFtAndIfft)
kern.just.x <- fftR.make.kern.debug(-control.fft.32$dt * pars.fft.just.x$hi$drift, sqrt(control.fft.32$dt * pars.fft.just.x$hi$diffusion), control.fft.32$nx, control.fft.32$dx, pars.fft.just.x$hi$padding[1], pars.fft.just.x$hi$padding[2])

paste(kern.just.x, collapse=", ")
## 0.986703287028858, 0.00664835445182386, 2.03374705433156e-09, 2.82445649260927e-20, 1.78085279698565e-35, 5.09772422059472e-55, 6.62490770689586e-79, 3.90875553004076e-107, 1.04701370374391e-139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.04701370374391e-139, 3.90875553004076e-107, 6.62490770689586e-79, 5.09772422059472e-55, 1.78085279698565e-35, 2.82445649260927e-20, 2.03374705433156e-09, 0.00664835445182386

# checking fft-ed fY matches
fy.just.x <- fft(kern.just.x)

## alternating elements with 0's to match SST's JavaFftService requirements
exp.fft.fy <- rep(0, length(fy.just.x)*2)
idxs <- seq(1, length(fy.just.x)*2, by=2)
j <- 1
for (i in fy.just.x) {
    exp.fft.fy[idxs[j]] = as.numeric(i)
    j = j+1
}

paste(exp.fft.fy, collapse=", ")
## 1, 0, 0.999744507157237, 0, 0.998987847110852, 0, 0.997759097982437, 0, 0.996105480062091, 0, 0.994090541136289, 0, 0.991791714355113, 0, 0.989297342492751, 0, 0.986703282961364, 0, 0.984109224049216, 0, 0.981614853950298, 0, 0.979316029808302, 0, 0.977301093995625, 0, 0.975647479188405, 0, 0.97441873269917, 0, 0.973662074416229, 0, 0.973406582192705, 0, 0.973662074416229, 0, 0.97441873269917, 0, 0.975647479188405, 0, 0.977301093995625, 0, 0.979316029808302, 0, 0.981614853950298, 0, 0.984109224049216, 0, 0.986703282961364, 0, 0.989297342492751, 0, 0.991791714355113, 0, 0.994090541136289, 0, 0.996105480062091, 0, 0.997759097982437, 0, 0.998987847110852, 0, 0.999744507157237, 0

## paste(Re(fy.just.x), collapse=", ")
## (if using JTransforms) expectedFFTedfY = 1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216, 0.981614853950298, 0.979316029808302, 0.977301093995625, 0.975647479188405, 0.97441873269917, 0.973662074416229, 0.973406582192705, 0.973662074416229, 0.97441873269917, 0.975647479188405, 0.977301093995625, 0.979316029808302, 0.981614853950298, 0.984109224049216, 0.986703282961364, 0.989297342492751, 0.991791714355113, 0.994090541136289, 0.996105480062091, 0.997759097982437, 0.998987847110852, 0.999744507157237

# checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x, 2, fft)
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x
Re(cstep.2[,2])[1:10] # matches with prints
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # matches with prints

ds.prop.x <- fftR.propagate.x.debug(vars.fft.just.x, control.fft.32$nx, fy.just.x, pars.fft.just.x$lo$padding[1], pars.fft.just.x$lo$padding[2])

## alternating elements with 0's to match SST's JavaFftService requirements
exp.ds.prop.x <- rep(0, length(ds.prop.x[,2])*2)
idxs <- seq(1, length(ds.prop.x[,2])*2, by=2)
j <- 1
for (i in ds.prop.x[,2]) {
    exp.ds.prop.x[idxs[j]] = i
    j = j+1
}

paste(exp.ds.prop.x, collapse=", ")
## expectedSp1DsAfterPropX = 0.709491856924629, 0, 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.99530083804265, 0, 3.88416335936269, 0, 4.83940600155165, 0, 5.79327421787607, 0, 6.66336349661662, 0, 7.36377090119626, 0, 7.81887626199731, 0, 7.97674483551017, 0, 7.81887626199731, 0, 7.36377090119625, 0, 6.66336349661661, 0, 5.79327421787607, 0, 4.83940600155166, 0, 3.88416335936269, 0, 2.99530083804265, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## paste(ds.prop.x[,1], collapse=", ")
## (if using JTransforms) expectedSp1EsAfterPropX = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## paste(ds.prop.x[,2], collapse=", ")
## (if using JTransforms) expectedSp1DsAfterPropX = 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.99530083804265, 3.88416335936269, 4.83940600155165, 5.79327421787607, 6.66336349661662, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119625, 6.66336349661661, 5.79327421787607, 4.83940600155166, 3.88416335936269, 2.99530083804265, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0



## (5) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX (and 10.1 testPropagateChOneCh1024QuaSSETest)

dt <- 0.01
dx <- 0.0005
w <- 10
drift <- 0.0
diffusion <- 0.001
sd <- 0.05
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y
lambda <- sigmoid.x # lambda is now a function
mu <- constant.x # mu is now a function

control.fft.just.x.1024 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=dt, # dt
                    nx=1024, # number of X bins
                    dx=dx, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

args.fft.just.x.1024 <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x.1024 <- c(y0, y1, xmid, r, death, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x.1024 <- quasse.extent.debug(control.fft.just.x.1024, drift, diffusion) # prepares X axis stuff
# ext.fft.just.x.1024$padding

pars.fft.just.x.1024 <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.just.x.1024, ext.fft.just.x.1024, pars.just.x.1024) # adds lambda and mu vectors to ext.fft.just.x

## initialization
vars.fft.just.x.1024 <- matrix(0, control.fft.just.x.1024$nx, 2) # low resolution

## populating (see previous test for explanation if necessary at this point)
vars.fft.just.x.1024[seq_len(ext.fft.just.x.1024$ndat[2]),2] <- dnorm(ext.fft.just.x.1024$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.1024$lo$padding should be 0.0

## alternating elements with 0's to match SST's JavaFftService requirements
exp.sp1.initial.ds <- rep(0, length(vars.fft.just.x.1024[,2])*2)
idxs <- seq(1, length(vars.fft.just.x.1024[,2])*2, by=2)
j <- 1
for (i in vars.fft.just.x.1024[,2]) {
    exp.sp1.initial.ds[idxs[j]] = i
    j = j+1
}

print(paste(exp.sp1.initial.ds[1:96], collapse=","))
## expectedInitialDs = 0.000365714984190948,0,0.000382414193856355,0,0.000399835934138456,0,0.00041800955800901,0,0.000436965523466329,0,0.000456735431302938,0,0.000477352064003592,0,0.000498849425801072,0,0.000521262783917567,0,0.000544628711019852,0,0.000568985128916886,0,0.000594371353528846,0,0.000620828141157005,0,0.000648397736084276,0,0.000677123919536558,0,0.000707052060035464,0,0.000738229165173324,0,0.000770703934841743,0,0.000804526815945299,0,0.000839750058632349,0,0.000876427774075162,0,0.000914615993832026,0,0.000954372730824099,0,0.000995758041960245,0,0.0010388340924432,0,0.0010836652217908,0,0.00113031801160615,0,0.0011788613551308,0,0.00122936652861539,0,0.00128190726454212,0,0.00133655982673381,0,0.00139340308738429,0,0.00145251860604505,0,0.00151399071060322,0,0.00157790658028586,0,0.00164435633072572,0,0.00171343310112364,0,0.00178523314354266,0,0.00185985591436892,0,0.00193740416797439,0,0.00201798405261629,0,0.00210170520860801,0,0.00218868086879601,0,0.00227902796137729,0,0.00237286721509132,0,0.00247032326682047,0,0.00257152477163242,0,0.00267660451529771,0

## print(paste(vars.fft.just.x.1024[1:48,2], collapse=","))
## (if using JTransforms) expectedInitialDs = 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771

kern.just.x.1024 <- fftR.make.kern.debug(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$lo$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$lo$diffusion), control.fft.just.x.1024$nx, control.fft.just.x.1024$dx, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])

fy.just.x.1024 <- fft(kern.just.x.1024)

## alternating elements with 0's to match SST's JavaFftService requirements
exp.fft.fy.1024 <- rep(0, length(fy.just.x.1024)*2)
idxs <- seq(1, length(fy.just.x.1024)*2, by=2)
j <- 1
for (i in Re(fy.just.x.1024)) {
    exp.fft.fy.1024[idxs[j]] = i
    j = j+1
}

paste(Re(exp.fft.fy.1024)[1:96], collapse=", ")
## expectedFFTedfY = 1, 0, 0.999247292368191, 0, 0.996992567179915, 0, 0.993245992007481, 0, 0.988024427909734, 0, 0.981351303026827, 0, 0.973256437474405, 0, 0.963775821368549, 0, 0.952951348293497, 0, 0.940830506973933, 0, 0.92746603432656, 0, 0.912915533436717, 0, 0.897241060330147, 0, 0.880508683684162, 0, 0.862788021843104, 0, 0.844151761668242, 0, 0.824675163860569, 0, 0.804435559446082, 0, 0.783511842107391, 0, 0.761983960984176, 0, 0.739932418450195, 0, 0.717437777209001, 0, 0.694580180837756, 0, 0.671438891652661, 0, 0.64809184947515, 0, 0.624615254550175, 0, 0.601083177512084, 0, 0.577567198915383, 0, 0.554136080452835, 0, 0.530855469577788, 0, 0.507787638837061, 0, 0.484991260810779, 0, 0.462521219151724, 0, 0.440428455824044, 0, 0.418759854264325, 0, 0.39755815783139, 0, 0.376861922578424, 0, 0.356705503075537, 0, 0.3371190697353, 0, 0.318128655850257, 0, 0.299756232341542, 0, 0.282019808042489, 0, 0.264933553200873, 0, 0.248507943778199, 0, 0.232749924053504, 0, 0.217663085001486, 0, 0.203247855908884, 0, 0.189501706717032, 0

## print(paste(Re(fy.just.x.1024)[1:48], collapse=", "))
## (if using JTransforms) expectedFFTedfY = 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368548, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577788, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032

## checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x.1024, 2, fft)
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x.1024
Re(cstep.2[,2])[1:10] # matches with prints
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # matches with prints (depending on CPU architecture, this differs if numbers get too small!)

ds.prop.x.1024 <- fftR.propagate.x.debug(vars.fft.just.x.1024, control.fft.just.x.1024$nx, fy.just.x.1024, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])

## alternating elements with 0's to match SST's JavaFftService requirements
exp.es.x.1024 <- rep(0, length(ds.prop.x.1024[1:48,1])*2)
idxs <- seq(1, length(ds.prop.x.1024[1:48,1])*2, by=2)
j <- 1
for (i in ds.prop.x.1024[1:48,1]) {
    exp.es.x.1024[idxs[j]] = as.numeric(i)
    j = j+1
}

print(paste(exp.es.x.1024, collapse=", "))
## expectedSp1EsAfterPropX = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## alternating elements with 0's to match SST's JavaFftService requirements
exp.ds.x.1024.first48 <- rep(0, length(ds.prop.x.1024[1:48,2])*2)
idxs <- seq(1, length(ds.prop.x.1024[1:48,2])*2, by=2)
j <- 1
for (i in ds.prop.x.1024[1:48,2]) {
    exp.ds.x.1024.first48[idxs[j]] = as.numeric(i)
    j = j+1
}

print(paste(exp.ds.x.1024.first48, collapse=", "))
## expectedSp1DsFirst48AfterPropX = 0.000365714984190948, 0, 0.000382414193856355, 0, 0.000399835934138456, 0, 0.00041800955800901, 0, 0.000436965523466329, 0, 0.000456735431302938, 0, 0.000477352064003592, 0, 0.000498849425801072, 0, 0.000521262783917567, 0, 0.000544628711019852, 0, 0.000568985128916886, 0, 0.000594371353528846, 0, 0.000620828141157005, 0, 0.000648397736084276, 0, 0.000677123919536558, 0, 0.000707052060035464, 0, 0.000738229165173324, 0, 0.000770703934841743, 0, 0.000804526815945299, 0, 0.000839750058632349, 0, 0.000876427774075162, 0, 0.000914615993832026, 0, 0.000954372730824099, 0, 0.000995758041960245, 0, 0.0010388340924432, 0, 0.0010836652217908, 0, 0.00113031801160615, 0, 0.0011788613551308, 0, 0.00122936652861539, 0, 0.00128190726454212, 0, 0.00133655982673381, 0, 0.00139340308738429, 0, 0.00145251860604505, 0, 0.00151399071060322, 0, 0.00157790658028586, 0, 0.00164435633072572, 0, 0.00171343310112364, 0, 0.00178523314354266, 0, 0.00185985591436892, 0, 0.00193740416797439, 0, 0.00201798405261629, 0, 0.00210170520860801, 0, 0.00218868086879601, 0, 0.00227902796137729, 0, 0.00237286721509132, 0, 0.00247032326682047, 0, 0.00257152477163242, 0, 0.00267660451529771, 0

## print(paste(ds.prop.x.1024[1:48,2], collapse=", "))
## (if using JTransforms) expectedSp1DsFirst48AfterPropX = 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771

## alternating elements with 0's to match SST's JavaFftService requirements
exp.ds.x.1024.mid48 <- rep(0, length(ds.prop.x.1024[848:895,2])*2)
idxs <- seq(1, length(ds.prop.x.1024[848:895,2])*2, by=2)
j <- 1
for (i in ds.prop.x.1024[848:895,2]) {
    exp.ds.x.1024.mid48[idxs[j]] = as.numeric(i)
    j = j+1
}

print(paste(exp.ds.x.1024.mid48, collapse=", "))
## expectedSp1DsLater48AfterPropX = 0.00267660451529771, 0, 0.00257152477163242, 0, 0.00247032326682047, 0, 0.00237286721509132, 0, 0.0022790279613773, 0, 0.00218868086879601, 0, 0.00210170520860801, 0, 0.0020179840526163, 0, 0.00193740416797439, 0, 0.00185985591436892, 0, 0.00178523314354266, 0, 0.00171343310112364, 0, 0.00164435633072573, 0, 0.00157790658028586, 0, 0.00151399071060322, 0, 0.00145251860604505, 0, 0.00139340308738429, 0, 0.00133655982673381, 0, 0.00128190726454212, 0, 0.00122936652861539, 0, 0.0011788613551308, 0, 0.00113031801160615, 0, 0.0010836652217908, 0, 0.0010388340924432, 0, 0.000995758041960245, 0, 0.000954372730824099, 0, 0.000914615993832026, 0, 0.000876427774075162, 0, 0.000839750058632349, 0, 0.000804526815945299, 0, 0.000770703934841743, 0, 0.000738229165173324, 0, 0.000707052060035464, 0, 0.000677123919536558, 0, 0.000648397736084276, 0, 0.000620828141157005, 0, 0.000594371353528846, 0, 0.000568985128916886, 0, 0.000544628711019852, 0, 0.000521262783917567, 0, 0.000498849425801072, 0, 0.000477352064003592, 0, 0.000456735431302938, 0, 0.000436965523466329, 0, 0.00041800955800901, 0, 0.000399835934138456, 0, 0.000382414193856355, 0, 0.000365714984190948, 0

## print(paste(ds.prop.x.1024[848:895,2], collapse=", "))
## (if using JTransforms) expectedSp1DsLast48AfterPropX = 0.00267660451529771, 0.00257152477163242, 0.00247032326682047, 0.00237286721509132, 0.0022790279613773, 0.00218868086879601, 0.00210170520860801, 0.0020179840526163, 0.00193740416797439, 0.00185985591436892, 0.00178523314354266, 0.00171343310112364, 0.00164435633072573, 0.00157790658028586, 0.00151399071060322, 0.00145251860604505, 0.00139340308738429, 0.00133655982673381, 0.00128190726454212, 0.00122936652861539, 0.0011788613551308, 0.00113031801160615, 0.0010836652217908, 0.0010388340924432, 0.000995758041960245, 0.000954372730824099, 0.000914615993832026, 0.000876427774075162, 0.000839750058632349, 0.000804526815945299, 0.000770703934841743, 0.000738229165173324, 0.000707052060035464, 0.000677123919536558, 0.000648397736084276, 0.000620828141157005, 0.000594371353528846, 0.000568985128916886, 0.000544628711019852, 0.000521262783917567, 0.000498849425801072, 0.000477352064003592, 0.000456735431302938, 0.000436965523466329, 0.00041800955800901, 0.000399835934138456, 0.000382414193856355, 0.000365714984190948



## (6) testIntegrateOneBranchHiRes4096BinsOutsideClassJustX (note that I'm using 4096 as low res here, but I'll co-opt 4096 high res in Java)  (and 11.1 testPropagateChOneCh1024QuaSSETest)
##
## (builds on [4] as it depends on control.fft.just.x.1024 and pars.fft.just.x.1024)

## initialization
vars.fft.just.x.4096 <- matrix(0, control.fft.just.x.1024$nx * control.fft.just.x.1024$r, 2) # low resolution

## populating (see previous test for explanation if necessary at this point)
vars.fft.just.x.4096[seq_len(ext.fft.just.x.1024$ndat[1]),2] <- dnorm(ext.fft.just.x.1024$x[[1]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.4096$lo$padding should be 0.0

exp.sp1.initial.ds <- rep(0, length(vars.fft.just.x.4096[1:48,2])*2)
idxs <- seq(1, length(vars.fft.just.x.4096[1:48,2])*2, by=2)
j <- 1
for (i in vars.fft.just.x.4096[1:48,2]) {
    exp.sp1.initial.ds[idxs[j]] = i
    j = j+1
}

print(paste(exp.sp1.initial.ds[1:96], collapse=","))
## expectedInitialDsFirst48 = 0.000353647683540531,0,0.000357627448646487,0,0.000361649739618714,0,0.000365714984190948,0,0.000369823614096463,0,0.000373976065102163,0,0.000378172777042955,0,0.000382414193856355,0,0.000386700763617376,0,0.000391032938573662,0,0.000395411175180895,0,0.000399835934138456,0,0.000404307680425366,0,0.000408826883336479,0,0.000413394016518956,0,0.00041800955800901,0,0.000422673990268911,0,0.00042738780022428,0,0.000432151479301657,0,0.000436965523466329,0,0.000441830433260469,0,0.000446746713841521,0,0.000451714875020898,0,0.000456735431302938,0,0.000461808901924177,0,0.00046693581089287,0,0.000472116687028841,0,0.000477352064003592,0,0.000482642480380733,0,0.000487988479656683,0,0.000493390610301674,0,0.000498849425801072,0,0.000504365484696964,0,0.000509939350630071,0,0.000515571592381964,0,0.000521262783917567,0,0.00052701350442799,0,0.000532824338373647,0,0.000538695875527712,0,0.000544628711019852,0,0.000550623445380317,0,0.000556680684584297,0,0.000562801040096648,0,0.000568985128916886,0,0.000575233573624552,0,0.000581547002424852,0,0.000587926049194663,0,0.000594371353528846,0

## paste(vars.fft.just.x.4096[1:48,2], collapse=", ")
## (if using JTransforms) expectedInitialDsFirst48 = 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846

kern.just.x.4096 <- fftR.make.kern.debug(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$diffusion), control.fft.just.x.1024$nx * control.fft.just.x.1024$r, control.fft.just.x.1024$dx / control.fft.just.x.1024$r, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2]) # note the division of dx / r to do high resolution!

fy.just.x.4096 <- fft(kern.just.x.4096)

## alternating elements with 0's to match SST's JavaFftService requirements
exp.fft.fy.4096 <- rep(0, length(fy.just.x.4096[1:96])*2)
idxs <- seq(1, length(fy.just.x.4096[1:96])*2, by=2)
j <- 1
for (i in Re(fy.just.x.4096[1:96])) {
    exp.fft.fy.4096[idxs[j]] = i
    j = j+1
}

paste(exp.fft.fy.4096[1:192], collapse=", ")
## expectedFFTedfY = 1, 0, 0.999247292368191, 0, 0.996992567179915, 0, 0.993245992007481, 0, 0.988024427909734, 0, 0.981351303026827, 0, 0.973256437474405, 0, 0.963775821368549, 0, 0.952951348293497, 0, 0.940830506973933, 0, 0.92746603432656, 0, 0.912915533436717, 0, 0.897241060330147, 0, 0.880508683684162, 0, 0.862788021843104, 0, 0.844151761668242, 0, 0.824675163860569, 0, 0.804435559446082, 0, 0.783511842107391, 0, 0.761983960984176, 0, 0.739932418450195, 0, 0.717437777209001, 0, 0.694580180837756, 0, 0.671438891652661, 0, 0.64809184947515, 0, 0.624615254550175, 0, 0.601083177512084, 0, 0.577567198915383, 0, 0.554136080452835, 0, 0.530855469577789, 0, 0.507787638837061, 0, 0.484991260810779, 0, 0.462521219151724, 0, 0.440428455824044, 0, 0.418759854264325, 0, 0.39755815783139, 0, 0.376861922578424, 0, 0.356705503075537, 0, 0.3371190697353, 0, 0.318128655850257, 0, 0.299756232341542, 0, 0.282019808042489, 0, 0.264933553200873, 0, 0.248507943778199, 0, 0.232749924053504, 0, 0.217663085001486, 0, 0.203247855908884, 0, 0.189501706717032, 0, 0.17641935863023, 0, 0.163993000606569, 0, 0.152212509447405, 0, 0.14106567132108, 0, 0.130538402692698, 0, 0.120614968781996, 0, 0.111278197832449, 0, 0.102509689644048, 0, 0.0942900169966733, 0, 0.0865989187681513, 0, 0.0794154837283513, 0, 0.0727183241657203, 0, 0.0664857386734408, 0, 0.0606958635869983, 0, 0.0553268127217662, 0, 0.0503568052068668, 0, 0.045764281348922, 0, 0.0415280065854677, 0, 0.0376271637021404, 0, 0.0340414335898102, 0, 0.0307510649074428, 0, 0.0277369330935997, 0, 0.0249805892342995, 0, 0.02246429934781, 0, 0.0201710746882624, 0, 0.0180846937003869, 0, 0.0161897162778361, 0, 0.0144714909882638, 0, 0.0129161559303678, 0, 0.0115106338823706, 0, 0.0102426223887588, 0, 0.00910057941344994, 0, 0.00807370516375022, 0, 0.0071519206614005, 0, 0.0063258436054819, 0, 0.00558676203776515, 0, 0.00492660628496708, 0, 0.00433791961500649, 0, 0.00381382800634398, 0, 0.00334800939141312, 0, 0.00293466269748287, 0, 0.0025684769714751, 0, 0.00224460083965355, 0, 0.00195861251901208, 0, 0.00170649056486081, 0, 0.001484585508736, 0, 0.00128959251247564, 0, 0.00111852513821068, 0

## print(paste(Re(fy.just.x.4096)[1:48], collapse=", "))
## (if using JTransforms) expectedFFTedfY = 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368549, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577789, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032

ds.prop.x.4096 <- fftR.propagate.x.debug(vars.fft.just.x.4096, control.fft.just.x.1024$nx * control.fft.just.x.1024$r, fy.just.x.4096, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2])

## alternating elements with 0's to match SST's JavaFftService requirements
exp.es.x.4096 <- rep(0, length(ds.prop.x.4096[1:48,1])*2)
idxs <- seq(1, length(ds.prop.x.4096[1:48,1])*2, by=2)
j <- 1
for (i in ds.prop.x.4096[1:48,1]) {
    exp.es.x.4096[idxs[j]] = as.numeric(i)
    j = j+1
}

print(paste(exp.es.x.4096[1:96], collapse=", "))
## expectedSp1EsAfterPropX = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## alternating elements with 0's to match SST's JavaFftService requirements
exp.ds.x.4096 <- rep(0, length(ds.prop.x.4096[1:48,2])*2)
idxs <- seq(1, length(ds.prop.x.4096[1:48,2])*2, by=2)
j <- 1
for (i in ds.prop.x.4096[1:48,2]) {
    exp.ds.x.4096[idxs[j]] = as.numeric(i)
    j = j+1
}

print(paste(exp.ds.x.4096[1:96], collapse=", "))
## expectedSp1DsFirst48AfterPropX = 0.000353647683540531, 0, 0.000357627448646487, 0, 0.000361649739618714, 0, 0.000365714984190948, 0, 0.000369823614096463, 0, 0.000373976065102163, 0, 0.000378172777042955, 0, 0.000382414193856355, 0, 0.000386700763617376, 0, 0.000391032938573662, 0, 0.000395411175180895, 0, 0.000399835934138456, 0, 0.000404307680425366, 0, 0.000408826883336479, 0, 0.000413394016518956, 0, 0.00041800955800901, 0, 0.000422673990268911, 0, 0.00042738780022428, 0, 0.000432151479301657, 0, 0.000436965523466329, 0, 0.000441830433260469, 0, 0.000446746713841521, 0, 0.000451714875020898, 0, 0.000456735431302938, 0, 0.000461808901924177, 0, 0.00046693581089287, 0, 0.000472116687028841, 0, 0.000477352064003592, 0, 0.000482642480380733, 0, 0.000487988479656683, 0, 0.000493390610301674, 0, 0.000498849425801072, 0, 0.000504365484696964, 0, 0.000509939350630071, 0, 0.000515571592381964, 0, 0.000521262783917567, 0, 0.00052701350442799, 0, 0.000532824338373647, 0, 0.000538695875527712, 0, 0.000544628711019852, 0, 0.000550623445380317, 0, 0.000556680684584297, 0, 0.000562801040096648, 0, 0.000568985128916886, 0, 0.000575233573624552, 0, 0.000581547002424852, 0, 0.000587926049194663, 0, 0.000594371353528846, 0

## print(paste(ds.prop.x.4096[1:48,2], collapse=", "))
## (if using JTransforms) expectedSp1DsFirst48AfterPropX = 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846

## alternating elements with 0's to match SST's JavaFftService requirements
exp.ds.x.4096.mid48 <- rep(0, length(ds.prop.x.4096[1001:1048,2])*2)
idxs <- seq(1, length(ds.prop.x.4096[1001:1048,2])*2, by=2)
j <- 1
for (i in ds.prop.x.4096[1001:1048,2]) {
    exp.ds.x.4096.mid48[idxs[j]] = as.numeric(i)
    j = j+1
}

print(paste(exp.ds.x.4096.mid48, collapse=", "))
## expectedSp1DsLater48AfterPropX = 1.13579683037393, 0, 1.14139979445304, 0, 1.14702325798403, 0, 1.15266725177838, 0, 1.15833180641492, 0, 1.16401695223792, 0, 1.16972271935522, 0, 1.17544913763623, 0, 1.18119623671012, 0, 1.18696404596383, 0, 1.19275259454019, 0, 1.19856191133597, 0, 1.20439202499999, 0, 1.21024296393119, 0, 1.21611475627669, 0, 1.22200742992987, 0, 1.22792101252845, 0, 1.23385553145257, 0, 1.23981101382282, 0, 1.24578748649837, 0, 1.25178497607497, 0, 1.25780350888308, 0, 1.26384311098589, 0, 1.26990380817741, 0, 1.27598562598051, 0, 1.28208858964501, 0, 1.28821272414573, 0, 1.29435805418053, 0, 1.30052460416843, 0, 1.30671239824758, 0, 1.31292146027342, 0, 1.31915181381665, 0, 1.32540348216137, 0, 1.33167648830306, 0, 1.33797085494669, 0, 1.34428660450477, 0, 1.35062375909539, 0, 1.35698234054031, 0, 1.36336237036299, 0, 1.36976386978664, 0, 1.37618685973233, 0, 1.38263136081698, 0, 1.38909739335148, 0, 1.39558497733871, 0, 1.40209413247163, 0, 1.40862487813131, 0, 1.41517723338502, 0, 1.42175121698428, 0

## print(paste(ds.prop.x.4096[1001:1048,2], collapse=", "))
## (if using JTransforms) expectedSp1DsLater48AfterPropX = 1.12964601442883, 1.13523957985681, 1.14085371385815, 1.14648844781881, 1.1521438128928, 1.15781984000018, 1.16351655982519, 1.16923400281429, 1.17497219917421, 1.18073117887008, 1.18651097162341, 1.19231160691021, 1.19813311395904, 1.20397552174904, 1.20983885900802, 1.21572315421049, 1.22162843557575, 1.22755473106588, 1.23350206838385, 1.23947047497152, 1.24545997800775, 1.25147060440637, 1.25750238081428, 1.26355533360948, 1.26962948889911, 1.27572487251749, 1.28184151002416, 1.28797942670193, 1.29413864755491, 1.30031919730656, 1.30652110039772, 1.31274438098464, 1.31898906293703, 1.32525516983611, 1.33154272497262, 1.33785175134486, 1.34418227165675, 1.35053430831583, 1.35690788343134, 1.36330301881221, 1.36971973596514, 1.37615805609261, 1.38261800009092, 1.38909958854822, 1.39560284174258, 1.40212777963999, 1.40867442189243, 1.41524278783588



## (7) testBothXandTPropagateMethodsInsideClassOneBranchLoRes32Bins

dt <- 0.01
dx <- 0.01
w <- 10
drift <- 0.0
diffusion <- 0.001
sd <- 0.05
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y
lambda <- sigmoid.x # lambda is now a function
mu <- constant.x # mu is now a function

control.fft.32 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=dt, # dt
                    nx=32, # number of X bins
                    dx=dx, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=w, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

args.fft.32.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.32.both <- c(y0, y1, xmid, r, death, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.32.both <- quasse.extent.debug(control.fft.32, drift, diffusion) # prepares X axis stuff

## initialization
vars.fft.32.both <- matrix(0, control.fft.32$nx, 2) # low resolution

## populating (see previous test for explanation if necessary at this point)
vars.fft.32.both[seq_len(ext.fft.32.both$ndat[2]),2] <- dnorm(ext.fft.32.both$x[[2]], 0.0, sd)

exp.sp1.initial.ds <- rep(0, length(vars.fft.32.both[,2])*2)
idxs <- seq(1, length(vars.fft.32.both[,2])*2, by=2)
j <- 1
for (i in vars.fft.32.both[,2]) {
    exp.sp1.initial.ds[idxs[j]] = i
    j = j+1
}

print(paste(exp.sp1.initial.ds, collapse=", "))
## expectedInitialDs = 0.709491856924629, 0, 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.9945493127149, 0, 3.88372109966426, 0, 4.83941449038287, 0, 5.79383105522966, 0, 6.66449205783599, 0, 7.36540280606647, 0, 7.82085387950912, 0, 7.97884560802865, 0, 7.82085387950912, 0, 7.36540280606647, 0, 6.66449205783599, 0, 5.79383105522966, 0, 4.83941449038287, 0, 3.88372109966426, 0, 2.9945493127149, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## (if using JTransforms) paste(c(vars.fft.32.both[seq_len(ext.fft.32.both$ndat[2]),2], rep(0, (32-length(ext.fft.32.both$x[[2]])))), collapse=", ")
## expectedInitialDs = 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0

pars.fft.32.both <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.32.both, ext.fft.32.both, pars.32.both)

pde.fftR.32.both <- with(control.fft.32, make.pde.quasse.fftR.debug2(nx, dx, dt.max, 2L)) # creates function

ans.fftR.32.both <- pde.fftR.32.both(vars.fft.32.both, control.fft.32$dt.max, pars.fft.32.both$lo, 0)

## alternating elements with 0's to match SST's JavaFftService requirements
fft.ft.32 <- c(1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216, 0.981614853950298, 0.979316029808302, 0.977301093995625, 0.975647479188405, 0.97441873269917, 0.973662074416229, 0.973406582192705, 0.973662074416229, 0.97441873269917, 0.975647479188405, 0.977301093995625, 0.979316029808302, 0.981614853950298, 0.984109224049216, 0.986703282961364, 0.989297342492751, 0.991791714355113, 0.994090541136289, 0.996105480062091, 0.997759097982437, 0.998987847110852, 0.999744507157237)

exp.fft.fy.32 <- rep(0, length(fft.ft.32)*2)
idxs <- seq(1, length(fft.ft.32)*2, by=2)
j <- 1
for (i in fft.ft.32) {
    exp.fft.fy.32[idxs[j]] = i
    j = j+1
}

paste(exp.fft.fy.32, collapse=", ")
## expectedFFTedfY = 1, 0, 0.999744507157237, 0, 0.998987847110852, 0, 0.997759097982437, 0, 0.996105480062091, 0, 0.994090541136289, 0, 0.991791714355113, 0, 0.989297342492751, 0, 0.986703282961364, 0, 0.984109224049216, 0, 0.981614853950298, 0, 0.979316029808302, 0, 0.977301093995625, 0, 0.975647479188405, 0, 0.97441873269917, 0, 0.973662074416229, 0, 0.973406582192705, 0, 0.973662074416229, 0, 0.97441873269917, 0, 0.975647479188405, 0, 0.977301093995625, 0, 0.979316029808302, 0, 0.981614853950298, 0, 0.984109224049216, 0, 0.986703282961364, 0, 0.989297342492751, 0, 0.991791714355113, 0, 0.994090541136289, 0, 0.996105480062091, 0, 0.997759097982437, 0, 0.998987847110852, 0, 0.999744507157237, 0

## (if using JTransforms)
## expectedFFTedfY = fy = 1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216, 0.981614853950298, 0.979316029808302, 0.977301093995625, 0.975647479188405, 0.97441873269917, 0.973662074416229, 0.973406582192705, 0.973662074416229, 0.97441873269917, 0.975647479188405, 0.977301093995625, 0.979316029808302, 0.981614853950298, 0.984109224049216, 0.986703282961364, 0.989297342492751, 0.991791714355113, 0.994090541136289, 0.996105480062091, 0.997759097982437, 0.998987847110852, 0.999744507157237

## E's (note how first 4 and last 4 do not change after propagate in t because they are flanking bins on the left and right)

es.32.t <- c(0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0)

exp.es.32.t <- rep(0, length(es.32.t)*2)
idxs <- seq(1, length(es.32.t)*2, by=2)
j <- 1
for (i in es.32.t) {
    exp.es.32.t[idxs[j]] = i
    j = j+1
}

paste(exp.es.32.t, collapse=", ")
## expectedSp1EsAfterPropT = 0.000299740440744498, 0, 0.000299739520470888, 0, 0.000299738597335782, 0, 0.00029973767161558, 0, 0.000299736743576341, 0, 0.000299735813529336, 0, 0.000299734881744068, 0, 0.000299733948507616, 0, 0.000299733014108822, 0, 0.000299732078837909, 0, 0.000299731142986375, 0, 0.000299730206846234, 0, 0.000299729270710011, 0, 0.000299728334870154, 0, 0.000299727399618708, 0, 0.000299726465247143, 0, 0.000299725532045626, 0, 0.000299724600302962, 0, 0.000299723670306136, 0, 0.000299722742324704, 0, 0.000299721816669763, 0, 0.000299720893607378, 0, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

exp.es.32.t.and.x <- rep(0, length(ans.fftR.32.both[[2]][,1])*2)
idxs <- seq(1, length(ans.fftR.32.both[[2]][,1])*2, by=2)
j <- 1
for (i in ans.fftR.32.both[[2]][,1]) {
    exp.es.32.t.and.x[idxs[j]] = i
    j = j+1
}

paste(exp.es.32.t.and.x, collapse=", ")
## expectedSp1EsAfterPropTandX = 0.000299740440744498, 0, 0.000299739520470888, 0, 0.000299738597335782, 0, 0.00029973767161558, 0, 0.000299736743576341, 0, 0.000299735813529336, 0, 0.000299734881744068, 0, 0.000299733948507616, 0, 0.000299733014108822, 0, 0.000299732078837909, 0, 0.000299731142986375, 0, 0.000299730206846234, 0, 0.000299729270710011, 0, 0.000299728334870154, 0, 0.000299727399618708, 0, 0.000299726465247143, 0, 0.000299725532045626, 0, 0.000299724600302962, 0, 0.000299723670306136, 0, 0.000299722742324704, 0, 0.000299721816669763, 0, 0.000299720893607378, 0, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## (if using JTransforms) paste(ans.fftR.32.both[[2]][,1], collapse=", ")
## expectedSp1EsAfterPropTandX = 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743576341, 0.000299735813529336, 0.000299734881744068, 0.000299733948507616, 0.000299733014108822, 0.000299732078837909, 0.000299731142986375, 0.000299730206846234, 0.000299729270710011, 0.000299728334870154, 0.000299727399618708, 0.000299726465247143, 0.000299725532045626, 0.000299724600302962, 0.000299723670306136, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0

## D's (note how first 4 and last 4 do not change after propagate in t because they are flanking bins on the left and right)

ds.32.t <- c(0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052118, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0)

exp.ds.32.t <- rep(0, length(ds.32.t)*2)
idxs <- seq(1, length(ds.32.t)*2, by=2)
j <- 1
for (i in ds.32.t) {
    exp.ds.32.t[idxs[j]] = i
    j = j+1
}

paste(exp.ds.32.t, collapse=", ")
## expectedSp1DsAfterPropT = 0.708264611249434, 0, 1.07794488915543, 0, 1.57625248872262, 0, 2.21453845459856, 0, 2.98929572353432, 0, 3.87688349797518, 0, 4.83086427268124, 0, 5.78355856417974, 0, 6.65263439389469, 0, 7.35225216820117, 0, 7.80684129110362, 0, 7.96450018578724, 0, 7.80674374052118, 0, 7.35206845754541, 0, 6.65238511637526, 0, 5.78326972079126, 0, 4.83056283595408, 0, 3.87659337350287, 0, 2.98903491521347, 0, 2.21431781312691, 0, 1.57607596745573, 0, 1.07781089197137, 0, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

exp.ds.32.t.and.x <- rep(0, length(ans.fftR.32.both[[2]][,2])*2)
idxs <- seq(1, length(ans.fftR.32.both[[2]][,2])*2, by=2)
j <- 1
for (i in (exp(ans.fftR.32.both[[1]]) * ans.fftR.32.both[[2]][,2])) {
    exp.ds.32.t.and.x[idxs[j]] = i
    j = j+1
}

paste(exp.ds.32.t.and.x, collapse=", ") ## note the normalization factor being added back in through multiplication (after un-log-ing)
## expectedSp1DsAfterPropTandX = 0.708264611249434, 0, 1.07794488915543, 0, 1.57625248872262, 0, 2.21453845459856, 0, 2.99004586159941, 0, 3.87732490267096, 0, 4.83085571964462, 0, 5.78300263831972, 0, 6.65150777531997, 0, 7.35062312893061, 0, 7.80486719135139, 0, 7.96240319031702, 0, 7.80476971649718, 0, 7.35043955518597, 0, 6.65125867066456, 0, 5.78271397425439, 0, 4.83055444185051, 0, 3.87703489789509, 0, 2.98978512541793, 0, 2.21431781312691, 0, 1.57607596745573, 0, 1.07781089197137, 0, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## (if using JTransforms) print(paste(exp(ans.fftR.32.both[[1]]) * ans.fftR.32.both[[2]][,2], collapse=", "))
## expectedSp1DsAfterPropTandX = 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.99004586159941, 3.87732490267096, 4.83085571964462, 5.78300263831972, 6.65150777531997, 7.35062312893061, 7.80486719135139, 7.96240319031702, 7.80476971649718, 7.35043955518597, 6.65125867066456, 5.78271397425439, 4.83055444185051, 3.87703489789509, 2.98978512541793, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0



## (8) testBothXandTPropagateMethodsInsideClassOneBranchLoRes32BinsDt002

dt <- 0.02
dx <- 0.01
w <- 10
drift <- 0.0
diffusion <- 0.001
sd <- 0.05
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y
lambda <- sigmoid.x # lambda is now a function
mu <- constant.x # mu is now a function

control.fft.32.dt002 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=dt, # dt
                    nx=32, # number of X bins
                    dx=dx, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=w, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

args.fft.32.dt002.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.32.dt002.both <- c(y0, y1, xmid, r, death, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.32.dt002.both <- quasse.extent.debug(control.fft.32.dt002, drift, diffusion) # prepares X axis stuff

## initialization
vars.fft.32.dt002.both <- matrix(0, control.fft.32.dt002$nx, 2) # low resolution

## populating (see previous tests for explanation if necessary at this point)
vars.fft.32.dt002.both[seq_len(ext.fft.32.dt002.both$ndat[2]),2] <- dnorm(ext.fft.32.dt002.both$x[[2]], 0.0, sd)

exp.sp1.initial.ds.dt002 <- rep(0, length(vars.fft.32.dt002.both[,2])*2)
idxs <- seq(1, length(vars.fft.32.dt002.both[,2])*2, by=2)
j <- 1
for (i in vars.fft.32.dt002.both[,2]) {
    exp.sp1.initial.ds.dt002[idxs[j]] = i
    j = j+1
}

print(paste(exp.sp1.initial.ds.dt002, collapse=", " ))
## expectedInitialDs = 1.07981933026376, 0, 1.57900316601788, 0, 2.21841669358911, 0, 2.9945493127149, 0, 3.88372109966426, 0, 4.83941449038287, 0, 5.79383105522965, 0, 6.66449205783599, 0, 7.36540280606647, 0, 7.82085387950912, 0, 7.97884560802865, 0, 7.82085387950912, 0, 7.36540280606647, 0, 6.66449205783599, 0, 5.79383105522965, 0, 4.83941449038287, 0, 3.88372109966426, 0, 2.9945493127149, 0, 2.21841669358911, 0, 1.57900316601788, 0, 1.07981933026376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## (if using JTransforms) print(paste(c(vars.fft.32.dt002.both[seq_len(ext.fft.32.dt002.both$ndat[2]),2], rep(0, (32-length(ext.fft.32.dt002.both$x[[2]])))), collapse=", "))
## expectedInitialDs = 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

pars.fft.32.dt002.both <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.32.dt002.both, ext.fft.32.dt002.both, pars.32.dt002.both)

pde.fftR.32.dt002.both <- with(control.fft.32.dt002, make.pde.quasse.fftR.debug2(nx, dx, dt.max, 2L)) # creates function

ans.fftR.32.dt002.both <- pde.fftR.32.dt002.both(vars.fft.32.dt002.both, control.fft.32.dt002$dt.max, pars.fft.32.dt002.both$lo, 0)

## alternating elements with 0's to match SST's JavaFftService requirements
fft.ft.32.dt002 <- c(1, 0.997284635663215, 0.989243568247244, 0.976187735593861, 0.958621745686731, 0.937224046438699, 0.91282033618305, 0.886351315420652, 0.858836097834951, 0.831332753409919, 0.80489754454845, 0.780544437378578, 0.759206428540079, 0.741700129042481, 0.728694899474875, 0.720687643959832, 0.717984152783717, 0.720687643959832, 0.728694899474875, 0.741700129042481, 0.759206428540079, 0.780544437378578, 0.80489754454845, 0.831332753409919, 0.858836097834951, 0.886351315420652, 0.91282033618305, 0.937224046438699, 0.958621745686731, 0.976187735593861, 0.989243568247244, 0.997284635663215)

exp.fft.fy.32.dt002 <- rep(0, length(fft.ft.32.dt002)*2)
idxs <- seq(1, length(fft.ft.32.dt002)*2, by=2)
j <- 1
for (i in fft.ft.32.dt002) {
    exp.fft.fy.32.dt002[idxs[j]] = i
    j = j+1
}

paste(exp.fft.fy.32.dt002, collapse=", ")
## expectedFFTedfY = 1, 0, 0.997284635663215, 0, 0.989243568247244, 0, 0.976187735593861, 0, 0.958621745686731, 0, 0.937224046438699, 0, 0.91282033618305, 0, 0.886351315420652, 0, 0.858836097834951, 0, 0.831332753409919, 0, 0.80489754454845, 0, 0.780544437378578, 0, 0.759206428540079, 0, 0.741700129042481, 0, 0.728694899474875, 0, 0.720687643959832, 0, 0.717984152783717, 0, 0.720687643959832, 0, 0.728694899474875, 0, 0.741700129042481, 0, 0.759206428540079, 0, 0.780544437378578, 0, 0.80489754454845, 0, 0.831332753409919, 0, 0.858836097834951, 0, 0.886351315420652, 0, 0.91282033618305, 0, 0.937224046438699, 0, 0.958621745686731, 0, 0.976187735593861, 0, 0.989243568247244, 0, 0.997284635663215, 0

## (if using JTransforms) expectedFFTedfY = fy = 1, 0.997284635663215, 0.989243568247244, 0.976187735593861, 0.958621745686731, 0.937224046438699, 0.91282033618305, 0.886351315420652, 0.858836097834951, 0.831332753409919, 0.80489754454845, 0.780544437378578, 0.759206428540079, 0.741700129042481, 0.728694899474875, 0.720687643959832, 0.717984152783717, 0.720687643959832, 0.728694899474875, 0.741700129042481, 0.759206428540079, 0.780544437378578, 0.80489754454845, 0.831332753409919, 0.858836097834951, 0.886351315420652, 0.91282033618305, 0.937224046438699, 0.958621745686731, 0.976187735593861, 0.989243568247244, 0.997284635663215

## E's (note how first 5 and last 5 do not change after propagate in t because they are flanking bins on the left and right)

es.32.t.dt002 <- c(0.000598958856744621, 0.000598955169219663, 0.000598951471383546, 0.000598947764353066, 0.000598944049256322, 0.000598940327231526, 0.000598936599425321, 0.000598932866991505, 0.000598929131089762, 0.000598925392884107, 0.000598921653541243, 0.000598917914229581, 0.000598914176117241, 0.000598910440370952, 0.000598906708154486, 0.000598902980627182, 0.000598899258942585, 0.000598895544247008, 0.000598891837677968, 0.000598888140363113, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

exp.es.32.t.dt002 <- rep(0, length(es.32.t.dt002)*2)
idxs <- seq(1, length(es.32.t.dt002)*2, by=2)
j <- 1
for (i in es.32.t.dt002) {
    exp.es.32.t.dt002[idxs[j]] = i
    j = j+1
}

paste(exp.es.32.t.dt002, collapse=", ")
## expectedSp1EsAfterPropT = 0.000598958856744621, 0, 0.000598955169219663, 0, 0.000598951471383546, 0, 0.000598947764353066, 0, 0.000598944049256322, 0, 0.000598940327231526, 0, 0.000598936599425321, 0, 0.000598932866991505, 0, 0.000598929131089762, 0, 0.000598925392884107, 0, 0.000598921653541243, 0, 0.000598917914229581, 0, 0.000598914176117241, 0, 0.000598910440370952, 0, 0.000598906708154486, 0, 0.000598902980627182, 0, 0.000598899258942585, 0, 0.000598895544247008, 0, 0.000598891837677968, 0, 0.000598888140363113, 0, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## (if using JTransforms) expectedSp1EsAfterPropT = E's after propagate in t = 0.000598958856744621, 0.000598955169219663, 0.000598951471383546, 0.000598947764353066, 0.000598944049256322, 0.000598940327231526, 0.000598936599425321, 0.000598932866991505, 0.000598929131089762, 0.000598925392884107, 0.000598921653541243, 0.000598917914229581, 0.000598914176117241, 0.000598910440370952, 0.000598906708154486, 0.000598902980627182, 0.000598899258942585, 0.000598895544247008, 0.000598891837677968, 0.000598888140363113, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

exp.es.32.t.and.x.dt002 <- rep(0, length(ans.fftR.32.dt002.both[[2]][,1])*2)
idxs <- seq(1, length(ans.fftR.32.dt002.both[[2]][,1])*2, by=2)
j <- 1
for (i in ans.fftR.32.dt002.both[[2]][,1]) {
    exp.es.32.t.and.x.dt002[idxs[j]] = i
    j = j+1
}

paste(exp.es.32.t.and.x.dt002, collapse=", ")
## expectedSp1EsAfterPropTandX = 0.000598958856744621, 0, 0.000598955169219663, 0, 0.000598951471383546, 0, 0.000598947764353066, 0, 0.000598944049256322, 0, 0.000598940326823012, 0, 0.000598936599098334, 0, 0.000598932866746461, 0, 0.000598929130926967, 0, 0.000598925392803752, 0, 0.000598921653543447, 0, 0.000598917914314325, 0, 0.000598914176284426, 0, 0.000598910440620369, 0, 0.000598906708485821, 0, 0.000598902981040026, 0, 0.000598899258942585, 0, 0.000598895544247008, 0, 0.000598891837677968, 0, 0.000598888140363113, 0, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## (if using JTransforms)  paste(ans.fftR.32.dt002.both[[2]][,1], collapse=", ")
## expectedSp1EsAfterPropTandX = E's after propagate in t and x = 0.000598958856744621, 0.000598955169219663, 0.000598951471383546, 0.000598947764353066, 0.000598944049256322, 0.000598940326823012, 0.000598936599098334, 0.000598932866746461, 0.000598929130926967, 0.000598925392803752, 0.000598921653543447, 0.000598917914314325, 0.000598914176284426, 0.000598910440620369, 0.000598906708485821, 0.000598902981040026, 0.000598899258942585, 0.000598895544247008, 0.000598891837677968, 0.000598888140363113, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## D's (note how first 5 and last 5 do not change after propagate in t because they are flanking bins on the left and right)

ds.32.t.dt002 <- c(1.07607462857111, 1.57350796410007, 2.2106689156898, 2.98405395410774, 3.87006132459304, 4.82233340368848, 5.77330938669843, 6.64080371901192, 7.33913154877464, 7.79286078082761, 7.95018769794391, 7.79266608851656, 7.33876489743091, 6.64030620874454, 5.77273291052126, 4.82173179361503, 3.86948229161728, 2.98353343056498, 2.21022855749815, 1.5731556614014, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

exp.ds.32.t.dt002 <- rep(0, length(ds.32.t.dt002)*2)
idxs <- seq(1, length(ds.32.t.dt002)*2, by=2)
j <- 1
for (i in ds.32.t.dt002) {
    exp.ds.32.t.dt002[idxs[j]] = i
    j = j+1
}

paste(exp.ds.32.t.dt002, collapse=", ")
## expectedSp1DsAfterPropT = 1.07607462857111, 0, 1.57350796410007, 0, 2.2106689156898, 0, 2.98405395410774, 0, 3.87006132459304, 0, 4.82233340368848, 0, 5.77330938669843, 0, 6.64080371901192, 0, 7.33913154877464, 0, 7.79286078082761, 0, 7.95018769794391, 0, 7.79266608851656, 0, 7.33876489743091, 0, 6.64030620874454, 0, 5.77273291052126, 0, 4.82173179361503, 0, 3.86948229161728, 0, 2.98353343056498, 0, 2.21022855749815, 0, 1.5731556614014, 0, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## (if using JTransforms) expectedSp1DsAfterPropT = D's after propagate in t = 1.07607462857111, 1.57350796410007, 2.2106689156898, 2.98405395410774, 3.87006132459304, 4.82233340368848, 5.77330938669843, 6.64080371901192, 7.33913154877464, 7.79286078082761, 7.95018769794391, 7.79266608851656, 7.33876489743091, 6.64030620874454, 5.77273291052126, 4.82173179361503, 3.86948229161728, 2.98353343056498, 2.21022855749815, 1.5731556614014, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

exp.ds.32.t.and.x.dt002 <- rep(0, length(ans.fftR.32.dt002.both[[2]][,2])*2)
idxs <- seq(1, length(ans.fftR.32.dt002.both[[2]][,2])*2, by=2)
j <- 1
for (i in (exp(ans.fftR.32.dt002.both[[1]]) * ans.fftR.32.dt002.both[[2]][,2])) {
    exp.ds.32.t.and.x.dt002[idxs[j]] = i
    j = j+1
}

paste(exp.ds.32.t.and.x.dt002, collapse=", ") ## note the normalization factor being added back in through multiplication (after un-log-ing)
## expectedSp1DsAfterPropTandX = 1.07607462857111, 0, 1.57350796410007, 0, 2.2106689156898, 0, 2.98405395410774, 0, 3.87006132459304, 0, 4.82224125131593, 0, 5.76741044168543, 0, 6.62885082364537, 0, 7.32184914756387, 0, 7.77191831137883, 0, 7.92794195899393, 0, 7.77172522522963, 0, 7.3214854000561, 0, 6.62835697979005, 0, 5.76683776896087, 0, 4.8216430122779, 0, 3.86948229161728, 0, 2.98353343056498, 0, 2.21022855749815, 0, 1.5731556614014, 0, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## (if using JTransforms) paste(exp(ans.fftR.32.dt002.both[[1]]) * ans.fftR.32.dt002.both[[2]][,2], collapse=", ")
## expectedSp1DsAfterPropTandX = D's after propagate t and x = 1.07607462857111, 1.57350796410007, 2.2106689156898, 2.98405395410774, 3.87006132459304, 4.82224125131593, 5.76741044168543, 6.62885082364537, 7.32184914756387, 7.77191831137883, 7.92794195899393, 7.77172522522963, 7.3214854000561, 6.62835697979005, 5.76683776896087, 4.8216430122779, 3.86948229161728, 2.98353343056498, 2.21022855749815, 1.5731556614014, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0



## (9) testIntegrateOneBranchHiRes32BinsInsideClassBothXandTDt002
##
## requires control.fft.32.dt002 from (8)

dx <- 0.01
hi.lo.ratio <- 4
drift <- 0.0
diffusion <- 0.001
sd <- 0.05
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y
lambda <- sigmoid.x # lambda is now a function
mu <- constant.x # mu is now a function

args.fft.32.dt002.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.32.dt002.both <- c(y0, y1, xmid, r, death, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.32.dt002.both <- quasse.extent.debug(control.fft.32.dt002, drift, diffusion) # prepares X axis stuff

## initialization
vars.fft.32.dt002.both <- matrix(0, control.fft.32.dt002$nx * 4, 2) # high resolution

## populating (see previous tests for explanation if necessary at this point)
vars.fft.32.dt002.both[seq_len(ext.fft.32.dt002.both$ndat[1]),2] <- dnorm(ext.fft.32.dt002.both$x[[1]], 0.0, sd)

exp.sp1.initial.ds.dt002 <- rep(0, length(vars.fft.32.dt002.both[,2])*2)
idxs <- seq(1, length(vars.fft.32.dt002.both[,2])*2, by=2)
j <- 1
for (i in vars.fft.32.dt002.both[,2]) {
    exp.sp1.initial.ds.dt002[idxs[j]] = i
    j = j+1
}

print(paste(exp.sp1.initial.ds.dt002[1:20], collapse=", " ))
## expectedInitialDsFirst10 = 0.791000831787404, 0, 0.879671919608544, 0, 0.975840371583655, 0, 1.07981933026376, 0, 1.19189412137632, 0, 1.31231629549353, 0, 1.44129748672436, 0, 1.57900316601788, 0, 1.72554637653023, 0, 1.88098154753774, 0

## (if using JTransforms) print(paste(vars.fft.32.dt002.both[1:10,2], collapse=", "))
## expectedInitialDsFirst10 = 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774

print(paste(exp.sp1.initial.ds.dt002[155:173], collapse=", " ))
## expectedInitialDsFirst10 = 1.88098154753774, 0, 1.72554637653023, 0, 1.57900316601788, 0, 1.44129748672436, 0, 1.31231629549353, 0, 1.19189412137632, 0, 1.07981933026376, 0, 0.975840371583657, 0, 0.879671919608545, 0, 0.791000831787405

## (if using JTransforms)  print(paste(vars.fft.32.dt002.both[78:87,2], collapse=", "))
## expectedInitialDsLater10 = 1.88098154753774, 1.72554637653023, 1.57900316601788, 1.44129748672436, 1.31231629549353, 1.19189412137632, 1.07981933026376, 0.975840371583657, 0.879671919608545, 0.791000831787405

## have to normalize it to match behavior inside likelihood class (4 because it's high resolution)
vars.fft.32.dt002.both[,2] <- vars.fft.32.dt002.both[,2] / (sum(vars.fft.32.dt002.both[,2]) * dx / hi.lo.ratio)

## not tests below, just checking
## print(paste(vars.fft.32.dt002.both[1:10,2], collapse=", "))
## 0.815139677630763, 0.906516726853789, 1.00561993609328, 1.11277200402141, 1.22826696360519, 1.35236404194541, 1.48528133155527, 1.62718935303434, 1.77820459925865, 1.93838316051407

## print(paste(vars.fft.32.dt002.both[78:87,2], collapse=", "))
## 1.93838316051407, 1.77820459925865, 1.62718935303434, 1.48528133155527, 1.35236404194541, 1.22826696360519, 1.11277200402141, 1.00561993609329, 0.90651672685379, 0.815139677630764

pars.fft.32.dt002.both <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.32.dt002.both, ext.fft.32.dt002.both, pars.32.dt002.both)

pde.fftR.32.dt002.both <- with(control.fft.32.dt002, make.pde.quasse.fftR.debug2(32 * 4, dx/4, dt.max, 2L)) # note the 32 * 4 to do high res, as well as the dx/4

ans.fftR.32.dt002.both <- pde.fftR.32.dt002.both(vars.fft.32.dt002.both, control.fft.32.dt002$dt.max, pars.fft.32.dt002.both$hi, 0)

## unnormalizing it from what make.pde.quasse.fftR.debug2 does
## log-normalization factor = ans.fftR.32.dt002.both[[1]], we un-log it
ans.fftR.32.dt002.both[[2]][,2] <- ans.fftR.32.dt002.both[[2]][,2] * exp(ans.fftR.32.dt002.both[[1]])

ans.fftR.32.dt002.both[[2]][,2] <- ans.fftR.32.dt002.both[[2]][,2] / (sum(ans.fftR.32.dt002.both[[2]][,2]) * dx / hi.lo.ratio) # re-normalize it to match the java class

## E's
exp.es.t.and.x.dt002 <- rep(0, length(ans.fftR.32.dt002.both[[2]][,1])*2)
idxs <- seq(1, length(ans.fftR.32.dt002.both[[2]][,1])*2, by=2)
j <- 1
for (i in ans.fftR.32.dt002.both[[2]][,1]) {
    exp.es.t.and.x.dt002[idxs[j]] = i
    j = j+1
}

print(paste(exp.es.t.and.x.dt002[1:20], collapse=", "))
## expectedSp1EsAfterPropTandXFirst10 = 0.000598961614956886, 0, 0.000598960696294835, 0, 0.000598959776885031, 0, 0.000598958856744621, 0, 0.00059895793589067, 0, 0.000598957014340337, 0, 0.000598956092110961, 0, 0.000598955169219663, 0, 0.000598954245683891, 0, 0.000598953321520975, 0

## (if using JTransforms) print(paste(ans.fftR.32.dt002.both[[2]][1:10,1], collapse=", "))
## expectedSp1EsAfterPropTandXFirst10 = first 10 after t and x = 0.000598961614956886, 0.000598960696294835, 0.000598959776885031, 0.000598958856744621, 0.00059895793589067, 0.000598957014340337, 0.000598956092110961, 0.000598955169219663, 0.000598954245683891, 0.000598953321520975

print(paste(exp.es.t.and.x.dt002[155:173], collapse=", " ))
## expectedSp1EsAfterPropTandXLater10 = 0.000598889987793744, 0, 0.000598889063763026, 0, 0.000598888140363113, 0, 0.000598887217611451, 0, 0.000598886295525308, 0, 0.000598885374121974, 0, 0.000598884453418597, 0, 0.0005988835334325, 0, 0.000598882614180599, 0, 0.00059888169568006

## (if using JTransforms) print(paste(ans.fftR.32.dt002.both[[2]][78:87,1], collapse=", "))
## expectedSp1EsAfterPropTandXLater10 = first 10 after t and x = 0.000598889987793744, 0.000598889063763026, 0.000598888140363113, 0.000598887217611451, 0.000598886295525308, 0.000598885374121974, 0.000598884453418597, 0.0005988835334325, 0.000598882614180599, 0.00059888169568006

## D's
exp.ds.t.and.x.dt002 <- rep(0, length(ans.fftR.32.dt002.both[[2]][,2])*2)
idxs <- seq(1, length(ans.fftR.32.dt002.both[[2]][,2])*2, by=2)
j <- 1
for (i in ans.fftR.32.dt002.both[[2]][,2]) {
    exp.ds.t.and.x.dt002[idxs[j]] = i
    j = j+1
}

print(paste(exp.ds.t.and.x.dt002[1:20], collapse=", "))
## expectedSp1DsAfterPropTandXFirst10 = 0.816826443231728, 0, 0.908389791093209, 0, 1.00769467468308, 0, 1.11506438544145, 0, 1.23079348607031, 0, 1.35514165820029, 0, 1.48832736122564, 0, 1.63052138243079, 0, 1.78184036871647, 0, 1.9423404395558, 0

print(paste(exp.ds.t.and.x.dt002[155:173], collapse=", "))
## expectedSp1DsAfterPropTandXLater10 = 1.94192952805763, 0, 1.78145241269027, 0, 1.63015631463554, 0, 1.4879849556747, 0, 1.35482154601338, 0, 1.23049517164847, 0, 1.11478726270473, 0, 1.00743804313279, 0, 0.908152871485605, 0, 0.816608392676392

## (if using JTransforms) print(paste(ans.fftR.32.dt002.both[[2]][1:10,2], collapse=", "))
## expectedSp1DsAfterPropTandXFirst10 = 0.816826443231728, 0.908389791093209, 1.00769467468308, 1.11506438544145, 1.23079348607031, 1.35514165820029, 1.48832736122564, 1.63052138243079, 1.78184036871647, 1.9423404395558

## (if using JTransforms) print(paste(ans.fftR.32.dt002.both[[2]][78:87,2], collapse=", "))
## expectedSp1DsAfterPropTandXLater10 = 1.94192952805763, 1.78145241269027, 1.63015631463554, 1.4879849556747, 1.35482154601338, 1.23049517164847, 1.11478726270473, 1.00743804313279, 0.908152871485605, 0.816608392676392



## Run for 10, 11, 12 and 13

tc <- 0.005
dt <- 0.005
dx <- 0.01
w <- 10
drift <- 0.0
diffusion <- 0.001
sd <- 0.05
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y

lambda <- function(x) sigmoid.x(x, y0, y1, xmid, r)
mu <- function(x) constant.x(x, death)
## char <- make.brownian.with.drift(0, 0.025) # not sure why 0.025
char <- make.brownian.with.drift(drift, diffusion) # not sure why 0.025

set.seed(1)
tr <- tree.quasse(c(lambda, mu, char), max.taxa=2, x0=0, single.lineage=FALSE, verbose=TRUE)
# simplifying tree for unit test
tr$tip.state[1] <- 0.0 # ch state
tr$tip.state[2] <- 0.1 # ch state
tr$edge.length <- c(0.01, 0.01)  # branch lengths
tr$orig[,2] <- c(0.01, 0.01) # branch lengths
tr$orig[,4] <- c(0.0, 0.1) # ch states

pars <- c(.1, .2, 0, 2.5, .03, drift, diffusion)

control.C.1 <- list(tc=0.005,
                    dt.max=0.005,
                    nx=32,
                    dx=0.01,
                    r=4L,
                    xmid=0,
                    w=10,
                    flags=0L,
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC") # single dt covering entire bifurcating tree height

control.R.1 <- list(tc=0.005,
                    dt.max=0.005,
                    nx=32,
                    dx=0.01,
                    r=4L,
                    xmid=0,
                    w=10,
                    flags=0L,
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftR") # single dt covering entire bifurcating tree height

## reverse-engineering code

## lik.C.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
## (ll.C.1 <- lik.C.1(pars)) ## -6.389708

lik.C.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1) # for some reason, this matches the R result

lik.R.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)

## prints it all on screen
(ll.C.1 <- lik.C.1(pars)) ## -6.389642, differs from make.quasse, within diversitree (see above)
(ll.R.1 <- lik.R.1(pars)) ## -6.389642

## output of above command, from make.quasse.debug
## for JTransforms
root.ds.consecutive <- c(4.66345926030482e-05, 0.000175320546519251, 0.000608432288010798, 0.00197994686277693, 0.00584255914263498, 0.0159175108393565, 0.0400377135639036, 0.0929790087002835, 0.199352433863763, 0.394619574541096, 0.72119984930653, 1.21689198263463, 1.89568999042485, 2.72647162610174, 3.62036642668937, 4.43834919455334, 5.02351058794779, 5.24939946304328, 5.06439897859485, 4.51088490297466, 3.70945996798424, 2.81627005522746, 1.97008358430641, 1.27377969405492, 0.760235257721598)

root.ds.interdig <- rep(0, (length(root.ds.consecutive)+14)*2)
idxs <- seq(1, length(root.ds.consecutive)*2, by=2)
j <- 1
for (i in root.ds.consecutive) {
    root.ds.interdig[idxs[j]] = i
    j = j+1
}

print(paste(root.ds.interdig, collapse=", "))
## input for Java unit test
## rootDs = 4.66345926030482e-05, 0, 0.000175320546519251, 0, 0.000608432288010798, 0, 0.00197994686277693, 0, 0.00584255914263498, 0, 0.0159175108393565, 0, 0.0400377135639036, 0, 0.0929790087002835, 0, 0.199352433863763, 0, 0.394619574541096, 0, 0.72119984930653, 0, 1.21689198263463, 0, 1.89568999042485, 0, 2.72647162610174, 0, 3.62036642668937, 0, 4.43834919455334, 0, 5.02351058794779, 0, 5.24939946304328, 0, 5.06439897859485, 0, 4.51088490297466, 0, 3.70945996798424, 0, 2.81627005522746, 0, 1.97008358430641, 0, 1.27377969405492, 0, 0.760235257721598, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

root.es.consecutive <- c(0.000299741357883797, 0.000299740440744553, 0.000299739520470931, 0.00029973859731951, 0.000299737671601055, 0.000299736743577048, 0.000299735813529964, 0.000299734881744545, 0.000299733948508, 0.000299733014109106, 0.000299732078838168, 0.000299731142986417, 0.000299730206846241, 0.00029972927070989, 0.000299728334869967, 0.000299727399618395, 0.000299726465246688, 0.000299725532045103, 0.000299724600302324, 0.000299723670305445, 0.000299722742339267, 0.000299721816686102, 0.00029972089360741, 0.000299719973413787, 0.000299719056361797)

root.es.interdig <- rep(0, (length(root.es.consecutive)+14)*2)
idxs <- seq(1, length(root.es.consecutive)*2, by=2)
j <- 1
for (i in root.es.consecutive) {
    root.es.interdig[idxs[j]] = i
    j = j+1
}

print(paste(root.es.interdig, collapse=", "))
## rootEs = 0.000299741357883797, 0, 0.000299740440744553, 0, 0.000299739520470931, 0, 0.00029973859731951, 0, 0.000299737671601055, 0, 0.000299736743577048, 0, 0.000299735813529964, 0, 0.000299734881744545, 0, 0.000299733948508, 0, 0.000299733014109106, 0, 0.000299732078838168, 0, 0.000299731142986417, 0, 0.000299730206846241, 0, 0.00029972927070989, 0, 0.000299728334869967, 0, 0.000299727399618395, 0, 0.000299726465246688, 0, 0.000299725532045103, 0, 0.000299724600302324, 0, 0.000299723670305445, 0, 0.000299722742339267, 0, 0.000299721816686102, 0, 0.00029972089360741, 0, 0.000299719973413787, 0, 0.000299719056361797, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0


## (10) testPriorProbAtRootObserved
root.p <- c(0.000101936764975178, 0.000383226449472287, 0.0013299487715952, 0.00432788980772701, 0.0127710256471166, 0.0347934756336926, 0.0875169004358447, 0.203239244270493, 0.435756829076961, 0.862583772672889, 1.57644305300714, 2.65995745025329, 4.14371594632748, 5.95968961768647, 7.9136199323705, 9.70161705011622, 10.9806988668854, 11.474460683722, 11.0700752296872, 9.86017006942937, 8.10836608266039, 6.15597655521472, 4.30633004611629, 2.7843061138813, 1.66176905311156)

print(paste(root.p, collapse=", "))
## expectedObservedPriorProbAtRoot = 0.000101936764975178, 0.000383226449472287, 0.0013299487715952, 0.00432788980772701, 0.0127710256471166, 0.0347934756336926, 0.0875169004358447, 0.203239244270493, 0.435756829076961, 0.862583772672889, 1.57644305300714, 2.65995745025329, 4.14371594632748, 5.95968961768647, 7.9136199323705, 9.70161705011622, 10.9806988668854, 11.474460683722, 11.0700752296872, 9.86017006942937, 8.10836608266039, 6.15597655521472, 4.30633004611629, 2.7843061138813, 1.66176905311156



## (11) testPriorProbAtRootFlat
##
## depends on [10], as it uses exp.root.p

(ll.C.1 <- lik.C.1(pars, root=ROOT.FLAT))

root.p.flat <- 3.2258064516129

exp.root.p.flat <- rep(root.p.flat, length(root.p)) ## we do not interdigitate the root prior, even when using JavaFftService

print(paste(exp.root.p.flat, collapse=", "))
## expectedFlatPriorProbAtRoot = 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129, 3.2258064516129



## (12) testRootCalcProcedure
(ll.R.1 <- lik.R.1(pars)) ## -6.389642

## logLik = -6.389642



## (13) testPruneBifTree32Bins

ds.initial.sp1 <- c(0.308986942687903, 0.350566009871371, 0.396747087835907, 0.447890605896858, 0.504364398303888,
                    0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544,
                    0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436,
                    1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911,
                    2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839,
                    3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101,
                    4.59764281368466, 4.83941449038287, 5.08118112938378, 5.3217049979751, 5.55969772261993,
                    5.79383105522966, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599,
                    6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876,
                    7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024,
                    7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778,
                    7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647,
                    7.21053924923296, 7.04130653528599, 6.85887710038768, 6.66449205783599, 6.45944719335828,
                    6.24507866733522, 6.02274864309609, 5.79383105522965, 5.55969772261993, 5.32170499797509,
                    5.08118112938378, 4.83941449038287, 4.59764281368466, 4.35704354065101, 4.11872537439949,
                    3.88372109966426, 3.65298170778044, 3.42737184095615, 3.20766654683839, 2.9945493127149,
                    2.7886113289072, 2.59035191331783, 2.40018001393971, 2.21841669358911, 2.04529849127956,
                    1.88098154753774, 1.72554637653023, 1.57900316601788, 1.44129748672436, 1.31231629549353,
                    1.19189412137632, 1.07981933026376, 0.975840371583655, 0.879671919608544, 0.791000831787404,
                    0.709491856924628, 0.634793036713349, 0.566540754832024, 0.504364398303888, 0.447890605896858,
                    0.396747087835907, 0.350566009871371, 0.308986942687903)

exp.ds.initial.sp1 <- rep(0, (length(ds.initial.sp1)+25)*2) ## note the +25 to add up to 128
idxs <- seq(1, (length(ds.initial.sp1)+25)*2, by=2)
j <- 1
for (i in ds.initial.sp1) {
    exp.ds.initial.sp1[idxs[j]] = i
    j = j+1
}

print(paste(exp.ds.initial.sp1, collapse=", "))
## expectedDsHiAtNodeInitialSp1 = 0.308986942687903, 0, 0.350566009871371, 0, 0.396747087835907, 0, 0.447890605896858, 0, 0.504364398303888, 0, 0.566540754832024, 0, 0.634793036713348, 0, 0.709491856924629, 0, 0.791000831787404, 0, 0.879671919608544, 0, 0.975840371583655, 0, 1.07981933026376, 0, 1.19189412137632, 0, 1.31231629549353, 0, 1.44129748672436, 0, 1.57900316601788, 0, 1.72554637653023, 0, 1.88098154753774, 0, 2.04529849127956, 0, 2.21841669358911, 0, 2.40018001393971, 0, 2.59035191331783, 0, 2.7886113289072, 0, 2.9945493127149, 0, 3.20766654683839, 0, 3.42737184095615, 0, 3.65298170778044, 0, 3.88372109966426, 0, 4.11872537439949, 0, 4.35704354065101, 0, 4.59764281368466, 0, 4.83941449038287, 0, 5.08118112938378, 0, 5.3217049979751, 0, 5.55969772261993, 0, 5.79383105522966, 0, 6.02274864309609, 0, 6.24507866733522, 0, 6.45944719335828, 0, 6.66449205783599, 0, 6.85887710038768, 0, 7.04130653528599, 0, 7.21053924923296, 0, 7.36540280606647, 0, 7.50480693833876, 0, 7.62775630921048, 0, 7.73336233605698, 0, 7.82085387950912, 0, 7.88958661815778, 0, 7.93905094954024, 0, 7.96887828189528, 0, 7.97884560802865, 0, 7.96887828189528, 0, 7.93905094954024, 0, 7.88958661815778, 0, 7.82085387950912, 0, 7.73336233605698, 0, 7.62775630921048, 0, 7.50480693833876, 0, 7.36540280606647, 0, 7.21053924923296, 0, 7.04130653528599, 0, 6.85887710038768, 0, 6.66449205783599, 0, 6.45944719335828, 0, 6.24507866733522, 0, 6.02274864309609, 0, 5.79383105522965, 0, 5.55969772261993, 0, 5.32170499797509, 0, 5.08118112938378, 0, 4.83941449038287, 0, 4.59764281368466, 0, 4.35704354065101, 0, 4.11872537439949, 0, 3.88372109966426, 0, 3.65298170778044, 0, 3.42737184095615, 0, 3.20766654683839, 0, 2.9945493127149, 0, 2.7886113289072, 0, 2.59035191331783, 0, 2.40018001393971, 0, 2.21841669358911, 0, 2.04529849127956, 0, 1.88098154753774, 0, 1.72554637653023, 0, 1.57900316601788, 0, 1.44129748672436, 0, 1.31231629549353, 0, 1.19189412137632, 0, 1.07981933026376, 0, 0.975840371583655, 0, 0.879671919608544, 0, 0.791000831787404, 0, 0.709491856924628, 0, 0.634793036713349, 0, 0.566540754832024, 0, 0.504364398303888, 0, 0.447890605896858, 0, 0.396747087835907, 0, 0.350566009871371, 0, 0.308986942687903, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

ds.initial.sp2 <- c(0.000254946647636669, 0.000319674822138109, 0.000399835934138456, 0.000498849425801072, 0.000620828141157003,
                    0.000770703934841743, 0.000954372730824099, 0.0011788613551308, 0.00145251860604505, 0.00178523314354266,
                    0.00218868086879601, 0.00267660451529771, 0.00326512817532484, 0.00397310942785545, 0.00482253160451986,
                    0.0058389385158292, 0.00705191364734891, 0.00849560541101504, 0.0102092994868837, 0.0122380386022755,
                    0.0146332892566062, 0.0174536539009152, 0.0207656259132282, 0.0246443833694604, 0.0291746160933349,
                    0.0344513787810736, 0.0405809611459954, 0.0476817640292969, 0.0558851682975889, 0.0653363811239983,
                    0.0761952419644361, 0.08863696823876, 0.102852818461079, 0.119050648395517, 0.137455333812279,
                    0.158309031659599, 0.181871250031821, 0.208418696288452, 0.238244872152104, 0.271659384673712,
                    0.308986942687903, 0.350566009871371, 0.396747087835906, 0.447890605896858, 0.504364398303888,
                    0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544,
                    0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436,
                    1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911,
                    2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839,
                    3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101,
                    4.59764281368466, 4.83941449038287, 5.08118112938378, 5.32170499797509, 5.55969772261993,
                    5.79383105522965, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599,
                    6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876,
                    7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024,
                    7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778,
                    7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647,
                    7.21053924923296, 7.04130653528599, 6.85887710038768)

exp.ds.initial.sp2 <- rep(0, (length(ds.initial.sp2)+25)*2) ## note the +25 to add up to 128
idxs <- seq(1, (length(ds.initial.sp2)+25)*2, by=2)
j <- 1
for (i in ds.initial.sp2) {
    exp.ds.initial.sp2[idxs[j]] = i
    j = j+1
}

print(paste(exp.ds.initial.sp2, collapse=", "))
## expectedDsHiAtNodeInitialSp2 = 0.000254946647636669, 0, 0.000319674822138109, 0, 0.000399835934138456, 0, 0.000498849425801072, 0, 0.000620828141157003, 0, 0.000770703934841743, 0, 0.000954372730824099, 0, 0.0011788613551308, 0, 0.00145251860604505, 0, 0.00178523314354266, 0, 0.00218868086879601, 0, 0.00267660451529771, 0, 0.00326512817532484, 0, 0.00397310942785545, 0, 0.00482253160451986, 0, 0.0058389385158292, 0, 0.00705191364734891, 0, 0.00849560541101504, 0, 0.0102092994868837, 0, 0.0122380386022755, 0, 0.0146332892566062, 0, 0.0174536539009152, 0, 0.0207656259132282, 0, 0.0246443833694604, 0, 0.0291746160933349, 0, 0.0344513787810736, 0, 0.0405809611459954, 0, 0.0476817640292969, 0, 0.0558851682975889, 0, 0.0653363811239983, 0, 0.0761952419644361, 0, 0.08863696823876, 0, 0.102852818461079, 0, 0.119050648395517, 0, 0.137455333812279, 0, 0.158309031659599, 0, 0.181871250031821, 0, 0.208418696288452, 0, 0.238244872152104, 0, 0.271659384673712, 0, 0.308986942687903, 0, 0.350566009871371, 0, 0.396747087835906, 0, 0.447890605896858, 0, 0.504364398303888, 0, 0.566540754832024, 0, 0.634793036713348, 0, 0.709491856924629, 0, 0.791000831787404, 0, 0.879671919608544, 0, 0.975840371583655, 0, 1.07981933026376, 0, 1.19189412137632, 0, 1.31231629549353, 0, 1.44129748672436, 0, 1.57900316601788, 0, 1.72554637653023, 0, 1.88098154753774, 0, 2.04529849127956, 0, 2.21841669358911, 0, 2.40018001393971, 0, 2.59035191331783, 0, 2.7886113289072, 0, 2.9945493127149, 0, 3.20766654683839, 0, 3.42737184095615, 0, 3.65298170778044, 0, 3.88372109966426, 0, 4.11872537439949, 0, 4.35704354065101, 0, 4.59764281368466, 0, 4.83941449038287, 0, 5.08118112938378, 0, 5.32170499797509, 0, 5.55969772261993, 0, 5.79383105522965, 0, 6.02274864309609, 0, 6.24507866733522, 0, 6.45944719335828, 0, 6.66449205783599, 0, 6.85887710038768, 0, 7.04130653528599, 0, 7.21053924923296, 0, 7.36540280606647, 0, 7.50480693833876, 0, 7.62775630921048, 0, 7.73336233605698, 0, 7.82085387950912, 0, 7.88958661815778, 0, 7.93905094954024, 0, 7.96887828189528, 0, 7.97884560802865, 0, 7.96887828189528, 0, 7.93905094954024, 0, 7.88958661815778, 0, 7.82085387950912, 0, 7.73336233605698, 0, 7.62775630921048, 0, 7.50480693833876, 0, 7.36540280606647, 0, 7.21053924923296, 0, 7.04130653528599, 0, 6.85887710038768, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

paste(seq(6,(2*(length(ds.initial.sp1)-2)),by=8), collapse=", ")
## expectedHiLoIdxs4Transfer = 6, 14, 22, 30, 38, 46, 54, 62, 70, 78, 86, 94, 102, 110, 118, 126, 134, 142, 150, 158, 166, 174, 182, 190, 198

## (if using JTransforms) expectedHiLoIdxs4Transfer = 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99

lik.C.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
(ll.C.1 <- lik.C.1(pars, root=ROOT.FLAT)) ## -6.389708



## (14) testPruneBifTree0025Height32Bins

tr$edge.length <- c(0.025, 0.025)  # branch lengths
tr$orig[,2] <- c(0.025, 0.025) # branch lengths

lik.C.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1) # for some reason, this matches the R result

lik.R.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)

## prints it all on screen
(ll.C.1 <- lik.C.1(pars)) ## -6.394299 differs from lik.R.1 below
(ll.R.1 <- lik.R.1(pars)) ## -6.389642

## logLik = -6.394235



## (15) testPruneTrifTree32Bins

nx <- 32
tc <- 0.005
dt <- 0.005
dx <- 0.01
drift <- 0.0
diffusion <- 0.001
sd <- 0.05
death <- 0.03 # for constant link function
r <- 2.5 # logistic link function
xmid <- 0.0
y0 <- 0.1 # base y ("y" shift)
y1 <- 0.2 # max y

lambda <- function(x) sigmoid.x(x, y0, y1, xmid, r)
mu <- function(x) constant.x(x, death)
char <- make.brownian.with.drift(drift, 0.025) # not sure why 0.025

set.seed(1)
tr <- tree.quasse(c(lambda, mu, char), max.taxa=3, x0=0, single.lineage=FALSE, verbose=TRUE)
# simplifying tree for unit test
tr$tip.state[1] <- 0.0 # ch state
tr$tip.state[2] <- 0.1 # ch state
tr$tip.state[3] <- 0.2 # ch state
tr$edge.length <- c(0.02, 0.01, 0.01, 0.01)  # branch lengths
tr$orig[,2] <- c(0.02, 0.01, 0.01, 0.01) # branch lengths
tr$orig[,4] <- c(0.0, NA, 0.1, 0.2) # ch states

pars <- c(y0, y1, xmid, r, death, drift, diffusion)

control.C.1 <- list(tc=tc,
                    dt.max=dt,
                    nx=32,
                    dx=dx,
                    r=4L,
                    xmid=0,
                    w=10,
                    flags=0L,
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC") # single dt covering entire bifurcating tree height

control.R.1 <- list(tc=tc,
                    dt.max=dt,
                    nx=32,
                    dx=dx,
                    r=4L,
                    xmid=0,
                    w=10,
                    flags=0L,
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftR") # single dt covering entire bifurcating tree height

## reverse-engineering code
## lik.C.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
lik.C.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
lik.R.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)
(ll.C.1 <- lik.C.1(pars)) ## -9.082284 differs from lik.R.1 below
(ll.C.1 <- lik.R.1(pars)) ## -9.085542 DIFFERENT WITHIN DIVERSITREE

## logLik = -9.085542


## Some random notes...
##
## In R, among the debugging prints
## (1) First normalization from initial D's: log-normalization factor starts at -0.332160418518853
## (2) Second normalization after integrating branch at (first) high-res part: added -0.000910605086357751 to log-normalization factor (SLIGHTLY DIFFERENT FROM JAVA)
## (3) We add them up: log-normalization factor is now -0.333071023605211 (SLIGHTLY DIFFERENT FROM JAVA))
## (4) Third normalization after integrating branch at (second) low-res part: log-normalization factor from integration only = -0.0372438735770201 (DLIGHTLY DIFFERENT FROM JAVA)
## (5) We add them up: final log-normalization factor at return = -0.370314897182231 (SLIGHTLY DIFFERENT FROM JAVA)

## In Java
## for sp2:
## (1) First normalization from initial D's: Log-normalization factor inside normalizeDs = -0.332160418518852
## (2) Second normalization after integrating branch at (first) high-res part: Log-normalization factor inside normalizeDs = -9.236533466328776E-4 (-0.0009236533466328776)
## (3) We add them up: Log-normalization factors after high-res part = -0.33308407186548555 ( = -0.332160418518852 + -9.236533466328776E-4)
## (4) Third normalization after integrating branch at (second) low-res part: -0.037287301482475846
## (5) We add them up: logNormalizationFactors = -0.3703713733479614 ( = -0.33308407186548555 + -0.037287301482475846 )
## this final logNormalizationFactor above from (5) is the second entry of the line below, one of the final prints of the BEAST 2 pruning debugging messages
## logNormalizationFactors = [-0.015948966160484787, -0.3703713733479614, -2.712667833864449, 0.5214770496715425, 0.0]



## (16) testPruneFifteenSpTree1024Bins

y0 <- 0.1
y1 <- 0.2
xmid <- 0
r <- 2.5
death <- 0.03
drift <- 0.0
diffusion <- 0.01

lambda <- function(x) sigmoid.x(x, y0, y1, xmid, r)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(drift, 0.025) ## second argument is not diffusion, but in fact leads to a Normal kernel sd of 0.00707106781186548, which matches Java (not sure how I arrived at 0.025)

set.seed(1)
tr <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0, single.lineage=FALSE, verbose=TRUE)

nodes <- c("nd13", "nd9", "nd5")
split.t <- Inf

pars <- c(y0, y1, xmid, r, death, drift, .01) # lambda y0, lambda y1, lambda xmid, lambda r, mu, drift, diffusion
pars4 <- unlist(rep(list(pars), 4))

sd <- 1/50

control.C.1 <- list(dt.max=1/200, method="fftC")
lik.C.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
(ll.C.1 <- lik.C.1(pars4)) # -61.27245

## logLik = -61.27245

control.R.1 <- list(dt.max=1/200, method="fftR")
lik.R.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)
(ll.R.1 <- lik.R.1(pars)) # -61.27245
