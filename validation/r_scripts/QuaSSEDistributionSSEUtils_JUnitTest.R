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
## (8) testIntegrateOneBranchLoRes32BinsOutsideClassJustT
## (9) testIntegrateOneBranchLoRes32BinsOutsideClassJustX (see test (4))
## (10) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX
## (11) testIntegrateOneBranchHiRes4096BinsOutsideClassJustX (note that I'm using this in 4096 low res, but I'll co-opt 4096 high res in Java)
## (12) testBothXandTPropagateMethodsInsideClassOneBranchLoRes32Bins
## (13) testBothXandTPropagateMethodsInsideClassOneBranchLoRes32BinsDt002
## (14) testIntegrateOneBranchHiRes32BinsInsideClassBothXandTDt002
## (15) testPriorProbAtRootObserved
## (16) testPriorProbAtRootFlat
## (17) testRootCalcProcedure
## (18) testPruneBifTree32Bins
## (19) testPruneTrifTree32Bins
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
## (1) testPropagateTimeOneChQuaSSETest

vars <- cbind(rep(0.0001,32),
              c(8.50976807665641, 10.103792974434, 9.08976088347418, 11.847721337896, 8.51254745547751, 9.91650581983555, 8.95019918832521, 9.30609468137578, 11.0775496365384, 10.7639029400606, 10.9164931483932, 9.83064984974005, 11.7045125626528, 11.3382431919839, 8.94185500956388, 7.30298759647754, 11.1167065386435, 9.76891399488789, 9.76676261926709, 9.10540040702707, 8.93655752085786, 10.2580116547857, 10.2552822093573, 8.85921172559191, 11.0314684537514, 10.8738197102994, 10.4638936963999, 9.68617874031991, 8.35885856359494, 10.8426829704597, 7.66894489549493, 8.23694434625264))
# vars <- data.frame(0.0, 1.0) # E and D
lambda <- 1.0
mu <- 0.5
dt <- 0.01
ndat <- 28 # just one useful quant ch bin

#### eAndD

res <- fftR.propagate.t.debug(vars, lambda, mu, dt, ndat)
print(paste(res[1:10,1], collapse=", "))
print(paste(res[1:10,2], collapse=", "))

##              [,1]      [,2]
##  [1,] 0.005061285  8.383508
##  [2,] 0.005061285  9.953882
##  [3,] 0.005061285  8.954895
##  [4,] 0.005061285 11.671936
##  [5,] 0.005061285  8.386246
##  [6,] 0.005061285  9.769374
##  [7,] 0.005061285  8.817404
##  [8,] 0.005061285  9.168019
##  [9,] 0.005061285 10.913191
## [10,] 0.005061285 10.604198



## (2) testMakeNormalKernInPlaceAndFFT
##
## The normal kernel (fy) gives the probability of a change in quantitative
## trait values; here, each x bin is the CHANGE in quantitative trait values
## (not the quantitative trait values themselves!)
##
## Under BM, for example, this kernel is centered at 0.0
##
## This kernel (yValues) is later FFT-ed in the likelihood class, and then used
## in the convolution function, and for reasons I do not fully understand,
## this bell-shaped kernel needs to be cuts in two, with the left half being
## placed at the (right-)tail end of the kernel, and the right half at the
## (left-)head of the kernel

drift <- 0.0
diffusion <- 0.001
nkl <- nkr <- 2
nx <- 32

# changed from 0.001 -> 0.01, updating!
dx <- 0.01 # (dx * nx) gives the range of x values we're considering

# '-' in -dt because we're going back in time
kern <- fftR.make.kern.debug(-dt * drift, sqrt(dt * diffusion), nx, dx, nkl, nkr)

paste(kern, collapse=", ") # expectedFy = 0.986703287028858, 0.00664835445182386, 2.03374705433156e-09, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.03374705433156e-09, 0.00664835445182386

## plot(kern)

## x: the substitution rate
## y: gives you the probability of changing from x to (x + delta.x) over time dt,
## with delta.x being (dx * nkr) on the right-hand side of x, and (dx * nkl) on
## the left-hand side of x
##
## (later, y will be fourier-transformed into fy)

## worked example
##
## min x = 0.0
## max x = 0.1 = 1/10
## start x at rate = 0.01 = 1/100
## dx = 0.001
## nx = 100
## nkl = nkr = 4
##
## x is at 10th bin, so we ignore anything beyond nkl and nkr bins away from 10th, meaning
## x will influence from 5th - 15th bins
##
## the first bin gives the probability that x stays where it is
## the 2nd - 6th bin give the probability x changes to 11th - 15th bins
## the 95th - 100th bins give the probability x changes to the 9th - 5th

fy <- fft(kern)

## We don't immediately ifft fy in the likelihood calculation, this is being
## done here for testing purposes
ify <- ifft(fy)

paste(Re(fy)[1:10], collapse=", ") # expectedFftFY = 1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216

paste(Re(ify)[1:10], collapse=", ") # expectedIfftFY = 31.5745051849235, 0.212747342458364, 6.50799056989047e-08, 1.52655665885959e-16, -4.44089209850063e-16, -2.35784152995104e-16, -4.15496150260171e-16, -2.91433543964104e-16, -8.88178419700125e-16, 8.15320033709099e-17



## (3) testConvolve

## see vars.fft.just.x below to see where this came from
vars <- cbind(rep(0.0, 32), c(0.0058389385158292, 0.0122380386022755, 0.0246443833694604, 0.0476817640292969, 0.0886369682387602, 0.158309031659599, 0.271659384673712, 0.447890605896858, 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 0, 0, 0, 0))


## convolve
paste(Re(apply(apply(vars, 2, fft) * fy, 2, ifft))[,2], collapse=", ") # expectedDsUnnormalized = 0.186965206235882, 0.392895262124355, 0.790881995397243, 1.52962842317328, 2.84249242686772, 5.0751814647876, 8.70647806374478, 14.3506616579459, 22.726870664534, 34.581632450077, 50.5579348423759, 71.0184208525211, 95.849626817365, 124.293227499606, 154.860992049653, 185.384774972034, 213.227631891732, 235.64066883828, 250.204040383914, 255.255834736325, 250.204040383914, 235.64066883828, 213.227631891732, 185.384774972034, 154.860992049653, 124.293227499606, 95.8496267146036, 70.6824920549429, 0.471962450911363, 1.44374375565803e-07, 3.80026676793932e-10, 0.00124221944851399

## same as above, but normalizing Es and Ds by the number of bins
vars.out <- Re(apply(apply(vars, 2, fft) * fy, 2, ifft)) / nx
paste(vars.out[,1], collapse=", ") # expectedEs = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(vars.out[,2], collapse=", ")
## expectedDs = 0.00584266269487133, 0.0122779769413861,
##
## 0.0247150623561638, 0.047800888224165, 0.0888278883396163, 0.158599420774613, 0.272077439492024, 0.448458176810811, 0.710214708266688, 1.0806760140649, 1.57993546382425, 2.21932565164129, 2.99530083804266, 3.8841633593627, 4.83940600155166, 5.79327421787606, 6.66336349661661, 7.36377090119626, 7.81887626199732, 7.97674483551017, 7.81887626199732, 7.36377090119626, 6.66336349661661, 5.79327421787606,
##
## (continuation of expectedDs)
## 4.83940600155166, 3.88416335936269, 2.99530083483136, 2.20882787671697, 0.0147488265909801, 4.51169923643135e-09, 1.18758336498104e-11, 3.88193577660623e-05



## (4) testPropagateChOneCh32QuaSSETest
##
## This test differs from (3) in that fftR.propagate.x.debug will manage flanking bins
## (ignoring the leftmost and rightmost 2 bins), as it should be done during likelihood
## calculation. Note how the middle bins match (3), though.
##
## nx = 32
## nkl and nkr = 2
res <- fftR.propagate.x.debug(vars, nx, fy, nkl, nkr)

paste(res[,1], collapse=", ") # 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(res[,2], collapse=", ")
## expectedDs = 0.0058389385158292, 0.0122380386022755,
##
## 0.0247150623561638, 0.047800888224165, 0.0888278883396163, 0.158599420774613, 0.272077439492024, 0.448458176810811, 0.710214708266688, 1.0806760140649, 1.57993546382425, 2.21932565164129, 2.99530083804266, 3.8841633593627, 4.83940600155166, 5.79327421787606, 6.66336349661661, 7.36377090119626, 7.81887626199732, 7.97674483551017, 7.81887626199732, 7.36377090119626, 6.66336349661661, 5.79327421787606,
##
## (continuation of expectedDs)
## 4.83940600155166, 3.88372109966426, 2.9945493127149, 0, 0, 0, 0, 0



# (5) testLogistic

## The following code is inside quasse.extent
drift <- 0.0; diffusion <- 0.01
nx <- 1024; dx <- 0.01; dt <- 0.05
y0 <- 0.1; y1 <- 0.2; xmid <- 0; r <- 2.5; death <- 0.03
hi.lo.ratio <- 4 # r=4L when making control.fft
w <- 5

mean4test <- drift * dt
sd4test <- sqrt(diffusion * dt)
nkleft <- max(ceiling(-(mean4test - w * sd4test)/dx)) * c(hi.lo.ratio, 1)
nkright <- max(ceiling((mean4test + w * sd4test)/dx)) * c(hi.lo.ratio, 1)

ndat <- nx * c(hi.lo.ratio, 1) - (nkleft + 1 + nkright)
ndat.lo <- ndat[2];
ndat.hi <- ndat[1];

xmin.lo <- xmid - dx * ceiling((ndat.lo - 1)/2) # x.02
xmin.hi <- xmin.lo - dx * (1 - 1/hi.lo.ratio) # x.01
x.lo <- seq(xmin.lo, length.out=ndat.lo, by = dx) # same as ext.fft$x[[2]]
x.hi <- seq(xmin.hi, length.out=ndat.hi, by = dx/hi.lo.ratio) # same as ext.fft$x[[1]]

ls.hi <- sigmoid.x(x.hi, y0, y1, xmid, r) # same as pars.fft$hi$lambda
ls.lo <- sigmoid.x(x.lo, y0, y1, xmid, r) # same as pars.fft$hi$lambda

paste(ls.hi[1:10], collapse=", ")
# "0.100000375000363, 0.100000377351446, 0.100000379717269, 0.100000382097925, 0.100000384493506, 0.100000386904106, 0.10000038932982, 0.100000391770742, 0.100000394226967, 0.100000396698591"

paste(ls.hi[2501:2510], collapse=", ")
# "0.195816352937447, 0.195841335181637, 0.195866174682834, 0.195890872179954, 0.195915428409005, 0.195939844103087, 0.19596411999239, 0.195988256804195, 0.196012255262877, 0.196036116089903"

paste(ls.hi[3989:3999], collapse=", ")
# "0.199999600814288, 0.199999603301409, 0.199999605773033, 0.199999608229258, 0.19999961067018, 0.199999613095894, 0.199999615506494, 0.199999617902075, 0.199999620282731, 0.199999622648554, 0.199999624999637"



# (6) testDimensions

## Requires some variables from (5)
# modifying some of them though...
dt <- 0.01; dx <- 0.0005; diffusion <- 0.001; w <- 10

sd4test <- sqrt(diffusion * dt);
nkleft <- max(ceiling(-(mean4test - w * sd4test)/dx)) * c(hi.lo.ratio, 1)
nkright <- max(ceiling((mean4test + w * sd4test)/dx)) * c(hi.lo.ratio, 1)

ndat <- nx * c(hi.lo.ratio, 1) - (nkleft + 1 + nkright)
ndat.lo <- ndat[2]; ndat.lo
ndat.hi <- ndat[1]; ndat.hi

xmin.lo <- xmid - dx * ceiling((ndat.lo - 1)/2) # x.02
xmin.hi <- xmin.lo - dx * (1 - 1/hi.lo.ratio) # x.01
x.lo <- seq(xmin.lo, length.out=ndat.lo, by = dx) # same as ext.fft$x[[2]]
x.hi <- seq(xmin.hi, length.out=ndat.hi, by = dx/hi.lo.ratio) # same as ext.fft$x[[1]]

paste(x.lo[1:10], collapse=", ") # expectedXLoFirst10 = -0.2235, -0.223, -0.2225, -0.222, -0.2215, -0.221, -0.2205, -0.22, -0.2195, -0.219
paste(x.lo[886:895], collapse=", ") # expectedXLoLast10 = 0.219, 0.2195, 0.22, 0.2205, 0.221, 0.2215, 0.222, 0.2225, 0.223, 0.2235

paste(x.hi[1:10], collapse=", ") # expectedXHiFirst10 = -0.223875, -0.22375, -0.223625, -0.2235, -0.223375, -0.22325, -0.223125, -0.223, -0.222875, -0.22275
paste(x.hi[3574:3583], collapse=", ") # expectedXHiFirst10 = 0.22275, 0.222875, 0.223, 0.223125, 0.22325, 0.223375, 0.2235, 0.223625, 0.22375, 0.223875

lambda <- sigmoid.x
mu <- constant.x

# lambdas
paste(do.call(lambda, c(list(x.lo), c(y0, y1, xmid, r)))[1:10], collapse=", ") # expectedLambdaLoFirt10 = 0.13638367349015, 0.136412610857295, 0.13644155805569, 0.136470515067719, 0.136499481875741, 0.136528458462084, 0.136557444809052, 0.13658644089892, 0.136615446713936, 0.136644462236322
paste(do.call(lambda, c(list(x.lo), c(y0, y1, xmid, r)))[886:895], collapse=", ") # expectedLambdaLoLast10 = 0.163355537763678, 0.163384553286064, 0.16341355910108, 0.163442555190948, 0.163471541537916, 0.163500518124259, 0.163529484932281, 0.16355844194431, 0.163587389142705, 0.16361632650985

paste(do.call(lambda, c(list(x.hi), c(y0, y1, xmid, r)))[1:10], collapse=", ") # expectedLambdaHiFirt10 = 0.136361976927131, 0.136369208498794, 0.136376440686558, 0.13638367349015, 0.136390906909294, 0.136398140943716, 0.136405375593141, 0.136412610857295, 0.136419846735901, 0.136427083228686
paste(do.call(lambda, c(list(x.hi), c(y0, y1, xmid, r)))[3574:3583], collapse=", ") # expectedLambdaHiLast10 = 0.163572916771314, 0.163580153264099, 0.163587389142705, 0.163594624406859, 0.163601859056284, 0.163609093090706, 0.16361632650985, 0.163623559313442, 0.163630791501206, 0.163638023072869

# mus
paste(do.call(mu, c(list(x.lo), death))[1:10], collapse=", ") # expectedMuLoFirst10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(do.call(mu, c(list(x.lo), death))[886:895], collapse=", ") # expectedMuLoLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03

paste(do.call(mu, c(list(x.hi), death))[1:10], collapse=", ") # expectedMuHiFirst10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03
paste(do.call(mu, c(list(x.hi), death))[3574:3583], collapse=", ") # expectedMuHiLast10 = 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03




# (7) testInitializationOfTips

drift <- 0.0
diffusion <- 0.001
lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)

char <- make.brownian.with.drift(drift, diffusion)
set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa=2, x0=0, single.lineage=FALSE, verbose=FALSE)
pars <- c(.1, .2, 0, 2.5, .03, 0, .001) # 6th and 7th elements are drift and diffusion
sd <- 1/20
control.C.1 <- list(dt.max=0.01) # dt = 0.01
phy$tip.state[1] <- 0.0
phy$tip.state[2] <- 0.1

cache <- diversitree:::make.cache.quasse(phy, phy$tip.state, sd, lambda, mu, control.C.1, NULL)
cache$height[1] <- cache$height[2] <- 0.01 # setting things manually
cache$depth[3] <- 0.01
cache$edge.length <- c(0.01, 0.01)
cache$len <- c(0.01, 0.01, NA)

## these four lines are to make this cache match the previous unit tests
cache$args$drift <- 6
cache$args$diffusion <- 7
cache$control$dx <- 0.0005 ## to match unit test
cache$control$xmid <- 0.0 ## to match unit test
cache$control$w <- 10

## taking following lines from R/model-quasse.R, make.quasse()
# all.branches <- make.all.branches.quasse(cache, cache$control) # all.branches is a function
f.pars <- diversitree:::make.pars.quasse(cache)
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

sp1.y.ds <- sp1.y[nx+1:(length(sp1.y)-nx)]
sp2.y.ds <- sp2.y[nx+1:(length(sp2.y)-nx)]

sp1.y.ds.exp <- sp1.y.ds[2001:2048]
sp2.y.ds.exp <- sp2.y.ds[2001:2048]

paste(sp1.y.ds.exp, collapse=", ") # expectedSp1Ds = 6.96077358436849, 6.95166528584696, 6.94252551477923, 6.93335442671583, 6.92415217758908, 6.91491892370871, 6.90565482175746, 6.89636002878667, 6.88703470221182, 6.87767899980818, 6.86829307970631, 6.85887710038768, 6.84943122068018, 6.83995559975374, 6.83045039711584, 6.82091577260705, 6.81135188639661, 6.80175889897795, 6.79213697116422, 6.78248626408384, 6.772806939176, 6.76309915818623, 6.75336308316186, 6.74359887644761, 6.73380670068105, 6.72398671878815, 6.71413909397875, 6.70426398974212, 6.69436156984244, 6.68443199831428, 6.67447543945815, 6.66449205783599, 6.65448201826663, 6.64444548582134, 6.63438262581928, 6.62429360382306, 6.61417858563415, 6.60403773728847, 6.59387122505179, 6.58367921541529, 6.57346187509105, 6.5632193710075, 6.55295187030495, 6.54265954033109, 6.53234254863645, 6.52200106296994, 6.51163525127429, 6.50124528168164
paste(sp2.y.ds.exp, collapse=", ") # expectedSp2Ds = 2.67859021074856, 2.68849414736158, 2.69841783805887, 2.70836123148143, 2.71832427571076, 2.72830691826807, 2.73830910611356, 2.74833078564564, 2.75837190270027, 2.76843240255022, 2.77851222990441, 2.78861132890721, 2.79872964313779, 2.80886711560949, 2.81902368876922, 2.82919930449678, 2.83939390410431, 2.84960742833573, 2.85983981736613, 2.87009101080125, 2.88036094767694, 2.89064956645866, 2.90095680504094, 2.91128260074695, 2.92162689032799, 2.93198960996305, 2.9423706952584, 2.95277008124712, 2.96318770238874, 2.97362349256884, 2.9840773850987, 2.9945493127149, 3.00503920757903, 3.01554700127736, 3.02607262482052, 3.03661600864323, 3.04717708260405, 3.05775577598507, 3.06835201749176, 3.07896573525268, 3.0895968568193, 3.10024530916586, 3.11091101868917, 3.12159391120842, 3.13229391196515, 3.14301094562307, 3.15374493626797, 3.16449580740766



# (8) testIntegrateOneBranchLoRes32BinsOutsideClassJustT

control.fft.32 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.01, # dt
                    nx=32, # number of X bins
                    dx=0.01, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

lambda <- sigmoid.x
mu <- constant.x
diffusion <- 0.001
sd <- 0.05

args.fft.just.t <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.t <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.t <- quasse.extent.debug(control.fft.32, drift, diffusion) # prepares X axis stuff
## ext.fft.just.t$padding

pars.fft.just.t <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.just.t, ext.fft.just.t, pars.just.t) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.t <- matrix(0, control.fft.32$nx, 2) # low resolution
vars.fft.just.t[seq_len(ext.fft.just.t$ndat[2]),2] <- dnorm(ext.fft.just.t$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x$lo$padding should be 0.0
paste(vars.fft.just.t[,2], collapse=", ") # expectedInitialDs = 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0

vars.fft.just.t.tmp <- fftR.propagate.t.debug(vars.fft.just.t, pars.fft.just.t$lo$lambda, pars.fft.just.t$lo$mu, control.fft.32$dt.max, ext.fft.just.t$ndat[2])
paste(vars.fft.just.t.tmp[,1], collapse=", ") # expectedSp1EsAfterPropT = 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(vars.fft.just.t.tmp[,2], collapse=", ") # expectedSp1DsAfterPropT = 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052118, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0

pde.fftR.just.t <- with(control.fft.32, make.pde.quasse.fftR.debug(nx, dx, dt.max, 2L))
ans.fftR.just.t <- pde.fftR.just.t(vars.fft.just.t, control.fft.32$dt.max, pars.fft.just.t$lo, 0)

## same as above (just comparing)
paste(ans.fftR.just.t[[2]][,1], collapse=", ") # (again, same expectedSp1EsAfterPropT =) 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743589792, 0.000299735813540893, 0.000299734881753716, 0.000299733948515344, 0.00029973301411462, 0.00029973207884177, 0.000299731142988294, 0.000299730206846208, 0.00029972927070804, 0.000299728334866242, 0.000299727399612858, 0.000299726465239364, 0.000299725532035927, 0.000299724600291355, 0.000299723670292635, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ans.fftR.just.t[[2]][,2], collapse=", ") # (again, same expectedSp1DsAfterPropT =) 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.98929572353432, 3.87688349797518, 4.83086427268124, 5.78355856417974, 6.65263439389469, 7.35225216820117, 7.80684129110362, 7.96450018578724, 7.80674374052118, 7.35206845754541, 6.65238511637526, 5.78326972079126, 4.83056283595408, 3.87659337350287, 2.98903491521347, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0



# (9) testIntegrateOneBranchLoRes32BinsOutsideClassJustX

args.fft.just.x <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x <- quasse.extent.debug(control.fft.32, drift, diffusion) # prepares X axis stuff
ext.fft.just.x$padding # 16/16, 4/4

pars.fft.just.x <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.just.x, ext.fft.just.x, pars.just.x) # adds lambda and mu vectors to ext.fft.just.x

# initialization
vars.fft.just.x <- matrix(0, control.fft.32$nx, 2) # low resolution
vars.fft.just.x[seq_len(ext.fft.just.x$ndat[2]),2] <- dnorm(ext.fft.just.x$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x$lo$padding should be 0.0
paste(vars.fft.just.x[,2], collapse=", ") # expectedInitialDs = 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0

# checking kernel matches (see PropagatesQuaSSETest -> testMakeNormalKernInPlaceAndFFtAndIfft)
kern.just.x <- fftR.make.kern.debug(-control.fft.32$dt * pars.fft.just.x$hi$drift, sqrt(control.fft.32$dt * pars.fft.just.x$hi$diffusion), control.fft.32$nx, control.fft.32$dx, pars.fft.just.x$hi$padding[1], pars.fft.just.x$hi$padding[2])
paste(kern.just.x, collapse=", ") # 0.986703287028858, 0.00664835445182386, 2.03374705433156e-09, 2.82445649260927e-20, 1.78085279698565e-35, 5.09772422059472e-55, 6.62490770689586e-79, 3.90875553004076e-107, 1.04701370374391e-139, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.04701370374391e-139, 3.90875553004076e-107, 6.62490770689586e-79, 5.09772422059472e-55, 1.78085279698565e-35, 2.82445649260927e-20, 2.03374705433156e-09, 0.00664835445182386

# checking fft-ed fY matches
fy.just.x <- fft(kern.just.x)
paste(Re(fy.just.x), collapse=", ") # expectedFFTedfY = 1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216, 0.981614853950298, 0.979316029808302, 0.977301093995625, 0.975647479188405, 0.97441873269917, 0.973662074416229, 0.973406582192705, 0.973662074416229, 0.97441873269917, 0.975647479188405, 0.977301093995625, 0.979316029808302, 0.981614853950298, 0.984109224049216, 0.986703282961364, 0.989297342492751, 0.991791714355113, 0.994090541136289, 0.996105480062091, 0.997759097982437, 0.998987847110852, 0.999744507157237

# checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x, 2, fft)
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x
Re(cstep.2[,2])[1:10] # matches with prints
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # matches with prints

ds.prop.x <- fftR.propagate.x.debug(vars.fft.just.x, control.fft.32$nx, fy.just.x, pars.fft.just.x$lo$padding[1], pars.fft.just.x$lo$padding[2])
paste(ds.prop.x[,1], collapse=", ") # expectedSp1EsAfterPropX = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ds.prop.x[,2], collapse=", ") # expectedSp1DsAfterPropX = 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.99530083804265, 3.88416335936269, 4.83940600155165, 5.79327421787607, 6.66336349661662, 7.36377090119626, 7.81887626199731, 7.97674483551017, 7.81887626199731, 7.36377090119625, 6.66336349661661, 5.79327421787607, 4.83940600155166, 3.88416335936269, 2.99530083804265, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0



## (10) testIntegrateOneBranchLoRes1024BinsOutsideClassJustX (and 10.1 testPropagateChOneCh1024QuaSSETest)

control.fft.just.x.1024 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.01, # dt
                    nx=1024, # number of X bins
                    dx=0.0005, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

sd <- 0.05
drift <- 0.0
diffusion <- 0.001
lambda <- sigmoid.x
mu <- constant.x

args.fft.just.x.1024 <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.just.x.1024 <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.just.x.1024 <- quasse.extent.debug(control.fft.just.x.1024, drift, diffusion) # prepares X axis stuff
# ext.fft.just.x.1024$padding

pars.fft.just.x.1024 <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.just.x.1024, ext.fft.just.x.1024, pars.just.x.1024) # adds lambda and mu vectors to ext.fft.just.x

## initialization
vars.fft.just.x.1024 <- matrix(0, control.fft.just.x.1024$nx, 2) # low resolution
vars.fft.just.x.1024[seq_len(ext.fft.just.x.1024$ndat[2]),2] <- dnorm(ext.fft.just.x.1024$x[[2]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.1024$lo$padding should be 0.0
print(paste(vars.fft.just.x.1024[1:48,2], collapse=",")) # expectedInitialDs = 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771

kern.just.x.1024 <- fftR.make.kern.debug(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$lo$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$lo$diffusion), control.fft.just.x.1024$nx, control.fft.just.x.1024$dx, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])

fy.just.x.1024 <- fft(kern.just.x.1024)
print(paste(Re(fy.just.x.1024)[1:48], collapse=", ")) # expectedFFTedfY = 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368548, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577788, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032

## checking part of convolution matches
cstep.1 <- apply(vars.fft.just.x.1024, 2, fft)
Re(cstep.1[,2])[1:10]
cstep.2 <- cstep.1 * fy.just.x.1024
Re(cstep.2[,2])[1:10] # matches with prints
cstep.3 <- apply(cstep.2, 2, ifft)
Re(cstep.3[,2])[1:10] # matches with prints (depending on CPU architecture, this differs if numbers get too small!)

ds.prop.x.1024 <- fftR.propagate.x.debug(vars.fft.just.x.1024, control.fft.just.x.1024$nx, fy.just.x.1024, pars.fft.just.x.1024$lo$padding[1], pars.fft.just.x.1024$lo$padding[2])
print(paste(ds.prop.x.1024[1:48,2], collapse=", ")) # expectedSp1DsFirst48AfterPropX = 0.000365714984190948, 0.000382414193856355, 0.000399835934138456, 0.00041800955800901, 0.000436965523466329, 0.000456735431302938, 0.000477352064003592, 0.000498849425801072, 0.000521262783917567, 0.000544628711019852, 0.000568985128916886, 0.000594371353528846, 0.000620828141157005, 0.000648397736084276, 0.000677123919536558, 0.000707052060035464, 0.000738229165173324, 0.000770703934841743, 0.000804526815945299, 0.000839750058632349, 0.000876427774075162, 0.000914615993832026, 0.000954372730824099, 0.000995758041960245, 0.0010388340924432, 0.0010836652217908, 0.00113031801160615, 0.0011788613551308, 0.00122936652861539, 0.00128190726454212, 0.00133655982673381, 0.00139340308738429, 0.00145251860604505, 0.00151399071060322, 0.00157790658028586, 0.00164435633072572, 0.00171343310112364, 0.00178523314354266, 0.00185985591436892, 0.00193740416797439, 0.00201798405261629, 0.00210170520860801, 0.00218868086879601, 0.00227902796137729, 0.00237286721509132, 0.00247032326682047, 0.00257152477163242, 0.00267660451529771

print(paste(ds.prop.x.1024[848:895,2], collapse=", ")) # expectedSp1DsLast48AfterPropX = 0.00267660451529771, 0.00257152477163242, 0.00247032326682047, 0.00237286721509132, 0.0022790279613773, 0.00218868086879601, 0.00210170520860801, 0.0020179840526163, 0.00193740416797439, 0.00185985591436892, 0.00178523314354266, 0.00171343310112364, 0.00164435633072573, 0.00157790658028586, 0.00151399071060322, 0.00145251860604505, 0.00139340308738429, 0.00133655982673381, 0.00128190726454212, 0.00122936652861539, 0.0011788613551308, 0.00113031801160615, 0.0010836652217908, 0.0010388340924432, 0.000995758041960245, 0.000954372730824099, 0.000914615993832026, 0.000876427774075162, 0.000839750058632349, 0.000804526815945299, 0.000770703934841743, 0.000738229165173324, 0.000707052060035464, 0.000677123919536558, 0.000648397736084276, 0.000620828141157005, 0.000594371353528846, 0.000568985128916886, 0.000544628711019852, 0.000521262783917567, 0.000498849425801072, 0.000477352064003592, 0.000456735431302938, 0.000436965523466329, 0.00041800955800901, 0.000399835934138456, 0.000382414193856355, 0.000365714984190948



## (11) testIntegrateOneBranchHiRes4096BinsOutsideClassJustX (note that I'm using 4096 as low res here, but I'll co-opt 4096 high res in Java)  (and 11.1 testPropagateChOneCh1024QuaSSETest)

## requires control.fft.just.x.1024 and pars.fft.just.x.1024 from (10)

## initialization
vars.fft.just.x.4096 <- matrix(0, control.fft.just.x.1024$nx * control.fft.just.x.1024$r, 2) # low resolution
vars.fft.just.x.4096[seq_len(ext.fft.just.x.1024$ndat[1]),2] <- dnorm(ext.fft.just.x.1024$x[[1]], 0.0, sd) # states are the qu character values observed at the tips (assuming observed state is 0.0), note that the last pars.fft.just.x.4096$lo$padding should be 0.0
paste(vars.fft.just.x.4096[1:48,2], collapse=", ") # expectedInitialDsFirst48 = 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846

kern.just.x.4096 <- fftR.make.kern.debug(-control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$drift, sqrt(control.fft.just.x.1024$dt * pars.fft.just.x.1024$hi$diffusion), control.fft.just.x.1024$nx * control.fft.just.x.1024$r, control.fft.just.x.1024$dx / control.fft.just.x.1024$r, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2]) # note the division of dx / r to do high resolution!

fy.just.x.4096 <- fft(kern.just.x.4096)
print(paste(Re(fy.just.x.4096)[1:48], collapse=", ")) # expectedFFTedfY = 1, 0.999247292368191, 0.996992567179915, 0.993245992007481, 0.988024427909734, 0.981351303026827, 0.973256437474405, 0.963775821368549, 0.952951348293497, 0.940830506973933, 0.92746603432656, 0.912915533436717, 0.897241060330147, 0.880508683684162, 0.862788021843104, 0.844151761668242, 0.824675163860569, 0.804435559446082, 0.783511842107391, 0.761983960984176, 0.739932418450195, 0.717437777209001, 0.694580180837756, 0.671438891652661, 0.64809184947515, 0.624615254550175, 0.601083177512084, 0.577567198915383, 0.554136080452835, 0.530855469577789, 0.507787638837061, 0.484991260810779, 0.462521219151724, 0.440428455824044, 0.418759854264325, 0.39755815783139, 0.376861922578424, 0.356705503075537, 0.3371190697353, 0.318128655850257, 0.299756232341542, 0.282019808042489, 0.264933553200873, 0.248507943778199, 0.232749924053504, 0.217663085001486, 0.203247855908884, 0.189501706717032

ds.prop.x.4096 <- fftR.propagate.x.debug(vars.fft.just.x.4096, control.fft.just.x.1024$nx * control.fft.just.x.1024$r, fy.just.x.4096, pars.fft.just.x.1024$hi$padding[1], pars.fft.just.x.1024$hi$padding[2])
print(paste(ds.prop.x.4096[1:48,2], collapse=", ")) # expectedSp1DsFirst48AfterPropX = 0.000353647683540531, 0.000357627448646487, 0.000361649739618714, 0.000365714984190948, 0.000369823614096463, 0.000373976065102163, 0.000378172777042955, 0.000382414193856355, 0.000386700763617376, 0.000391032938573662, 0.000395411175180895, 0.000399835934138456, 0.000404307680425366, 0.000408826883336479, 0.000413394016518956, 0.00041800955800901, 0.000422673990268911, 0.00042738780022428, 0.000432151479301657, 0.000436965523466329, 0.000441830433260469, 0.000446746713841521, 0.000451714875020898, 0.000456735431302938, 0.000461808901924177, 0.00046693581089287, 0.000472116687028841, 0.000477352064003592, 0.000482642480380733, 0.000487988479656683, 0.000493390610301674, 0.000498849425801072, 0.000504365484696964, 0.000509939350630071, 0.000515571592381964, 0.000521262783917567, 0.00052701350442799, 0.000532824338373647, 0.000538695875527712, 0.000544628711019852, 0.000550623445380317, 0.000556680684584297, 0.000562801040096648, 0.000568985128916886, 0.000575233573624552, 0.000581547002424852, 0.000587926049194663, 0.000594371353528846
print(paste(ds.prop.x.4096[1001:1048,2], collapse=", ")) # expectedSp1DsLater48AfterPropX = 1.12964601442883, 1.13523957985681, 1.14085371385815, 1.14648844781881, 1.1521438128928, 1.15781984000018, 1.16351655982519, 1.16923400281429, 1.17497219917421, 1.18073117887008, 1.18651097162341, 1.19231160691021, 1.19813311395904, 1.20397552174904, 1.20983885900802, 1.21572315421049, 1.22162843557575, 1.22755473106588, 1.23350206838385, 1.23947047497152, 1.24545997800775, 1.25147060440637, 1.25750238081428, 1.26355533360948, 1.26962948889911, 1.27572487251749, 1.28184151002416, 1.28797942670193, 1.29413864755491, 1.30031919730656, 1.30652110039772, 1.31274438098464, 1.31898906293703, 1.32525516983611, 1.33154272497262, 1.33785175134486, 1.34418227165675, 1.35053430831583, 1.35690788343134, 1.36330301881221, 1.36971973596514, 1.37615805609261, 1.38261800009092, 1.38909958854822, 1.39560284174258, 1.40212777963999, 1.40867442189243, 1.41524278783588



## (12) testBothXandTPropagateMethodsInsideClassOneBranchLoRes32Bins

## requires control.fft.32 from (8)

lambda <- sigmoid.x
mu <- constant.x
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

args.fft.32.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.32.both <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.32.both <- quasse.extent.debug(control.fft.32, drift, diffusion) # prepares X axis stuff

vars.fft.32.both <- matrix(0, control.fft.32$nx, 2) # low resolution
vars.fft.32.both[seq_len(ext.fft.32.both$ndat[2]),2] <- dnorm(ext.fft.32.both$x[[2]], 0.0, sd)
paste(c(vars.fft.32.both[seq_len(ext.fft.32.both$ndat[2]),2], rep(0, (32-length(ext.fft.32.both$x[[2]])))), collapse=", ") # expectedInitialDs = 0.709491856924629, 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522966, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522966, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0.709491856924629, 0, 0, 0, 0, 0, 0, 0, 0, 0

pars.fft.32.both <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.32.both, ext.fft.32.both, pars.32.both)

pde.fftR.32.both <- with(control.fft.32, make.pde.quasse.fftR.debug2(nx, dx, dt.max, 2L)) # creates function

ans.fftR.32.both <- pde.fftR.32.both(vars.fft.32.both, control.fft.32$dt.max, pars.fft.32.both$lo, 0)

## expectedFFTedfY = fy = 1, 0.999744507157237, 0.998987847110852, 0.997759097982437, 0.996105480062091, 0.994090541136289, 0.991791714355113, 0.989297342492751, 0.986703282961364, 0.984109224049216, 0.981614853950298, 0.979316029808302, 0.977301093995625, 0.975647479188405, 0.97441873269917, 0.973662074416229, 0.973406582192705, 0.973662074416229, 0.97441873269917, 0.975647479188405, 0.977301093995625, 0.979316029808302, 0.981614853950298, 0.984109224049216, 0.986703282961364, 0.989297342492751, 0.991791714355113, 0.994090541136289, 0.996105480062091, 0.997759097982437, 0.998987847110852, 0.999744507157237

## E's (note how first 4 and last 4 do not change after propagate in t because they are flanking bins on the left and right)
paste(ans.fftR.32.both[[2]][,1], collapse=", ") # expectedSp1EsAfterPropTandX = 0.000299740440744498, 0.000299739520470888, 0.000299738597335782, 0.00029973767161558, 0.000299736743576341, 0.000299735813529336, 0.000299734881744068, 0.000299733948507616, 0.000299733014108822, 0.000299732078837909, 0.000299731142986375, 0.000299730206846234, 0.000299729270710011, 0.000299728334870154, 0.000299727399618708, 0.000299726465247143, 0.000299725532045626, 0.000299724600302962, 0.000299723670306136, 0.000299722742324704, 0.000299721816669763, 0.000299720893607378, 0.00029971997341377, 0, 0, 0, 0, 0, 0, 0, 0, 0

## D's (note how first 4 and last 4 do not change after propagate in t because they are flanking bins on the left and right)
print(paste(exp(ans.fftR.32.both[[1]]) * ans.fftR.32.both[[2]][,2], collapse=", ")) # expectedSp1DsAfterPropTandX = 0.708264611249434, 1.07794488915543, 1.57625248872262, 2.21453845459856, 2.99004586159941, 3.87732490267096, 4.83085571964462, 5.78300263831972, 6.65150777531997, 7.35062312893061, 7.80486719135139, 7.96240319031702, 7.80476971649718, 7.35043955518597, 6.65125867066456, 5.78271397425439, 4.83055444185051, 3.87703489789509, 2.98978512541793, 2.21431781312691, 1.57607596745573, 1.07781089197137, 0.708167869605178, 0, 0, 0, 0, 0, 0, 0, 0, 0



## (13) testBothXandTPropagateMethodsInsideClassOneBranchLoRes32BinsDt002

control.fft.32.dt002 <- list(tc=100.0, # time point at which we go from high -> low resolution of X
                    dt.max=0.02, # dt
                    nx=32, # number of X bins
                    dx=0.01, # size of each X bin
                    r=4L, # high res of X = nx * r
                    xmid=0, # sp rate ~ X is a logistic regression, and xmid is the value of X at the inflection point
                    w=10, # used to determine nkl and nkr (see quasse.extent)
                    flags=0L, # FFT stuff below
                    verbose=0L,
                    atol=1e-6,
                    rtol=1e-6,
                    eps=1e-3,
                    method="fftC")

lambda <- sigmoid.x
mu <- constant.x
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

args.fft.32.dt002.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.32.dt002.both <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.32.dt002.both <- quasse.extent.debug(control.fft.32.dt002, drift, diffusion) # prepares X axis stuff

vars.fft.32.dt002.both <- matrix(0, control.fft.32.dt002$nx, 2) # low resolution
vars.fft.32.dt002.both[seq_len(ext.fft.32.dt002.both$ndat[2]),2] <- dnorm(ext.fft.32.dt002.both$x[[2]], 0.0, sd)
print(paste(c(vars.fft.32.dt002.both[seq_len(ext.fft.32.dt002.both$ndat[2]),2], rep(0, (32-length(ext.fft.32.dt002.both$x[[2]])))), collapse=", ")) # expectedInitialDs = 1.07981933026376, 1.57900316601788, 2.21841669358911, 2.9945493127149, 3.88372109966426, 4.83941449038287, 5.79383105522965, 6.66449205783599, 7.36540280606647, 7.82085387950912, 7.97884560802865, 7.82085387950912, 7.36540280606647, 6.66449205783599, 5.79383105522965, 4.83941449038287, 3.88372109966426, 2.9945493127149, 2.21841669358911, 1.57900316601788, 1.07981933026376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

pars.fft.32.dt002.both <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.32.dt002.both, ext.fft.32.dt002.both, pars.32.dt002.both)

pde.fftR.32.dt002.both <- with(control.fft.32.dt002, make.pde.quasse.fftR.debug2(nx, dx, dt.max, 2L)) # creates function

ans.fftR.32.dt002.both <- pde.fftR.32.dt002.both(vars.fft.32.dt002.both, control.fft.32.dt002$dt.max, pars.fft.32.dt002.both$lo, 0)

## expectedFFTedfY = fy = 1, 0.997284635663215, 0.989243568247244, 0.976187735593861, 0.958621745686731, 0.937224046438699, 0.91282033618305, 0.886351315420652, 0.858836097834951, 0.831332753409919, 0.80489754454845, 0.780544437378578, 0.759206428540079, 0.741700129042481, 0.728694899474875, 0.720687643959832, 0.717984152783717, 0.720687643959832, 0.728694899474875, 0.741700129042481, 0.759206428540079, 0.780544437378578, 0.80489754454845, 0.831332753409919, 0.858836097834951, 0.886351315420652, 0.91282033618305, 0.937224046438699, 0.958621745686731, 0.976187735593861, 0.989243568247244, 0.997284635663215

## E's (note how first 5 and last 5 do not change after propagate in t because they are flanking bins on the left and right)
## expectedSp1EsAfterPropT = E's after propagate in t = 0.000598958856744621, 0.000598955169219663, 0.000598951471383546, 0.000598947764353066, 0.000598944049256322, 0.000598940327231526, 0.000598936599425321, 0.000598932866991505, 0.000598929131089762, 0.000598925392884107, 0.000598921653541243, 0.000598917914229581, 0.000598914176117241, 0.000598910440370952, 0.000598906708154486, 0.000598902980627182, 0.000598899258942585, 0.000598895544247008, 0.000598891837677968, 0.000598888140363113, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(ans.fftR.32.dt002.both[[2]][,1], collapse=", ") # expectedSp1EsAfterPropTandX = E's after propagate in t and x = 0.000598958856744621, 0.000598955169219663, 0.000598951471383546, 0.000598947764353066, 0.000598944049256322, 0.000598940326823012, 0.000598936599098334, 0.000598932866746461, 0.000598929130926967, 0.000598925392803752, 0.000598921653543447, 0.000598917914314325, 0.000598914176284426, 0.000598910440620369, 0.000598906708485821, 0.000598902981040026, 0.000598899258942585, 0.000598895544247008, 0.000598891837677968, 0.000598888140363113, 0.000598884453418597, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0

## D's (note how first 5 and last 5 do not change after propagate in t because they are flanking bins on the left and right)
## expectedSp1DsAfterPropT = D's after propagate in t = 1.07607462857111, 1.57350796410007, 2.2106689156898, 2.98405395410774, 3.87006132459304, 4.82233340368848, 5.77330938669843, 6.64080371901192, 7.33913154877464, 7.79286078082761, 7.95018769794391, 7.79266608851656, 7.33876489743091, 6.64030620874454, 5.77273291052126, 4.82173179361503, 3.86948229161728, 2.98353343056498, 2.21022855749815, 1.5731556614014, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
paste(exp(ans.fftR.32.dt002.both[[1]]) * ans.fftR.32.dt002.both[[2]][,2], collapse=", ") # expectedSp1DsAfterPropTandX = D's after propagate t and x = 1.07607462857111, 1.57350796410007, 2.2106689156898, 2.98405395410774, 3.87006132459304, 4.82224125131593, 5.76741044168543, 6.62885082364537, 7.32184914756387, 7.77191831137883, 7.92794195899393, 7.77172522522963, 7.3214854000561, 6.62835697979005, 5.76683776896087, 4.8216430122779, 3.86948229161728, 2.98353343056498, 2.21022855749815, 1.5731556614014, 1.07580719581128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0


## (14) testIntegrateOneBranchHiRes32BinsInsideClassBothXandTDt002

## requires control.fft.48.dt002 from (13)

lambda <- sigmoid.x
mu <- constant.x
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

args.fft.32.dt002.both <- list(lambda=1:4, mu=5, drift=6, diffusion=7) # index of each parameter in args
pars.32.dt002.both <- c(.1, .2, 0, 2.5, .03, drift, diffusion) # specifies parameter values
# note that y0=0.1, y1=0.2, xmid=0.0, and r=2.5 for the sigmoid function for lambda

ext.fft.32.dt002.both <- quasse.extent.debug(control.fft.32.dt002, drift, diffusion) # prepares X axis stuff

vars.fft.32.dt002.both <- matrix(0, control.fft.32.dt002$nx * 4, 2) # high resolution
vars.fft.32.dt002.both[seq_len(ext.fft.32.dt002.both$ndat[1]),2] <- dnorm(ext.fft.32.dt002.both$x[[1]], 0.0, sd)

print(paste(vars.fft.32.dt002.both[1:10,2], collapse=", ")) # expectedInitialDsFirst10 = 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774
print(paste(vars.fft.32.dt002.both[78:87,2], collapse=", ")) # expectedInitialDsLater10 = 1.88098154753774, 1.72554637653023, 1.57900316601788, 1.44129748672436, 1.31231629549353, 1.19189412137632, 1.07981933026376, 0.975840371583657, 0.879671919608545, 0.791000831787405

## have to normalize it to match behavior inside likelihood class
vars.fft.32.dt002.both[,2] <- vars.fft.32.dt002.both[,2] / (sum(vars.fft.32.dt002.both[,2]) * 0.01 / 4)
print(paste(vars.fft.32.dt002.both[1:10,2], collapse=", ")) # 0.815139677630763, 0.906516726853789, 1.00561993609328, 1.11277200402141, 1.22826696360519, 1.35236404194541, 1.48528133155527, 1.62718935303434, 1.77820459925865, 1.93838316051407
print(paste(vars.fft.32.dt002.both[78:87,2], collapse=", ")) # 1.93838316051407, 1.77820459925865, 1.62718935303434, 1.48528133155527, 1.35236404194541, 1.22826696360519, 1.11277200402141, 1.00561993609329, 0.90651672685379, 0.815139677630764

pars.fft.32.dt002.both <- diversitree:::expand.pars.quasse(lambda, mu, args.fft.32.dt002.both, ext.fft.32.dt002.both, pars.32.dt002.both)

pde.fftR.32.dt002.both <- with(control.fft.32.dt002, make.pde.quasse.fftR.debug2(32 * 4, dx/4, dt.max, 2L)) # note the 32 * 4 to do high res, as well as the dx/4

ans.fftR.32.dt002.both <- pde.fftR.32.dt002.both(vars.fft.32.dt002.both, control.fft.32.dt002$dt.max, pars.fft.32.dt002.both$hi, 0)
ans.fftR.32.dt002.both[[2]][,2] <- ans.fftR.32.dt002.both[[2]][,2] * exp(ans.fftR.32.dt002.both[[1]]) # unnormalizing it from what make.pde.quasse.fftR.debug2 does
ans.fftR.32.dt002.both[[2]][,2] <- ans.fftR.32.dt002.both[[2]][,2] / (sum(ans.fftR.32.dt002.both[[2]][,2]) * 0.0025) # re-normalize it to match the java class

## E's
print(paste(ans.fftR.32.dt002.both[[2]][1:10,1], collapse=", ")) # expectedSp1EsAfterPropTandXFirst10 = first 10 after t and x = 0.000598961614956886, 0.000598960696294835, 0.000598959776885031, 0.000598958856744621, 0.00059895793589067, 0.000598957014340337, 0.000598956092110961, 0.000598955169219663, 0.000598954245683891, 0.000598953321520975
print(paste(ans.fftR.32.dt002.both[[2]][78:87,1], collapse=", ")) # expectedSp1EsAfterPropTandXLater10 = first 10 after t and x = 0.000598889987793744, 0.000598889063763026, 0.000598888140363113, 0.000598887217611451, 0.000598886295525308, 0.000598885374121974, 0.000598884453418597, 0.0005988835334325, 0.000598882614180599, 0.00059888169568006

## D's
print(paste(ans.fftR.32.dt002.both[[2]][1:10,2], collapse=", ")) # expectedSp1DsAfterPropTandXFirst10 = 0.816826443231728, 0.908389791093209, 1.00769467468308, 1.11506438544145, 1.23079348607031, 1.35514165820029, 1.48832736122564, 1.63052138243079, 1.78184036871647, 1.9423404395558
print(paste(ans.fftR.32.dt002.both[[2]][78:87,2], collapse=", ")) # expectedSp1DsAfterPropTandXLater10 = 1.94192952805763, 1.78145241269027, 1.63015631463554, 1.4879849556747, 1.35482154601338, 1.23049517164847, 1.11478726270473, 1.00743804313279, 0.908152871485605, 0.816608392676392



## (15, 16, 17 and 18)

lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

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
lik.C.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1) # for some reason, this matches the R result
lik.C.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
lik.R.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)
(ll.C.1 <- lik.C.1(pars)) ## -6.389708

## (for 17 and 18)
(ll.R.1 <- lik.R.1(pars)) ## -6.389642, also with make.quasse.debug for some reason # DIFFERS WITHIN DIVERSITREE

## run with make.quasse.debug

## (for 15) expectedObservedPriorProbAtRoot = root p = 0.000101936764975178, 0.000383226449472287, 0.0013299487715952, 0.00432788980772701, 0.0127710256471166, 0.0347934756336926, 0.0875169004358447, 0.203239244270493, 0.435756829076961, 0.862583772672889, 1.57644305300714, 2.65995745025329, 4.14371594632748, 5.95968961768647, 7.9136199323705, 9.70161705011622, 10.9806988668854, 11.474460683722, 11.0700752296872, 9.86017006942937, 8.10836608266039, 6.15597655521472, 4.30633004611629, 2.7843061138813, 1.66176905311156

## (for 18) esDsHiAtNodeInitial0[1] = 0.308986942687903, 0.350566009871371, 0.396747087835907, 0.447890605896858, 0.504364398303888, 0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911, 2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839, 3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101, 4.59764281368466, 4.83941449038287, 5.08118112938378, 5.3217049979751, 5.55969772261993, 5.79383105522966, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599, 6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876, 7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024, 7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778, 7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647, 7.21053924923296, 7.04130653528599, 6.85887710038768, 6.66449205783599, 6.45944719335828, 6.24507866733522, 6.02274864309609, 5.79383105522965, 5.55969772261993, 5.32170499797509, 5.08118112938378, 4.83941449038287, 4.59764281368466, 4.35704354065101, 4.11872537439949, 3.88372109966426, 3.65298170778044, 3.42737184095615, 3.20766654683839, 2.9945493127149, 2.7886113289072, 2.59035191331783, 2.40018001393971, 2.21841669358911, 2.04529849127956, 1.88098154753774, 1.72554637653023, 1.57900316601788, 1.44129748672436, 1.31231629549353, 1.19189412137632, 1.07981933026376, 0.975840371583655, 0.879671919608544, 0.791000831787404, 0.709491856924628, 0.634793036713349, 0.566540754832024, 0.504364398303888, 0.447890605896858, 0.396747087835907, 0.350566009871371, 0.308986942687903

## (for 18) esDsHiAtNodeInitial1[1] = 0.000254946647636669, 0.000319674822138109, 0.000399835934138456, 0.000498849425801072, 0.000620828141157003, 0.000770703934841743, 0.000954372730824099, 0.0011788613551308, 0.00145251860604505, 0.00178523314354266, 0.00218868086879601, 0.00267660451529771, 0.00326512817532484, 0.00397310942785545, 0.00482253160451986, 0.0058389385158292, 0.00705191364734891, 0.00849560541101504, 0.0102092994868837, 0.0122380386022755, 0.0146332892566062, 0.0174536539009152, 0.0207656259132282, 0.0246443833694604, 0.0291746160933349, 0.0344513787810736, 0.0405809611459954, 0.0476817640292969, 0.0558851682975889, 0.0653363811239983, 0.0761952419644361, 0.08863696823876, 0.102852818461079, 0.119050648395517, 0.137455333812279, 0.158309031659599, 0.181871250031821, 0.208418696288452, 0.238244872152104, 0.271659384673712, 0.308986942687903, 0.350566009871371, 0.396747087835906, 0.447890605896858, 0.504364398303888, 0.566540754832024, 0.634793036713348, 0.709491856924629, 0.791000831787404, 0.879671919608544, 0.975840371583655, 1.07981933026376, 1.19189412137632, 1.31231629549353, 1.44129748672436, 1.57900316601788, 1.72554637653023, 1.88098154753774, 2.04529849127956, 2.21841669358911, 2.40018001393971, 2.59035191331783, 2.7886113289072, 2.9945493127149, 3.20766654683839, 3.42737184095615, 3.65298170778044, 3.88372109966426, 4.11872537439949, 4.35704354065101, 4.59764281368466, 4.83941449038287, 5.08118112938378, 5.32170499797509, 5.55969772261993, 5.79383105522965, 6.02274864309609, 6.24507866733522, 6.45944719335828, 6.66449205783599, 6.85887710038768, 7.04130653528599, 7.21053924923296, 7.36540280606647, 7.50480693833876, 7.62775630921048, 7.73336233605698, 7.82085387950912, 7.88958661815778, 7.93905094954024, 7.96887828189528, 7.97884560802865, 7.96887828189528, 7.93905094954024, 7.88958661815778, 7.82085387950912, 7.73336233605698, 7.62775630921048, 7.50480693833876, 7.36540280606647, 7.21053924923296, 7.04130653528599, 6.85887710038768

## (for 18) expectedHiLoIdxs4Transfer = 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99

lik.C.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
(ll.C.1 <- lik.C.1(pars, root=ROOT.FLAT)) ## -6.389708

## (for 16) expectedFlatPriorProbAtRoot = root p = 3.2258064516129

## (19) testPruneTrifTree32Bins

set.seed(1)
tr <- tree.quasse(c(lambda, mu, char), max.taxa=3, x0=0, single.lineage=FALSE, verbose=TRUE)
# simplifying tree for unit test
tr$tip.state[1] <- 0.0 # ch state
tr$tip.state[2] <- 0.1 # ch state
tr$tip.state[3] <- 0.2 # ch state
tr$edge.length <- c(0.02, 0.01, 0.01, 0.01)  # branch lengths
tr$orig[,2] <- c(0.02, 0.01, 0.01, 0.01) # branch lengths
tr$orig[,4] <- c(0.0, NA, 0.1, 0.2) # ch states

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
## lik.C.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
lik.C.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
lik.R.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)
(ll.C.1 <- lik.C.1(pars)) ## -9.082284
(ll.C.1 <- lik.R.1(pars)) ## -9.085542 DIFFERENT WITHIN DIVERSITREE

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




## (21) testPruneBifTree004Height32Bins

lambda <- function(x) sigmoid.x(x, 0.1, 0.2,  0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)

set.seed(1)
tr <- tree.quasse(c(lambda, mu, char), max.taxa=15, x0=0, single.lineage=FALSE, verbose=TRUE)

nodes <- c("nd13", "nd9", "nd5")
split.t <- Inf

pars <- c(.1, .2, 0, 2.5, .03, 0, .01) # lambda y0, lambda y1, lambda xmid, lambda r, mu, drift, diffusion
pars4 <- unlist(rep(list(pars), 4))

sd <- 1/50

control.C.1 <- list(dt.max=1/200, method="fftC")
lik.C.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
(ll.C.1 <- lik.C.1(pars4)) # -61.27245

control.R.1 <- list(dt.max=1/200, method="fftR")
lik.R.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)
(ll.R.1 <- lik.R.1(pars)) # -61.27245

## I looked closely at this example:
## control.C.1 <- list(tc=1,
##                     dt.max=0.1,
##                     nx=32,
##                     dx=0.1,
##                     r=4L,
##                     xmid=0.0,
##                     w=5,
##                     flags=0L,
##                     verbose=0L,
##                     atol=1e-6,
##                     rtol=1e-6,
##                     eps=1e-3,
##                     method="fftC")
## Diversitree starts by pruning sp16 (nodeIdx = 7 in Java, for the XYZ unit test).
## Initial D's are very close, and so are the kernel and fy = FFTed kernel
## After propagating in t, the D's that enter propagate in X are close.
##
##
## Things start getting more substantially different after we do ((FFTed D's) * fy) --> into iFFT.
## When running the above, it's the line that prints "after FFT -> *fy -> ifft (before scaling).
## You'll notice that just the elements in the middle of the D vector match the Java implementation
## (which happen to be the larger D's). The initial and final D's, which are much smaller, will be
## different between the Java and the C implementation.
##
## The next steps keep or increase the differences wrt the Java implementation, as after scaling
## and inverse FFT, one implementation might yield a positive value, the other a negative value.
## These differences add up, also because negative values are zero-ed before the next round
## of propagate in T and X.
##
## I still need to run the Java unit test on a linux machine to see if the numbers
## remain different (in which case this is a language-dependent difference), or if the numbers
## match (in which case, this is again an architecture-dependent difference).

## (15, 16, 17 and 18)

lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)
drift <- 0.0
diffusion <- 0.001
sd <- 0.05

set.seed(1)
tr <- tree.quasse(c(lambda, mu, char), max.taxa=2, x0=0, single.lineage=FALSE, verbose=TRUE)
# simplifying tree for unit test
tr$tip.state[1] <- 0.0 # ch state
tr$tip.state[2] <- 0.1 # ch state
tr$edge.length <- c(0.025, 0.025)  # branch lengths
tr$orig[,2] <- c(0.04, 0.04) # branch lengths
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
## lik.C.1 <- make.quasse.debug(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1) # for some reason, this matches the R result
lik.C.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.C.1)
(ll.C.1 <- lik.C.1(pars)) ## -6.394299

lik.R.1 <- make.quasse(tr, tr$tip.state, sd, sigmoid.x, constant.x, control.R.1)
(ll.R.1 <- lik.R.1(pars)) ## -6.394235, also with make.quasse.debug for some reason # DIFFERS WITHIN DIVERSITREE
