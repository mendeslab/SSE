library(diversitree)

normalise <- function(x) x / sum(x)

ifft <- function(x) fft(x, inverse=TRUE)

sigmoid.x <- function (x, y0, y1, xmid, r) {
    to.add = (y1 - y0) / (1 + exp(r * (xmid - x)))
    y0 + to.add
}

quasse.extent.debug <- function (control, drift, diffusion) {
    nx <- control$nx
    dx <- control$dx
    dt <- control$dt.max
    xmid <- control$xmid
    r <- control$r
    w <- control$w
    if (control$method == "mol") {
        ndat <- nx * c(r, 1)
        padding <- NULL
    }
    else {
        mean <- drift * dt
        sd <- sqrt(diffusion * dt)
        nkl <- max(ceiling(-(mean - w * sd)/dx)) * c(r, 1)
        nkr <- max(ceiling((mean + w * sd)/dx)) * c(r, 1)

        print(paste("nx =", nx))
        print(paste("diffusion = ", diffusion, " drift = ", drift))
        print(paste("mean =", mean, "sd =", sd, "w =", w, "r =", r))
        print(paste("nkl =", nkl, "nkr =", nkr))

        ndat <- nx * c(r, 1) - (nkl + 1 + nkr)

        print(paste("ndat =", ndat))

        padding <- cbind(nkl, nkr)
        storage.mode(padding) <- "integer"
    }
    x0.2 <- xmid - dx * ceiling((ndat[2] - 1)/2)
    x0.1 <- x0.2 - dx * (1 - 1/r)
    print(paste("xMinLo =", x0.2, " xMinHi =", x0.1))
    x <- list(seq(x0.1, length.out = ndat[1], by = dx/r), seq(x0.2,
        length.out = ndat[2], by = dx))
    tr <- seq(r, length.out = ndat[2], by = r)
    list(x = x, padding = padding, ndat = ndat, tr = tr, nx = c(nx *
        r, nx))
}

fftR.make.kern.debug <- function(mean, sd, nx, dx, nkl, nkr) {
  kern <- rep(0, nx)
  xkern <- (-nkl:nkr)*dx
  ikern <- c((nx - nkl + 1):nx, 1:(nkr + 1))

  print("inside fftR.make.kern")
  print(paste0("dx = ", dx))
  print("delta.x = ")
  print(xkern)
  print(paste0("mean = ", mean, " sd = ", sd))
  print("dnorm = ")
  print(dnorm(xkern, mean, sd))
  print(paste0("sum = ", sum(dnorm(xkern, mean, sd))))

  kern[ikern] <- normalise(dnorm(xkern, mean, sd))
  print("kern = ")
  print(kern)
  kern
}

fftR.propagate.t.debug <- function(vars, lambda, mu, dt, ndat) {
  i <- seq_len(ndat) # creates a seq that starts with 1 and ends in ndat, with increments of 1
  r <- lambda - mu
  z <- exp(dt * r)
  e0 <- vars[i,1]
  d0 <- vars[i,-1]
  ## print(paste("mu =", mu, "lambda =", lambda, "e0 =", e0, "z =", z))
  vars[i,1] <- (mu + z*(e0 - 1)*mu - lambda*e0) / (mu + z*(e0 - 1)*lambda - lambda*e0) # E's
  ## print(paste("e =", vars[i,1]))
  dd <- (z * r * r)/(z * lambda - mu + (1-z)*lambda*e0)^2

  # for debugging
  ## print(paste("numerator", (z * r * r)))
  ## print(paste("denominator", (z * lambda - mu + (1-z)*lambda*e0)))
  ## print(paste("numerator/denominator", dd))

  ## print(paste("scratch", dd))
  ## print(paste("d0", d0))
  vars[i,-1] <- dd * d0 # D's are set in place (note that E's are not!!! can only capture E if we throw the return of this function into a variable)
  ## print(paste("dd * d0", vars[i,-1]))
  vars
}

fftR.propagate.x.debug <- function(vars, nx, fy, nkl, nkr) {

    ## print(paste("nkl=", nkl, " nkr=", nkl, " nx=", nx))
    ## print("vars="); print(vars)

    ## print("first FFT of quant trait = "); print(apply(vars, 2, fft))
    ## print("kernel = "); print(fy)
    ## print("first FFT of quant trait * kernel = "); print(apply(vars, 2, fft) * fy)
    ## print("after FFT -> *fy -> ifft (before scaling)"); print(apply(apply(vars, 2, fft) * fy, 2, ifft))
    ## print(paste0("scaling uses 1/nx, with nx being ", nx))
    vars.out <- Re(apply(apply(vars, 2, fft) * fy, 2, ifft))/nx

    print("vars after fft-ing, multiplying by fy kernel, i-ffting, and scaling"); print(vars.out)
    ndat <- nx - (nkl + 1 + nkr) # this plus 1 here is confusing, need to ask Xia
    i.prev.l <- 1:nkl
    i.prev.r <- (ndat-nkr+1):ndat
    i.zero <- (ndat+1):nx

    # print(c(i.prev.l, i.prev.r))

    vars.out[c(i.prev.l, i.prev.r),] <- vars[c(i.prev.l, i.prev.r),]
    vars.out[i.zero,] <- 0
    vars.out[vars.out < 0] <- 0

    ## print("vars after rearranging it"); print(vars.out)

    vars.out
}

root.p.quasse.debug <- function(d.root, pars, root, root.f) {
    if (!is.null(root.f) && root != ROOT.GIVEN)
        warning("Ignoring specified root state")
    x <- pars$x
    dx <- x[2] - x[1]
    if (root == ROOT.FLAT) {
        print("root d's =")
        print(paste(d.root, collapse=", "))
        p <- 1/((pars$nx - 1) * dx)
        print("root p =")
        print(paste(p, collapse=", "))
    }
    else if (root == ROOT.OBS) {
        print("getting prior at root")
        print(paste0("denominator = ", (sum(d.root) * dx)))
        print(paste0("dx = ", dx))
        p <- d.root/(sum(d.root) * dx)
        print(paste0("root d's (size=", length(d.root), ") ="))
        print(paste(d.root, collapse=", "))
        print("root p =")
        print(paste(p, collapse=", "))
    }
    else if (root == ROOT.GIVEN) {
        p <- root.f(x)
    }
    else {
        stop("Unsupported root mode")
    }
    p
}

make.rootfunc.quasse.debug <- function(cache) {
    root.idx <- cache$root
    nx <- cache$control$nx
    dx <- cache$control$dx

    function(res, pars, condition.surv, root, root.f, intermediates) {
        ## vals and lq come from all.branches.list
        vals <- matrix(res$vals, nx, 2)[seq_len(pars$lo$ndat),] # just making matrix
        lq <- res$lq # all normalizing factors for all nodes in tree
        d.root <- vals[, 2] # normalized D's at root
        root.p <- root.p.quasse.debug(d.root, pars$lo, root, root.f) # prior at root for each trait value bin (several options, user-defined); the default is the weighted sum, which is (d.root/sum(d.root) ==> returns a vector)

        if (condition.surv) {
            lambda <- pars$lo$lambda
            e.root <- vals[, 1]

            print("root e's = ")
            print(paste(e.root, collapse=", "))

            ## param = lambda, mu
            ## data D = tree (fixed)
            ## tip states are boundary conditions
            ## P(D|param) = L(param;D)
            ## at this point, d.root = P(D,tree_exists|params) = P(D|params,tree_exists) * P(tree_exists)
            ## but what we want is P(D|params,tree_exists)
            ## so we need to d.root / P(tree_exists)

            print(paste0("lambda (size=", length(lambda), ") = "))
            print(paste(lambda, collapse=", "))
            ## Root-prior values are densities: multiply the denominator by dx before dividing D.
            ## This conditions on survival without introducing the historical extra dx^2 factor.
            p.surv <- sum(root.p * lambda * (1 - e.root)^2) * dx
            d.root <- d.root / p.surv
            ## e.root is the probability of extinction of a lineage starting at the root
            ## (1 - e.root) is the probability that lineage survives until today
            ## we square it because there are two lineages coming from the root (i.e., the definition of a tree)
            ## multiply by lambda because we need a sp'n event at the root generating those two lineages
            ## multiply by root.p to account for the trait value states (those influence the extinction probability)
            ## * dx: multiplying by bin size converts densities (median of bin is a density) into probabilities

            print(paste0("Denominator sum for conditioning on survival = ", sum(root.p * lambda * (1 - e.root)^2)))
            print(paste0("dx = ", dx))
            print(paste0("D's after conditioning on survival (size=", length(d.root), ") = "))
            print(paste(d.root, collapse=", "))
        }

        print(paste0("sum(lq) = ", sum(lq)))

        loglik <- log(sum(root.p * d.root) * dx) + sum(lq) # goes back to actual likelihood ("unnormalizing it")

        ## doesn't seem important (?)
        if (intermediates) {
            attr(loglik, "intermediates") <- res
            attr(loglik, "vals") <- vals
        }

        ## returns final result
        loglik
    }
}

make.initial.conditions.quasse.debug <- function(control) {
    tc <- control$tc
    r <- control$r
    nx.lo <- control$nx
    nx.hi <- nx.lo * r
    eps <- 1e-08
    function(init, pars, t, idx) {
        if (length(init[[1]]) != length(init[[2]]))
            stop("Data have incompatible length")
        if (t < tc) {
            nx <- nx.hi
            lambda <- pars[[1]]$lambda
        }
        else {
            nx <- nx.lo
            lambda <- pars[[2]]$lambda
        }
        ndat <- length(lambda)
        i <- seq_len(nx)
        j <- seq.int(nx + 1, nx + ndat)

        ## return
        ## init[[1]][i] is the E from child 1 (= E from child 2)
        ## init[[1]][j] is the D from child 1
        ## init[[2]][j] is the D from child 2
        ## and 0's for padding
        print(paste0("nx = ", nx))
        print(paste0("ndat = ", ndat))
        print("merging:")
        print(paste0("E's at merge = ", paste(init[[1]][i], collapse=", ")))
        print(paste0("Left D's at merge = ", paste(init[[1]][j], collapse=", ")))
        print(paste0("Right D's at merge = ", paste(init[[2]][j], collapse=", ")))
        print(paste0("Lambdas at merge = ", paste(lambda, collapse=", ")))
        print(c(init[[1]][i], init[[1]][j] * init[[2]][j] * lambda,
              rep.int(0, nx - ndat)))

        c(init[[1]][i], init[[1]][j] * init[[2]][j] * lambda,
            rep.int(0, nx - ndat))
    }
}

combine.branches.quasse.debug <- function (f.hi, f.lo, control) {

    nx <- control$nx
    dx <- control$dx
    tc <- control$tc
    r <- control$r
    eps <- log(control$eps)
    dt.max <- control$dt.max

    careful <- function(f, y, len, pars, t0, dt.max) {

        ## if low res, f() == f.lo
        ## if high res, f() == f.hi
        ans <- f(y, len, pars, t0)

        if (ans[[1]] > eps) {
            ans
        }
        else {
            if (control$method == "fftC" || control$method ==
                "fftR")
                dt.max <- dt.max/2
            len2 <- len/2
            ans1 <- Recall(f, y, len2, pars, t0, dt.max)
            ans2 <- Recall(f, ans1[[2]], len2, pars, t0 + len2,
                dt.max)
            ans2[[1]][[1]] <- ans1[[1]][[1]] + ans2[[1]][[1]]
            ans2
        }
    }

    function(y, len, pars, t0, idx) {

        print(paste0("len = ", len))

        if (t0 < tc) {
            dx0 <- dx/r
            nx0 <- nx * r

            ## print(paste("dx0 = ", dx0))
            ## print(paste("nx0 = ", nx0))

        }
        else {
            dx0 <- dx
            nx0 <- nx
        }
        if (any(y < -1e-08))
            stop("Actual negative D value detected -- calculation failure")

        y[y < 0] <- 0
        y <- matrix(y, nx0, 2) # if in low-res, this re-sizes it

        print("y inside combine.branches before normalization = ")
        print(y)

        q0 <- sum(y[, 2]) * dx0 # calculating normalization factor from D's
        if (q0 <= 0)
            stop("No positive D values")

        y[, 2] <- y[, 2]/q0 # normalizing D's, before we integrate (NOT LOG-fied!!!!)

        print("y inside combine.branches after normalization = ")
        print(paste(y[,2], collapse=", "))

        lq0 <- log(q0) # will be added at the very end
        print(paste0("log-normalization factor starts at ", log(q0)))

        ## option 1: initial time is already past tc, low-res
        if (t0 >= tc) {

            print("y inside combine.branches for lo = ")
            print(y)

            ## ans[[1]] carries the log(sum(D's * dx)), the log-normalization factor
            ans <- careful(f.lo, y, len, pars$lo, t0, dt.max) # integrate

            print("ans at lo, prior to adding log-normalization factor")
            print(ans)
        }

        ## option 2: initial time + dt is still lower than tc, high-res
        else if (t0 + len < tc) {

            print("y inside combine.branches for hi = ")
            print(y)

            ans <- careful(f.hi, y, len, pars$hi, t0, dt.max) # integrate

            print("option 2: ans at hi, prior to adding log-normalization factor")
            print(ans)
        }

        ## option 3: first segment in high-resolution (up until tc), then remaining in low-res
        ## note that option 3 will happen before option 1 -- in option 1, y has already been
        ## resized to low-res
        else {

            print("option 3: y inside combine.branches for hi = ")
            print(y)

            len.hi <- tc - t0 # this first stretch in high-res

            print(paste0("len.hi = ", len.hi))

            ## ans.hi[[1]] carries the log(sum(D's * dx)), the log-normalization factor
            ans.hi <- careful(f.hi, y, len.hi, pars$hi, t0, dt.max) # integrate first stretch at high-res

            ## but then here we grab just every r-th element in ans.hi (as given by a vector of indices separated by r in pars$tr)
            ## y.lo will be the first elements inside the low-res y -- the rest will be 0's (see below)
            y.lo <- ans.hi[[2]][pars$tr, ]

            ## I believe ans.hi[[1]] here was normalized and log-fied inside the function that generated ans.hi
            lq0 <- lq0 + ans.hi[[1]] # we add the high-res normalization factor to the total normalization factor
            print(paste0("added ", ans.hi[[1]], " to log-normalization factor"))
            print(paste0("log-normalization factor is now ", lq0))

            ## why would nrow(y.lo) < nx (nx here is the low-res nx)?
            ## we're putting 0's after the y.lo
            print("pars$tr = ")
            print(paste(pars$tr, collapse=", "))
            print("expectedHiLoIdxs4Transfer = ")
            print(paste(pars$tr-1, collapse=", "))
            if (nrow(y.lo) < nx) y.lo <- rbind(y.lo, matrix(0, nx - length(pars$tr), 2))

            print(paste0("len - len.hi = ", len - len.hi))

            ## ans[[1]] carries the log(sum(D's * dx)), the log-normalization factor
            ans <- careful(f.lo, y.lo, len - len.hi, pars$lo, tc, dt.max) # integrate remaining stretch at low-res
            ## Unlike with the high-res part, there's no updating of lq0 here...?

            print("ans at hi, but interval < dt, prior to adding log-normalization factor")
            print(ans)
        }

        print("esDs at return of combine.branches.quasse, after adding log-normalization factor")
        print(ans[[2]]) # lq0 here can be either high or low-res depending on which of the if/else if/else block executed, i.e., whether the branch is < >= tc

        print(paste0("log-normalization factor from integration only = ", ans[[1]]))
        print(paste0("adding to log-normalization factor = ", lq0))
        print(paste0("final log-normalization factor at return = ", ans[[1]]+lq0))
        ## ans[[1]] is the log((sum over all D's) * dx), the input of which was normalized; we need to keep adding to it; so at the very end of the pruning (at the root) we can unnormalize the likelihood
        ## ans[[2]] is the dataframe with the E's and D's
        c(ans[[1]] + lq0, ans[[2]])
    }
}

quasse.integrate.fftR.debug <- function (vars, lambda, mu, drift, diffusion, nstep, dt, nx,
    ndat, dx, nkl, nkr) {
    kern <- diversitree:::fftR.make.kern(-dt * drift, sqrt(dt * diffusion),
        nx, dx, nkl, nkr)
    fy <- fft(kern)
    for (i in seq_len(nstep)) {
        vars <- diversitree:::fftR.propagate.t(vars, lambda, mu, dt, ndat)
        # print(paste("nx = ", nx, " length(fy)=", length(fy)))
        # vars <- fftR.propagate.x(vars, nx, fy, nkl, nkr) # ignoring propagate X
    }

    vars
}

make.pde.quasse.fftR.debug <- function (nx, dx, dt.max, nd) {
    function(y, len, pars, t0) {
        padding <- pars$padding
        ndat <- length(pars$lambda)
        nt <- as.integer(ceiling(len/dt.max))
        dt <- len/nt

        if (!(length(y) %in% (nd * nx)))
            stop("Wrong size y")
        if (length(pars$lambda) != length(pars$mu) || length(pars$lambda) >
            (nx - 3))
            stop("Incorrect length pars")
        if (pars$diffusion <= 0)
            stop("Invalid diffusion parameter")
        if (!is.matrix(y))
            y <- matrix(y, nx, nd)
        ans <- quasse.integrate.fftR.debug(y, pars$lambda, pars$mu,
            pars$drift, pars$diffusion, nt, dt, nx, ndat, dx,
            padding[1], padding[2])
        q <- sum(ans[, 2]) * dx

        # print(ans[,1]) # E's
        # print(ans[,2]) ## FKM: D's that are returned from fftR.propagate.t

        # ans[,2] <- ans[,2]/q ## FKM: D's are normalized before returning

        list(log(q), ans)
    }
}

quasse.integrate.fftR.debug2 <- function (vars, lambda, mu, drift, diffusion, nstep, dt, nx,
    ndat, dx, nkl, nkr) {

    print("vars"); print(vars)
    print("lambda"); print(lambda)
    print("mu"); print(mu)
    print("drift"); print(drift)
    print("diffusion"); print(diffusion)
    print("nstep"); print(nstep)
    print("dt"); print(dt)
    print("nx"); print(nx)
    print("dx"); print(dx)
    print("nkl"); print(nkl)
    print("nkr"); print(nkr)

    print("before fftR.make.kern.debug")
    print(paste0("sd = sqrt(dt * diffusion) = ", sqrt(dt * diffusion)))
    kern = fftR.make.kern.debug(-dt * drift, sqrt(dt * diffusion),
        nx, dx, nkl, nkr)

    print("kern = "); print(kern)

    fy = fft(kern)

    print("fy")
    print(paste0(Re(fy), collapse=", "))

    for (i in seq_len(nstep)) {
        print("E's before propagate in t")
        print(paste(vars[,1], collapse=", "))
        print("D's before propagate in t")
        print(paste(vars[,2], collapse=", "))

        vars = diversitree:::fftR.propagate.t(vars, lambda, mu, dt, ndat)

        print("E's after propagate in t")
        print(paste(vars[,1], collapse=", "))
        print("D's after propagate in t")
        print(paste(vars[,2], collapse=", "))

        vars = fftR.propagate.x.debug(vars, nx, fy, nkl, nkr)

        print("E's after propagate in t and x")
        print(paste(vars[,1], collapse=", "))
        print("D's after propagate in t and x")
        print(paste(vars[,2], collapse=", "))
    }

    vars
}

make.pde.quasse.fftR.debug2 <- function (nx, dx, dt.max, nd) {
    function(y, len, pars, t0) {
        padding <- pars$padding

        ndat <- length(pars$lambda)
        nt <- as.integer(ceiling(len/dt.max)) # number of time intervals, rounded up
        dt <- len/nt # size of the interval;

        ## branch length = 10
        ## dt.max = 0.45
        ## 0.1 length left after 2 * dt.max
        ## when he ceiling, he gets nt = 3
        ## then dt = 10 / 3 = 0.33333333
        ## then there's no left over, and dt is by definition < dt.max

        print("y"); print(y)
        print("length(y)"); print(length(y))
        print("nd"); print(nd)
        print("nd * nx"); print(nd * nx)
        print("ndat"); print(ndat)

        if (!(length(y) %in% (nd * nx)))
            stop("Wrong size y")
        if (length(pars$lambda) != length(pars$mu) || length(pars$lambda) >
            (nx - 3))
            stop("Incorrect length pars")
        if (pars$diffusion <= 0)
            stop("Invalid diffusion parameter")
        if (!is.matrix(y))
            y <- matrix(y, nx, nd)

        ## y comes in normalized from combine.branches.quasse
        ## line: y[, 2] <- y[, 2]/q0 (the y here is normalized by q0 in that function)
        ans = quasse.integrate.fftR.debug2(y, pars$lambda, pars$mu,
            pars$drift, pars$diffusion, nt, dt, nx, ndat, dx,
            padding[1], padding[2])

        q = sum(ans[, 2]) * dx

        print("Normalization factor in make.pde.quasse.fftR = ")
        print(q)
        print("Returning log-Normalization factor in make.pde.quasse.fftR = ")
        print(log(q))

        # print(ans[,1]) # E's
        # print(ans[,2]) ## FKM: D's that are returned from fftR.propagate.t

        ans[,2] <- ans[,2]/q ## FKM: D's are normalized again before returning

        list(log(q), ans)
    }
}

make.branches.quasse.fftR.debug <- function(control) {
    nx <- control$nx
    dx <- control$dx
    r <- control$r
    dt.max <- control$dt.max
    tc <- control$tc

    ## make.pde.quasse.fftR returns a function, which in turn runs quasse.integrate.fftR, returning c(log(q), ans)
    f.hi <- make.pde.quasse.fftR.debug2(nx * r, dx/r, dt.max, 2L)
    f.lo <- make.pde.quasse.fftR.debug2(nx, dx, dt.max, 2L)

    combine.branches.quasse.debug(f.hi, f.lo, control)
}

initial.tip.quasse.debug <- function(cache, control, x) {
    nx <- control$nx * control$r
    npad <- nx - length(x)
    e0 <- 1 - cache$sampling.f
    if (control$tips.combined) {
        tips <- cache$tips
        t <- cache$len[tips]
        i <- order(t)
        target <- tips[i]
        states <- cache$states[i]
        states.sd <- cache$states.sd[i]
        y <- mapply(function(mean, sd) c(dnorm(x, mean, sd),
            rep(0, npad)), states, states.sd, SIMPLIFY = FALSE)
        y <- matrix(c(rep(e0, nx), unlist(y)), nx, length(target) + 1)

        ## print("if: y inside initial.tip.quasse = ");  print(y)

        list(target = target, y = y, t = t[i])
    }
    else {

        print("cache$states = "); print(cache$states)
        print("cache$states.sd = "); print(cache$states.sd)

        y <- mapply(function(mean, sd) c(rep(e0, nx), dnorm(x,
            mean, sd), rep(0, npad)), cache$states, cache$states.sd,
            SIMPLIFY = FALSE)

        print("else: y inside initial.tip.quasse = ");
        print(y)
        ## print(paste(y, collapse=", "))
        ## print("esDsHiAtNodeInitial0[1] = ")
        ## print(paste(y$sp1[129:231], collapse=", "))
        ## print("esDsHiAtNodeInitial1[1] = ")
        ## print(paste(y$sp2[129:231], collapse=", "))

        diversitree:::dt.tips.ordered(y, cache$tips, cache$len[cache$tips])
    }
}

all.branches.list.debug <- function (pars, cache, initial.conditions, branches, preset) {
    len <- cache$len
    depth <- cache$depth
    children <- cache$children
    order <- cache$order[-length(cache$order)]
    root <- cache$root
    n <- length(len)
    lq <- rep(0, n)
    n.tip <- cache$n.tip
    y <- cache$y

    branch.init <- branch.base <- vector("list", n)

    ## not sure what this is for
    if (!is.null(preset)) {
        lq[preset$target] <- preset$lq
        branch.base[preset$target] <- preset$base
    }

    ## doing all tips first
    if (is.null(names(y))) {
        ## x is a tip
        for (x in y) {
            if (!is.null(x)) {
                idx <- x$target
                branch.init[idx] <- list(x$y)
                ans <- branches(x$y, x$t, pars, 0, idx) # return from combine.branches.quasse()
                lq[idx] <- unlist(lapply(ans, "[[", 1)) # getting the log((sum over all D's) * dt) for this node
                branch.base[idx] <- lapply(ans, "[", -1)
            }
        }
    }

    ## tip (not used)
    else {
        tip.t <- y$t
        tip.target <- y$target
        tip.y <- branch.init[tip.target] <- y$y
        for (i in seq_along(tip.t)) {
            idx <- tip.target[i]
            ans <- branches(tip.y[[i]], tip.t[i], pars, 0, idx) # return from combine.branches.quasse()
            lq[idx] <- ans[[1]] # getting the log((sum over all D's) * dt) for this node
            branch.base[[idx]] <- ans[-1]
        }
    }

    ## internal nodes done here
    for (i in order) {
        y.in <- initial.conditions(branch.base[children[i, ]], pars, depth[i], i) # compute initial D for i-th internal node from children
        if (!is.list(y.in) && !all(is.finite(y.in)))
            stop("Bad initial conditions: calculation failure along branches?")
        branch.init[[i]] <- y.in # book keeping, and MLE at each internal node
        ans <- branches(y.in, len[i], pars, depth[i], i) # return from combine.branches.quasse()
        lq[i] <- ans[[1]] # getting the log((sum over all D's) * dt) for this node
        branch.base[[i]] <- ans[-1] # E's and normalized D's
    }

    y.in <- initial.conditions(branch.base[children[root, ]],
        pars, depth[root], root)

    branch.init[[root]] <- y.in

    print(paste0("Log-normalizing factors = ", paste(lq, collapse=", ")))

    ## return with answer we care: lq (normalizing factors for all nodes in the tree) and y.in (normalized D of root)
    list(init = branch.init, base = branch.base, lq = lq, vals = y.in)
}

make.all.branches.quasse.debug <- function (cache, control)  {

    branches <- make.branches.quasse.fftR.debug(control)

    initial.conditions <- make.initial.conditions.quasse.debug(control)

    function(pars, intermediates, preset = NULL) {
        cache$y <- initial.tip.quasse.debug(cache, cache$control, pars[[1]]$x)

        ## prunes tree and returns list with likelihood
        all.branches.list.debug(pars, cache, initial.conditions, branches,
            preset)
    }
}

make.pars.quasse.debug <- function (cache) {
    args <- cache$args
    function(pars) {
        names(pars) <- NULL
        drift <- pars[args$drift]
        diffusion <- pars[args$diffusion]
        ext <- quasse.extent.debug(cache$control, drift, diffusion)
        pars <- diversitree:::expand.pars.quasse(cache$lambda, cache$mu, args,
            ext, pars)
        diversitree:::check.pars.quasse(pars$hi$lambda, pars$hi$mu, drift,
            diffusion)
        pars
    }
}

make.cache.quasse.debug <- function (tree, states, states.sd, lambda, mu, control, sampling.f, for.split = FALSE) {
    tree <- diversitree:::check.tree(tree)
    tmp <- diversitree:::check.states.quasse(tree, states, states.sd)
    states <- tmp$states
    states.sd <- tmp$states.sd
    control <- diversitree:::check.control.quasse(control, tree, states)
    cache <- diversitree:::make.cache(tree)
    cache$states <- states
    cache$states.sd <- states.sd
    cache$control <- control
    if (!for.split) {
        n.lambda <- diversitree:::check.f.quasse(lambda)
        n.mu <- diversitree:::check.f.quasse(mu)
        n.args <- n.lambda + n.mu + 2
        args <- list(lambda = seq_len(n.lambda), mu = seq_len(n.mu) +
            n.lambda, drift = n.lambda + n.mu + 1, diffusion = n.lambda +
            n.mu + 2)
        cache$lambda <- lambda
        cache$mu <- mu
        cache$args <- args
        sampling.f <- diversitree:::check.sampling.f(sampling.f, 1)
        cache$sampling.f <- sampling.f
    }
    cache$info <- diversitree:::make.info.quasse(lambda, mu, tree)
    cache
}

make.quasse.debug <- function (tree, states, states.sd, lambda, mu, control = NULL, sampling.f = NULL) {
    cache <- make.cache.quasse.debug(tree, states, states.sd, lambda, mu, control, sampling.f)

    print("cache$control")
    print(cache$control)

    ## We go function to function to function to function...
    ## make.all.branches.quasse (1) initializes things with make.initial.conditions.quasse(),
    ## and then (2) computes the likelihood with all.branches.list(), which in turn calls make.branches.quasse.fftR()
    ## then make.branches.quasse.fftR() calls combine.branches.quasse()
    ## combine.branches.quasse() makes use of make.pde.quasse.fftR(), which returns a function that actually does the integration
    ## in quasse.integrate.fftR()
    all.branches <- make.all.branches.quasse.debug(cache, cache$control)

    rootfunc <- make.rootfunc.quasse.debug(cache)
    f.pars <- make.pars.quasse.debug(cache)
    ll <- function(pars, condition.surv = TRUE, root = ROOT.OBS,
        root.f = NULL, intermediates = FALSE) {
        pars2 <- f.pars(pars)
        ans <- all.branches(pars2, intermediates) # list containing numbers we care about
        rootfunc(ans, pars2, condition.surv, root, root.f, intermediates)
    }
    class(ll) <- c("quasse", "dtlik", "function")
    ll
}
