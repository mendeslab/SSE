#!/usr/bin/env Rscript
# Independent reference for the constant-rate, fixed-origin coverage-study model.
# Requires a prepared SSE build (ant compile native), BEAST2_HOME, and R packages
# ape, coda, ggplot2, patchwork. Run --help for the four required options.

# Parse the diagnostic's explicit input/output paths; burn-in uses the logged Sample column.
options_from_args <- function(args) {
    usage <- paste("Usage: Rscript check.R --xml FILE --log FILE --burnin ITERATION --output DIR",
                   "Compare analytic birth-death inference with a fixed-tree QuaSSE run.", sep = "\n")
    if (any(args %in% c("-h", "--help"))) { cat(usage, "\n"); return(NULL) }
    if (length(args) != 8 || anyDuplicated(args[c(TRUE, FALSE)]) ||
        !setequal(args[c(TRUE, FALSE)], c("--xml", "--log", "--burnin", "--output"))) stop(usage)
    opt <- as.list(setNames(args[c(FALSE, TRUE)], substring(args[c(TRUE, FALSE)], 3)))
    opt$burnin <- as.numeric(opt$burnin)
    if (!is.finite(opt$burnin) || opt$burnin < 0) stop("Invalid burn-in cutoff")
    opt
}

# For one starting lineage and complete extant sampling, p1(t) is the probability of exactly
# one surviving descendant, and S(t) is survival probability. A branch from a to b contributes
# p1(b)/p1(a), and each observed bifurcation contributes lambda. Multiplying cancels all but
# p1(origin) and one p1 per internal node, including the root; divide once by S(origin).
# This is independent of the trait likelihood. Only parameter-independent labeling factors
# are omitted. Source: Stadler (2010), Sampling-through-time in birth-death trees.
# Work with d=lambda-mu>0 and r=mu/lambda in [0,1); log1p avoids cancellation near r=0.
tree_loglik <- function(d, r, ages, origin) {
    log_p1 <- function(t) 2 * log1p(-r) - d * t - 2 * log1p(-r * exp(-d * t))
    ans <- length(ages) * (log(d) - log1p(-r)) + log_p1(origin)
    for (age in ages) ans <- ans + log_p1(age)
    ans - (log1p(-r) - log1p(-r * exp(-d * origin)))
}

# This grid integrates in u=log(d), v=logit(r). Its Jacobian is d*r*(1-r).
# Trapezoidal endpoint weights complete the quadrature; uniform du*dv cancels in normalization.
posterior_grid <- function(n, tail, p, ages, origin) {
    bounds <- qlnorm(c(tail, 1 - tail), p$d_meanlog, p$d_sdlog)
    u <- seq(log(bounds[1]), log(bounds[2]), length.out = n)
    v <- seq(qlogis(qbeta(tail, p$r_alpha, p$r_beta)),
             qlogis(qbeta(1 - tail, p$r_alpha, p$r_beta)), length.out = n)
    z <- expand.grid(d = exp(u), r = plogis(v))
    lp <- dlnorm(z$d, p$d_meanlog, p$d_sdlog, log = TRUE) +
        dbeta(z$r, p$r_alpha, p$r_beta, log = TRUE) + log(z$d) + log(z$r) + log1p(-z$r)
    lp <- lp + tree_loglik(z$d, z$r, ages, origin)
    ends <- c(0.5, rep(1, n - 2), 0.5)
    w <- exp(lp - max(lp)) * as.vector(outer(ends, ends))
    z$weight <- w / sum(w)
    z$lambda <- z$d / (1 - z$r)
    z$mu <- z$r * z$lambda
    z
}

# Quadrature weights approximate a continuous distribution. Merge ties, then interpolate
# through the middle of each weight's CDF interval to avoid grid-sized quantile jumps.
weighted_quantiles <- function(x, w, probabilities = c(.025, .5, .975)) {
    a <- aggregate(w, list(x = x), sum)
    a <- a[order(a$x), ]
    cdf <- cumsum(a[, 2]) - a[, 2] / 2
    # Negligible tail weights can round to identical cumulative probabilities.
    keep <- !duplicated(cdf)
    approx(cdf[keep], a$x[keep], xout = probabilities, rule = 2)$y
}

# Report moments and equal-tailed quantiles on each biological parameter's original scale.
grid_summary <- function(z) {
    do.call(rbind, lapply(c("d", "r", "lambda", "mu"), function(name) {
        x <- z[[name]]; mean <- sum(x * z$weight)
        q <- weighted_quantiles(x, z$weight)
        data.frame(parameter = name, mean = mean, variance = sum((x - mean)^2 * z$weight),
                   lower = q[1], median = q[2], upper = q[3])
    }))
}

# Independent lognormal/Beta moments also give the moments of lambda=d/(1-r) and mu=d*r/(1-r).
prior_sds <- function(p) {
    a <- p$r_alpha; b <- p$r_beta
    if (b <= 2) stop("This diagnostic's prior-SD convergence criterion needs turnover beta > 2")
    ed <- exp(p$d_meanlog + p$d_sdlog^2 / 2); ed2 <- exp(2 * p$d_meanlog + 2 * p$d_sdlog^2)
    sqrt(c(d = ed2 - ed^2, r = a * b / ((a + b)^2 * (a + b + 1)),
           lambda = ed2 * (a + b - 1) * (a + b - 2) / ((b - 1) * (b - 2)) -
               (ed * (a + b - 1) / (b - 1))^2,
           mu = ed2 * a * (a + 1) / ((b - 1) * (b - 2)) - (ed * a / (b - 1))^2))
}

# Refinement and wider bounds must both stabilize summaries before the grid is used as a reference.
converged <- function(old, new, prior_sd) {
    changes <- abs(as.matrix(old[c("mean", "lower", "median", "upper")]) -
                   as.matrix(new[c("mean", "lower", "median", "upper")])) / prior_sd
    variance_change <- abs(old$variance / new$variance - 1)
    message(sprintf("Grid check: max location change / prior SD %.4g; max variance change %.4g",
                    max(changes), max(variance_change)))
    all(changes < .005) && all(variance_change < .01)
}

# Compare empirical and reference CDFs without density-smoothing choices. Joint contours are
# relative density levels (not credible regions). The lambda/mu density includes Jacobian 1/lambda.
comparison_plots <- function(z, draws, p, ages, origin) {
    labels <- c(d = "Diversification", r = "Turnover", lambda = "Speciation (lambda)", mu = "Extinction (mu)")
    cdfs <- lapply(names(labels), function(name) {
        a <- aggregate(z$weight, list(value = z[[name]]), sum)
        a <- a[order(a$value), ]; a$cdf <- cumsum(a[, 2])
        # Display at most 2000 points; full-resolution weights are retained for all calculations.
        a <- a[unique(round(seq(1, nrow(a), length.out = min(2000, nrow(a))))), ]
        ggplot2::ggplot() +
            ggplot2::geom_line(data = a, ggplot2::aes(value, cdf, colour = "Analytic")) +
            ggplot2::stat_ecdf(data = draws, ggplot2::aes(x = .data[[name]], colour = "BEAST")) +
            ggplot2::labs(title = labels[[name]], x = "Parameter value", y = "Cumulative probability") +
            ggplot2::scale_colour_manual(values = c(Analytic = "#0072B2", BEAST = "#C62828"), name = NULL) +
            ggplot2::theme_bw()
    })
    joints <- lapply(list(c("d", "r"), c("lambda", "mu")), function(names) {
        bounds <- lapply(names, function(name) weighted_quantiles(z[[name]], z$weight, c(.001, .999)))
        g <- expand.grid(x = seq(bounds[[1]][1], bounds[[1]][2], length.out = 151),
                         y = seq(bounds[[2]][1], bounds[[2]][2], length.out = 151))
        transformed <- names[1] == "lambda"
        d <- if (transformed) g$x - g$y else g$x
        r <- if (transformed) g$y / g$x else g$y
        valid <- d > 0 & r > 0 & r < 1
        g$log_density <- NA_real_
        g$log_density[valid] <- tree_loglik(d[valid], r[valid], ages, origin) +
            dlnorm(d[valid], p$d_meanlog, p$d_sdlog, log = TRUE) +
            dbeta(r[valid], p$r_alpha, p$r_beta, log = TRUE) -
            if (transformed) log(g$x[valid]) else 0
        g$relative_density <- exp(g$log_density - max(g$log_density, na.rm = TRUE))
        ggplot2::ggplot() +
            ggplot2::geom_point(data = draws, ggplot2::aes(x = .data[[names[1]]], y = .data[[names[2]]]),
                                alpha = .08, size = .5) +
            ggplot2::geom_contour(data = g, ggplot2::aes(x, y, z = relative_density),
                                  breaks = exp(-c(.5, 2, 4.5)), colour = "#0072B2") +
            ggplot2::labs(x = labels[[names[1]]], y = labels[[names[2]]]) + ggplot2::theme_bw()
    })
    list(patchwork::wrap_plots(cdfs, ncol = 2) + patchwork::plot_annotation(title = "Analytic and BEAST posteriors"),
         patchwork::wrap_plots(joints, ncol = 2) + patchwork::plot_annotation(
             title = "Joint posteriors: BEAST samples and analytic density contours",
             caption = "Contours: exp(-0.5), exp(-2), exp(-4.5) times the displayed grid's maximum density."))
}

# Run the standalone Java driver with explicit paths and capture its verbose initialization output.
run_command <- function(command, args, logfile) {
    status <- system2(command, shQuote(args), stdout = logfile, stderr = logfile)
    if (status != 0) stop(command, " failed; see ", logfile)
}

# Compile only the diagnostic driver, export/validate inputs, integrate, then compare production
# likelihoods and the existing chain. No model source, input XML, or MCMC log is modified.
main <- function() {
    opt <- options_from_args(commandArgs(trailingOnly = TRUE))
    if (is.null(opt)) return(invisible(NULL))
    for (pkg in c("ape", "coda", "ggplot2", "patchwork"))
        if (!requireNamespace(pkg, quietly = TRUE)) stop("Missing R package: ", pkg)
    script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)))
    root <- normalizePath(file.path(dirname(script), "../.."))
    beast <- Sys.getenv("BEAST2_HOME")
    if (!dir.exists(beast)) stop("Set BEAST2_HOME to the BEAST installation used for the prepared SSE build")
    dir.create(opt$output, recursive = TRUE, showWarnings = FALSE)
    out <- normalizePath(opt$output)
    build <- tempfile("birth-death-reference-", tmpdir = file.path(root, "build"))
    dir.create(build)
    on.exit(unlink(build, recursive = TRUE), add = TRUE)
    cp <- paste(c(build, file.path(root, c("build", "lib/*")),
                  file.path(beast, c("lib/packages/BEAST.base.jar", "lib/packages/BEAST.app.jar", "lib/launcher.jar"))),
                collapse = .Platform$path.sep)
    java_home <- Sys.getenv("JAVA_HOME")
    java <- if (nzchar(java_home)) file.path(java_home, "bin/java") else "java"
    javac <- if (nzchar(java_home)) file.path(java_home, "bin/javac") else "javac"
    run_command(javac, c("-cp", cp, "-d", build, file.path(dirname(script), "BirthDeathReference.java")),
                file.path(out, "compile.log"))
    jargs <- c("--enable-native-access=ALL-UNNAMED", "-Xmx3g", paste0("-Djava.library.path=", root, "/build/native"),
               "-cp", cp, "BirthDeathReference")
    run_command(java, c(jargs, "inputs", normalizePath(opt$xml), out, root), file.path(out, "inputs.log"))
    metadata <- read.delim(file.path(out, "inputs.tsv"))
    p <- as.list(setNames(metadata$value, metadata$name))
    if (any(!is.finite(unlist(p))) || min(p$d_sdlog, p$r_alpha, p$r_beta, p$sd_sdlog, p$origin) <= 0)
        stop("Invalid origin/prior parameters")
    tree <- ape::read.tree(file.path(out, "tree.nwk"))
    if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0) ||
        !ape::is.ultrametric(tree) || !ape::is.binary(tree) || !ape::is.rooted(tree))
        stop("Expected a rooted, bifurcating contemporaneous tree with finite nonnegative branch lengths")
    depth <- ape::node.depth.edgelength(tree)
    height <- max(depth[seq_along(tree$tip.label)])
    ages <- height - depth[length(tree$tip.label) + seq_len(tree$Nnode)]
    if (p$origin < height) stop("Origin is younger than the root")
    message(length(tree$tip.label), " tips; root age ", height, "; origin ", p$origin)
    prior_sd <- prior_sds(p)
    z <- posterior_grid(201, 1e-6, p, ages, p$origin); previous <- grid_summary(z)
    stable <- FALSE
    for (n in c(401, 801)) {
        z <- posterior_grid(n, 1e-6, p, ages, p$origin); current <- grid_summary(z)
        stable <- converged(previous, current, prior_sd)
        if (stable) break
        previous <- current
    }
    if (!stable) stop("Posterior grid has not converged at 801 points per dimension")
    expanded <- posterior_grid(n, 1e-8, p, ages, p$origin)
    if (!converged(current, grid_summary(expanded), prior_sd)) stop("Posterior bounds have not stabilized")
    z <- expanded; reference <- grid_summary(z)
    reference$prior_sd <- prior_sd
    write.table(reference, file.path(out, "reference.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    message("Reference posterior converged with ", n, " points per dimension and expanded bounds.")

    points <- expand.grid(d = qlnorm(c(.1, .5, .9), p$d_meanlog, p$d_sdlog),
                          r = qbeta(c(.1, .5, .9), p$r_alpha, p$r_beta),
                          sd = qlnorm(c(.25, .75), p$sd_meanlog, p$sd_sdlog))
    configs <- data.frame(config = c("original", "space", "time", "both", "wide"),
                          resolution = c(1, 2, 1, 2, 2), time = c(1, 1, .5, .5, .5), width = c(1, 1, 1, 1, 2))
    requests <- do.call(rbind, lapply(seq_len(nrow(configs)), function(i)
        cbind(point = seq_len(nrow(points)), config = configs$config[i], points,
              configs[rep(i, nrow(points)), c("resolution", "time", "width")])))
    write.table(requests, file.path(out, "points.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    message("Evaluating 18 parameter points at each of five numerical configurations...")
    run_command(java, c(jargs, "evaluate", normalizePath(opt$xml), file.path(out, "points.tsv"),
                       file.path(out, "quasse.tsv"), root, beast), file.path(out, "quasse.log"))
    likelihoods <- read.delim(file.path(out, "quasse.tsv"))
    likelihoods$analytic <- tree_loglik(likelihoods$d, likelihoods$r, ages, p$origin)
    likelihoods$error <- NA_real_
    for (group in split(seq_len(nrow(likelihoods)), interaction(likelihoods$config, likelihoods$sd))) {
        baseline <- group[likelihoods$point[group] %in% c(5, 14)]
        stopifnot(length(baseline) == 1)
        likelihoods$error[group] <- (likelihoods$logL[group] - likelihoods$logL[baseline]) -
            (likelihoods$analytic[group] - likelihoods$analytic[baseline])
    }
    write.table(likelihoods, file.path(out, "likelihood-comparison.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    print(aggregate(abs(error) ~ config, likelihoods, max))
    if (any(abs(likelihoods$error) > .001)) message("Some likelihood differences exceed the 0.001 target; inspect refinement.")

    draws <- read.table(opt$log, header = TRUE, comment.char = "#", check.names = FALSE)
    if (!all(c("Sample", "diversification", "turnover", "lambda", "mu") %in% names(draws))) stop("Missing log columns")
    if (any(diff(draws$Sample) <= 0)) stop("Sample must increase strictly")
    draws <- draws[draws$Sample > opt$burnin, ]
    if (nrow(draws) < 20 || any(diff(draws$Sample) != diff(draws$Sample)[1])) stop("Need >=20 equally spaced samples")
    draws$d <- draws$diversification; draws$r <- draws$turnover
    if (any(!is.finite(as.matrix(draws[c("d", "r", "lambda", "mu")])))) stop("Nonfinite log values")
    if (any(draws$d <= 0 | draws$r < 0 | draws$r >= 1) ||
        !isTRUE(all.equal(draws$lambda, draws$d / (1 - draws$r))) ||
        !isTRUE(all.equal(draws$mu, draws$r * draws$lambda))) stop("Logged rate parameterization does not match the XML")
    chain_summary <- do.call(rbind, lapply(c("d", "r", "lambda", "mu"), function(name) {
        x <- draws[[name]]; q <- quantile(x, c(.025, .5, .975), type = 1)
        halves <- split(x, seq_along(x) > floor(length(x) / 2))
        data.frame(parameter = name, mean = mean(x), variance = var(x), lower = q[1], median = q[2], upper = q[3],
                   ess = as.numeric(coda::effectiveSize(x)), first_mean = mean(halves[[1]]), second_mean = mean(halves[[2]]),
                   first_median = median(halves[[1]]), second_median = median(halves[[2]]),
                   first_variance = var(halves[[1]]), second_variance = var(halves[[2]]))
    }))
    chain_summary$mean_mcse <- sqrt(chain_summary$variance / chain_summary$ess)
    chain_summary$mean_difference <- chain_summary$mean - reference$mean
    chain_summary$variance_ratio <- chain_summary$variance / reference$variance
    write.table(chain_summary, file.path(out, "chain.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    print(reference); print(chain_summary)
    plots <- comparison_plots(z, draws, p, ages, p$origin)
    grDevices::pdf(file.path(out, "posterior-comparison.pdf"), width = 11, height = 8, onefile = TRUE)
    tryCatch(for (plot in plots) print(plot), finally = grDevices::dev.off())
    message("Wrote diagnostic tables and posterior-comparison.pdf to ", out)
}

if (sys.nframe() == 0L) main()
