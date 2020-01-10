
#' @param x A list of \code{m} \code{n x p} numeric matrices. No matrix should
#'     contain an intercept.
#' @param y A list of \code{m} length n numeric response vectors.
#' @param pf Penalty factor.
#' @param adWeight TODO
#' @param lambda TODO
#' @param nlambda Length of the generated lambda sequence if 'lambda' was
#'     not specified
#' @param maxit Maximum number of iterations to run. Default is 1000
#' @param eps Tolerance for convergence. Default is 1e-5
#' @export
galasso <- function(x, y, pf, adWeight, lambda = NULL, nlambda = 100,
                    lambda.min.ratio = 1e-3, maxit = 1000, eps = 1e-5)
{
    if (!is.list(x))
        stop("'x' should be a list of numeric matrices.")
    if (any(sapply(x, function(.x) !is.matrix(.x) || !is.numeric(.x))))
        stop("Every 'x' should be a numeric matrix.")

    dim <- dim(x[[1]])
    n <- dim[1]
    p <- dim[2]
    m <- length(x)

    if (any(sapply(x, function(.x) any(dim(x) != dim))))
        stop("Every matrix in 'x' must have the same dimensions.")

    if (!is.list(y))
        stop("'y' should be a list of numeric vectors.")
    if (length(y) != m)
        stop("'y' should should have the same length as 'x'.")
    if (any(sapply(y, function(y) !is.numeric(y) || !is.vector(y))))
        stop("Every 'y' should be a numeric vector.")
    if (any(sapply(y, function(y) !is.numeric(y) || !is.vector(y))))
            stop("Every 'y' should be a numeric vector.")

    x <- lapply(x, function(x) scale(x))

    v <- log(p * m) / log(n * m)
    gamma <- ceiling(2 * v / (1 - v)) + 1
    adWeight.power <- (gamma + 1) / 2

    Y_X <- matrix(0, m, p)
    for (i in 1:m)
        Y_X[i,] <- t(y[[i]]) %*% x[[i]]

    if (is.null(lambda)) {
        norm <- sqrt(apply(Y_X ^ 2, 2, sum))
        lambda.max <- max(pf * norm / (n * adWeight))

        if (all(adWeight == rep(1, p)))
            lambda <- exp(seq(log(lambda.max),
                              log(lambda.max * lambda.min.ratio),
                              length.out = nlambda))
        else
            lambda <- exp(seq(log(lambda.max * (n * D) ^ (-adWeight.power)),
                              log(lambda.max * (n * D) ^ (-0.5)),
                              length.out = nlambda))
    }

    beta <- matrix(0, nlambda, p + 1)
    dev <- rep(0, nlambda)
    df <- rep(0, nlambda)
    for (i in seq(nlambda)) {
        L <- lambda[i] * adWeight * pf
        fit <- fit.galasso.binomial(x, y, L, maxit, eps)

        dev[i]   <- mean(fit$dev)
        beta[i,] <- rowMeans(fit$coef)
        df[i]    <- sum(beta[i,] != 0)
    }


    structure(list(beta = beta, dev = dev, lambda = lambda, df = df),
              class = "galasso")
}


fit.galasso.binomial <- function(x, y, L, maxit, eps)
{

    n <- nrow(x[[1]])
    p <- ncol(x[[1]])
    m <- length(x)

    beta  <- matrix(1, p, m)
    beta0 <- rep(1, m)

    eta <- matrix(0, n, m)
    pi  <- matrix(0, n, m)
    res <- matrix(0, n, m)
    y.tilde <- y

    it <- 0
    beta.old  <- beta - 1
    beta0.old <- beta0 - 1
    comp.set  <- 1:p
    while (max(abs(beta - beta.old)) > eps && it < maxit) {

        it <- it + 1
        beta.old  <- beta
        beta0.old <- beta0
        for (i in 1:m) {
            eta[, i] <- x[[i]] %*% beta[, i]
            eta[, i] <- eta[, i] + beta0[i]

            pi[, i] <- exp(eta[, i]) / (1 + exp(eta[, i]))
            pi[abs(1 - pi[, i]) <= 1e-5, i] <- 1
            pi[abs(pi[, i]) <= 1e-5, i]     <- 0

            y.tilde[[i]] <- eta[, i] + 4 * (y[[i]] - pi[, i])
            res[, i]     <- y.tilde[[i]] - eta[, i]
        }

        # update intercept, res
        for (i in 1:m) {
            beta0[i] <- mean(y.tilde[[i]])
            res[, i] <- res[, i] - (beta0[i] - beta0.old[i])
        }

        # soft threshold each beta and update residuals
        for (j in comp.set) {

            z <- rep(0, m)
            for (i in 1:m)
                z[i] <- t(x[[i]][, j]) %*% res[, i] + n * beta[j, i]

            #soft threshold beta
            beta[j,] <- t(threshold.gaenet(z / (4 * n), L[j], 0))

            # update res with thresholded beta
            for (i in 1:m)
                res[, i] <- res[, i] - x[[i]][, j] * (beta[j, i] - beta.old[j, i])
        }

        comp.set <- which(beta[, 1] != 0)
    }
    if (max(abs(beta - beta.old)) >= eps)
        warning("gaenet did not converge.", max(abs(beta - beta.old)))

    # tranform back into scale
    coef.naive <- matrix(NA, p + 1, m)
    for (i in 1:m) {
        mu <- attr(x[[i]], "scaled:center")
        sd <- attr(x[[i]], "scaled:scale")

        coef.naive[1, i] <- beta0[i] - sum(mu / sd * beta[, i])
        coef.naive[-1, i] <- beta[, i] / sd
    }

    # tranform back into scale
    coef <- matrix(NA, p + 1, m)
    for (i in 1:m) {
        mu <- attr(x[[i]], "scaled:center")
        sd <- attr(x[[i]], "scaled:scale")

        # intercept
        coef[1, i] <- beta0[i] - sum(mu / sd *  beta[, i])
        coef[2:(p + 1), i] <- beta[, i] / sd
    }

    dev <- rep(0, m)
    for (i in seq(m)) {
        eta <- x[[i]] %*% beta[, i] + beta0[i]
        dev[i] <- -2 * mean(y[[i]] * eta - log(1 + exp(eta)))
    }

    list(coef.naive = coef.naive, coef = coef, dev = dev)
}
