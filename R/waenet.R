#' @export
waenet <- function(x, y, pf, weight, adWeight, nlambda = 100,
                   lambda.min.ratio = 1e-3, lambda = NULL, alpha = 1,
                   maxit = 1000, eps = 1e-5)
{
    if (!is.list(x))
        stop("'x' should be a list of numeric matrices.")
    if (any(sapply(x, function(.x) !is.matrix(.x) || !is.numeric(.x))))
        stop("Every 'x' should be a numeric matrix.")

    n <- nrow(x[[1]])
    p <- ncol(x[[1]])
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

    X <- do.call("rbind", x)
    Y <- do.call("c", y)
    X <- scale(X, scale = apply(X, 2, function(.X) sd(.X) * sqrt(m)))

    v <- log(p) / log(n * m)
    gamma <- ceiling(2 * v / (1 - v)) + 1
    adWeight.power <- (gamma + 1) / 2

    # adaptive ENet update
    if (is.null(lambda)) {
        wY_X <- t(Y * weight) %*% X / (n * adWeight * alpha) * pf
        lambda.max <- max(abs(wY_X))

        if (all(adWeight == rep(1, p)))
            lambda <- exp(seq(log(lambda.max),
                              log(lambda.max * lambda.min.ratio),
                              length.out = nlambda))
        else
            lambda <- exp(seq(log(lambda.max * (n * D) ^ (-adWeight.power)),
                              log(lambda.max * (n * D) ^ (-0.5)),
                              length.out = nlambda))
    }

    df   <- numeric(nlambda)
    dev  <- numeric(nlambda)
    beta <- matrix(0, nlambda, p + 1)
    for (i in seq(nlambda)) {
        L2 <- lambda[i] * (1 - alpha) * pf
        L1 <- lambda[i] * alpha * adWeight * pf

        fit <- fit.waenet.binomial(X, Y, n, p, m, weight, L1, L2, maxit, eps)
        beta[i, ] <- fit$coef
        dev[i] <- fit$dev[i]
        df[i] <- sum(beta[i,] != 0)
    }
    structure(list(beta = beta, dev = dev, lambda = lambda, df = df),
              class = "waenet")
}

fit.waenet.binomial <- function(X, Y, n, p, m, weight, L1, L2, maxit, eps)
{
    beta  <- rep(0, p)
    beta0 <- 0
    beta.old  <- beta - 1
    beta0.old <- beta0 - 1
    comp.set <- seq(p)

    eta <- matrix(0, n * m, 1)
    pi  <- numeric(n * m)
    res <- numeric(n * m)

    it <- 0
    while (max(abs(beta - beta.old)) > eps && it < maxit) {
        it <- it + 1

        beta.old  <- beta
        beta0.old <- beta0

        eta <- X %*% beta + beta0
        pi <- exp(eta) / (1 + exp(eta))
        pi[abs(1 - pi) <= 1e-5] <- 1
        pi[abs(pi) <= 1e-5] <- 0
        hessian <- pi * (1 - pi)
        hessian[hessian <= 1e-5] <- 1e-5
        res <- (Y - pi) / hessian

        # update beta0
        z2 <- sum(weight * hessian)
        z1 <- sum(weight * hessian * res)
        beta0 <- (z1 + beta0 * z2) / z2
        res <- res - (beta0 - beta0.old)

        # update beta
        for (j in comp.set) {
            z2 <- sum(weight * hessian * X[, j] * X[, j])
            z1 <- sum(weight * hessian * X[, j] * res)
            z <- z1 + beta[j] * z2

            beta[j] <- threshold.waenet(z / n, z2 / n, L1[j], L2[j])
            res <- res - X[, j] * (beta[j] - beta.old[j])
        }
        comp.set <- which(beta != 0)
    }

    mu <- attr(X, "scaled:center")
    sd <- attr(X, "scaled:scale")
    # transform back into scale
    coef <- rep(0, p + 1)
    coef[1]  <- beta0 - sum(mu / sd * beta.naive) # intercept
    coef[-1] <- beta / sd

    eta <- X %*% beta + beta0
    dev <- -2 * mean(Y * eta - log(1 + exp(eta)))
    list(coef = coef, dev = dev)
}
