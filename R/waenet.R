#' @export
waenet <- function(x, y, pf, adWeight, mids, alpha = 1, nlambda = 100,
                   lambda.min.ratio = 1e-3, lambda = NULL, maxit = 1000,
                   eps = 1e-5)
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

    if (class(mids) != "mids")
        stop("'mids' is not a mice 'mids' object.")

    weight <- (1 - rowMeans(mids$where))
    weight <- rep(weight / m, m)

    if (!is.numeric(nlambda) || length(nlambda) > 1 || nlambda < 1)
        stop("'nlambda' should be an integer >= 1.")

    if (!is.numeric(lambda.min.ratio) ||
        length(lambda.min.ratio) > 1  ||
        lambda.min.ratio < 0)
        stop("'lambda.min.ratio' should be an number >= 0.")

    if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 1)
        stop("'maxit' should be an integer >= 1.")

    if (!is.numeric(eps) || length(eps) > 1 || eps <= 0)
        stop("'eps' should be a postive number.")

    X <- do.call("rbind", x)
    Y <- do.call("c", y)

    X <- scale(X, scale = apply(X, 2, function(.X) sd(.X) * sqrt(m)))

    if (is.null(lambda)) {
        wY_X <- (t(Y * weight) %*% X) / (n * adWeight * alpha) * pf
        lambda.max <- max(abs(wY_X))

        if (all(adWeight == rep(1, p)))
            lambda <- exp(seq(log(lambda.max),
                              log(lambda.max * lambda.min.ratio),
                              length.out = nlambda))
        else
            lambda <- exp(seq(log(lambda.max),
                          log(lambda.max * lambda.min.ratio / 100),
                          length.out = nlambda))
    } else {
        if (!is.numeric(lambda) || !is.vector(lambda))
            stop("'lambda' must be a numeric vector.")
        if (any(lambda <= 0))
            stop("every 'lambda' must be positive.")
        nlambda <- length(lambda)
    }

    df   <- matrix(0, nlambda, length(alpha))
    dev  <- matrix(0, nlambda, length(alpha))
    beta <- array(0, c(nlambda, length(alpha), p + 1))
    for (j in seq(alpha)) {

    for (i in seq(nlambda)) {
            L2 <- lambda[i] * (1 - alpha[j]) * pf
            L1 <- lambda[i] * alpha[j] * adWeight * pf
            fit <- fit.waenet.binomial(X, Y, n, p, m, weight, L1, L2, maxit, eps)
            beta[i, j, ] <- fit$coef
            dev[i, j]    <- fit$dev
            df[i, j]     <- sum(beta[i, j,] != 0)
        }
    }
    structure(list(beta = beta, dev = dev, lambda = lambda, df = df), class = "waenet")
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
    coef[1]  <- beta0 - sum(mu / sd * beta) # intercept
    coef[-1] <- beta / sd

    eta <- X %*% beta + beta0
    dev <- -2 * m * mean(weight * (Y * eta - log(1 + exp(eta))))
    list(coef = coef, dev = dev)
}
