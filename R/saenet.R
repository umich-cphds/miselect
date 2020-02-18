#' Multiple Imputation Stacked Adaptive Elastic Net
#'
#' Fits an adaptive elastic net for multiply imputed data. The data is stacked
#' and is penalized that each imputation selects the same betas at each value
#' of lambda. "saenet" supports both continuous and binary responses.
#'
#' \code{saenet} works by stacking the multiply imputed data into a single
#' matrix and running a weighted adaptive elastic net on it. Simulations suggest
#' that the "stacked" objective function approaches tend to be more
#' computationally efficient and have better estimation and selection
#' properties.
#' @param x A list of \code{m} \code{n x p} numeric matrices. No matrix should
#'     contain an intercept, or any missing values
#' @param y A list of \code{m} length n numeric response vectors. No vector
#'     should contain missing values
#' @param pf Penalty factor. Can be used to differentially penalize certain
#'     variables
#' @param adWeight Numeric vector of length p representing the adaptive weights
#'     for the L1 penalty
#' @param weights Numeric vector of length n containing the proportion observed
#'     (non-missing) for each row in the un-imputed data.
#' @param family The type of response. "gaussian" implies a continuous response
#'     and "binomial" implies a binary response. Default is "gaussian".
#' @param alpha Elastic net parameter. Can be a vector to cross validate over.
#'     Default is 1
#' @param lambda Optional numeric vector of lambdas to fit. If NULL,
#'    \code{galasso} will automatically generate a lambda sequence based off
#'    of \code{nlambda} and code{lambda.min.ratio}. Default is NULL
#' @param nlambda Length of automatically generated 'lambda' sequence. If
#'     lambda' is non NULL, 'nlambda' is ignored. Default is 100
#' @param lambda.min.ratio Ratio that determines the minimum value of 'lambda'
#'     when automatically generating a 'lambda' sequence. If 'lambda' is not
#'     NULL, 'lambda.min.ratio' is ignored. Default is 1e-3
#' @param maxit Maximum number of iterations to run. Default is 1000
#' @param eps Tolerance for convergence. Default is 1e-5
#' @return
#' An object with type "saenet" and subtype
#' "saenet.gaussian" or saenet.binomial", depending on which family was used.
#' Both subtypes have 4 elements:
#' \describe{
#' \item{lambda}{Sequence of lambda fit.}
#' \item{beta}{nlambda x nalpha x p + 1 tensor representing the estimated betas
#'             at each value of lambda and alpha.}
#' \item{df}{Number of nonzero betas at each value of lambda and alpha.}
#' \item{mse}{For objects with subtype "saenet.gaussian", the training MSE for
#'            each value of lambda and alpha.}
#' \item{dev}{For objects with subtype "saenet.binomial", the training deviance
#'            for each value of lambda and alpha.}
#' }
#' @examples
#' library(mianet)
#' library(mice)
#'
#' mids <- mice(mianet.df, m = 5, printFlag = FALSE)
#' dfs <- lapply(1:5, function(i) complete(mids, action = i))
#'
#' # Generate list of imputed design matrices and imputed responses
#' x <- list()
#' y <- list()
#' for (i in 1:5) {
#'     x[[i]] <- as.matrix(dfs[[i]][, paste0("X", 1:20)])
#'     y[[i]] <- dfs[[i]]$Y
#' }
#'
#' # Calculate observational weights
#' weights  <- 1 - rowMeans(is.na(mianet.df))
#' pf       <- rep(1, 20)
#' adWeight <- rep(1, 20)
#'
#' # Since 'Y' is a binary variable, we use 'family = "binomial"'
#' \donttest{
#' fit <- saenet(x, y, pf, adWeight, weights, family = "binomial")
#' }
#' @references
#' TODO
#' @export
saenet <- function(x, y, pf, adWeight, weights, family = c("gaussian", "binomial"),
                   alpha = 1, nlambda = 100, lambda.min.ratio = 1e-3,
                   lambda = NULL, maxit = 1000, eps = 1e-5)
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

    if (!is.numeric(weights) || !is.vector(weights) || any(weights < 0) ||
        length(y[[1]]) != length(weights))
    {
        stop("'weights' must be a non negative numeric vector of length n.")
    }

    weights <- rep(weights / m, m)

    if (!is.numeric(nlambda) || length(nlambda) > 1 || nlambda < 1)
        stop("'nlambda' should be an integer >= 1.")

    if (!is.numeric(lambda.min.ratio) || length(lambda.min.ratio) > 1  ||
        lambda.min.ratio < 0)
    {
        stop("'lambda.min.ratio' should be a number >= 0.")
    }

    if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 1)
        stop("'maxit' should be an integer >= 1.")

    if (!is.numeric(eps) || length(eps) > 1 || eps <= 0)
        stop("'eps' should be a postive number.")

    X <- do.call("rbind", x)
    Y <- do.call("c", y)

    X <- scale(X, scale = apply(X, 2, function(.X) stats::sd(.X) * sqrt(m)))

    if (is.null(lambda)) {
        wY_X <- (t(Y * weights) %*% X) / (n * adWeight * mean(alpha)) * pf
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

    fit <- switch(match.arg(family),
        gaussian = fit.saenet.gaussian(X, Y, n, p, m, weights, nlambda, lambda,
                                       alpha, pf, adWeight, maxit, eps),
        binomial = fit.saenet.binomial(X, Y, n, p, m, weights, nlambda, lambda,
                                       alpha, pf, adWeight, maxit, eps)
    )

    cn <- colnames(x[[1]])
    if (is.null(cn))
        cn <- paste0("x", seq(p))
    dimnames(fit$beta) <- list(NULL, NULL, c("(Intercept)", cn))

    return(fit)
}

fit.saenet.binomial <- function(X, Y, n, p, m, weights, nlambda, lambda, alpha,
                                pf, adWeight, maxit, eps)
{

    eta <- matrix(0, n * m, 1)
    pi  <- numeric(n * m)
    res <- numeric(n * m)

    fit.model <- function(L1, L2)
    {
        beta  <- rep(0, p)
        beta0 <- 0
        beta.old  <- beta - 1
        beta0.old <- beta0 - 1
        comp.set <- seq(p)

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
            z2 <- sum(weights * hessian)
            z1 <- sum(weights * hessian * res)
            beta0 <- (z1 + beta0 * z2) / z2
            res <- res - (beta0 - beta0.old)

            # update beta
            for (j in comp.set) {
                z2 <- sum(weights * hessian * X[, j] * X[, j])
                z1 <- sum(weights * hessian * X[, j] * res)
                z <- z1 + beta[j] * z2

                beta[j] <- threshold.saenet(z / n, z2 / n, L1[j], L2[j])
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
        dev <- -2 * m * mean(weights * (Y * eta - log(1 + exp(eta))))
        list(coef = coef, dev = dev)
    }

    df   <- matrix(0, nlambda, length(alpha))
    dev  <- matrix(0, nlambda, length(alpha))
    beta <- array(0, c(nlambda, length(alpha), p + 1))
    for (j in seq(length(alpha))) {
        for (i in seq(nlambda)) {
            L1 <- lambda[i] * alpha[j] * adWeight * pf
            L2 <- lambda[i] * (1 - alpha[j]) * pf
            fit <- fit.model(L1, L2)
            beta[i, j, ] <- fit$coef
            dev[i, j]    <- fit$dev
            df[i, j]     <- sum(beta[i, j,] != 0)
        }
    }

    structure(list(lambda = lambda, alpha = alpha, beta = beta, df = df,
                   dev = dev), class = c("saenet.binomial", "saenet"))
}


fit.saenet.gaussian <- function(X, Y, n, p, m, weights, nlambda, lambda, alpha,
                                pf, adWeight, maxit, eps)
{
    res <- numeric(n * m)
    fit.model <- function(L1, L2)
    {
        beta  <- rep(0, p)
        beta0 <- mean(Y)
        beta.old  <- beta - 1
        beta0.old <- beta0 - 1
        comp.set <- seq(p)

        it <- 0
        while (max(abs(beta - beta.old)) > eps && it < maxit) {
            it <- it + 1

            beta.old  <- beta
            beta0.old <- beta0

            res <- Y - X %*% beta - beta0

            # update b
            for (j in comp.set) {
                z2 <- sum(weights * X[, j] * X[, j])
                z1 <- sum(weights * X[, j] * res)
                z  <- z1 + beta.old[j] * z2

                # soft threshold beta[j]
                beta[j] <- threshold.saenet(z / n, z2 / n, L1[j], L2[j])

                # update residuals
                res <- res - X[, j] * (beta[j] - beta.old[j])
            }
            comp.set <- which(beta != 0)
        }

        # transform back into scale
        mu <- attr(X, "scaled:center")
        sd <- attr(X, "scaled:scale")

        coef <- rep(0, p + 1)

        coef[1]  <- beta0 - sum(mu / sd * beta) # intercept
        coef[-1] <- beta / sd

        mse <- m * mean((Y - X %*% beta - beta0) ^ 2 * weights)
        list(coef = coef, mse = mse)
    }

    df   <- matrix(0, nlambda, length(alpha))
    mse  <- matrix(0, nlambda, length(alpha))
    beta <- array(0, c(nlambda, length(alpha), p + 1))
    for (j in seq(length(alpha))) {
        for (i in seq(nlambda)) {
            L1 <- lambda[i] * alpha[j] * adWeight * pf
            L2 <- lambda[i] * (1 - alpha[j]) * pf

            fit <- fit.model(L1, L2)
            beta[i, j, ] <- fit$coef
            mse[i, j]    <- fit$mse
            df[i, j]     <- sum(beta[i, j,] != 0)
        }
    }

    structure(list(lambda = lambda, alpha = alpha, beta = beta, df = df,
                   mse = mse), class = c("saenet.gaussian", "saenet"))
}
