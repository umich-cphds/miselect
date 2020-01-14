#' Cross Validated Multiple Imputation Grouped Adaptive LASSO
#' @param x A list of \code{m} \code{n x p} numeric matrices. No matrix should
#'     contain an intercept, or any missing values
#' @param y A list of \code{m} length n numeric response vectors. No vector
#'     should contain missing values
#' @param pf Penalty factor. TODO
#' @param adWeight TODO
#' @param lambda Optional numeric vector of lambdas to fit. If NULL,
#'    \code{galasso} will automatically generate a lambda sequence based off
#'    of \code{nlambda} and code{lambda.min.ratio}. Default is NULL
#' @param nlambda Length of automatically generated 'lambda' sequence. If
#'     lambda' is non NULL, 'nlambda' is ignored. Default is 100
#' @param lambda.min.ratio Ratio that determines the minimum value of 'lambda'
#'     when automatically generating a 'lambda' sequence. If 'lambda' is not
#'     NULL, 'lambda.min.ratio' is ignored. Default is 1e-3
#' @param nfolds Number of folds to use for cross validation. Default is 10
#' @param maxit Maximum number of iterations to run. Default is 1000
#' @param eps Tolerance for convergence. Default is 1e-5
#' @return TODO
#' @export
cv.galasso <- function(x, y, pf, adWeight, nlambda = 100, lambda.min.ratio =
                       1e-4, lambda = NULL, nfolds = 10, foldid = NULL, maxit =
                       1000, eps = 1e-5)
{

    if (!is.numeric(nfolds) || length(nfolds) > 1)
        stop("'nfolds' should a be single number.")

    if (!is.null(foldid))
        if (!is.numeric(foldid) || length(foldid) != length(y[[1]]))
            stop("'nfolds' should a be single number.")


    fit <- galasso(x, y, pf, adWeight, nlambda, lambda.min.ratio, lambda, maxit,
                   eps)

    n <- length(y[[1]])

    if (!is.null(foldid)) {
        stop("Not implemented")
    } else {
        r     <- n %% nfolds
        q     <- (n - r) / nfolds
        folds <- c(rep(seq(nfolds), q), seq(r))
        folds <- sample(folds, n)
    }

    x.scaled <- lapply(x, scale)
    m <- length(x)
    cvm  <- numeric(nlambda)
    cvse <- numeric(nlambda)
    for (i in seq(nlambda)) {
        L <- fit$lambda[i] * adWeight * pf
        cv.dev <- numeric(nfolds)
        for (j in seq(nfolds)) {
            x.test  <- lapply(x, function(.x) .x[folds == j, , drop = F])
            y.test  <- lapply(y, function(.y) .y[folds == j])

            y.train <- lapply(y, function(.y) .y[folds != j])
            x.train <- lapply(x.scaled, function(.x)
                              subset_scaled_matrix(.x, folds != j))


            cv.fit <- fit.galasso.binomial(x.train, y.train, L, maxit, eps)

            dev <- rep(0, m)
            for (k in seq(m)) {
                eta <- x.test[[k]] %*% cv.fit$coef[-1, k] + cv.fit$coef[1, k]
                dev[k] <- -2 * mean(y.test[[k]] * eta - log(1 + exp(eta)))
            }
            cv.dev[j] <- mean(dev)
        }
        cvm[i]  <- mean(cv.dev)
        cvse[i] <- sd(cv.dev) / sqrt(nfolds)
    }
    i <- which.min(cvm)
    lambda.min <- fit$lambda[i]
    j <- which((abs(cvm - cvm[i]) < cvse[i]))
    i <- which.min(fit$df[j])
    lambda.1se <- fit$lambda[j][i]

    structure(list(lambda = fit$lambda, cvm = cvm, cvse = cvse, galasso.fit =
                   fit, lambda.min = lambda.min, lambda.1se = lambda.1se, df =
                   fit$df), class = "cv.galasso")
}
