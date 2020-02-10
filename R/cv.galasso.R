#' Cross Validated Multiple Imputation Grouped Adaptive LASSO
#'
#' Does k-fold cross-validation for galasso, and returns an optimal value for
#' lambda.
#'
#' TODO
#' @param x A list of \code{m} \code{n x p} numeric matrices. No matrix should
#'     contain an intercept, or any missing values
#' @param y A list of \code{m} length n numeric response vectors. No vector
#'     should contain missing values
#' @param pf Penalty factor. Can be used to differentially penalize certain
#'     variables
#' @param adWeight Numeric vector of length p representing the adaptive weights
#'     for the L1 penalty
#' @param family The type of response. "gaussian" implies a continuous response
#'     and "binomial" implies a binary response. Default is "gaussian".
#' @param nlambda Length of automatically generated 'lambda' sequence. If
#'     lambda' is non NULL, 'nlambda' is ignored. Default is 100
#' @param lambda.min.ratio Ratio that determines the minimum value of 'lambda'
#'     when automatically generating a 'lambda' sequence. If 'lambda' is not
#'     NULL, 'lambda.min.ratio' is ignored. Default is 1e-4
#' @param lambda Optional numeric vector of lambdas to fit. If NULL,
#'    \code{galasso} will automatically generate a lambda sequence based off
#'    of \code{nlambda} and code{lambda.min.ratio}. Default is NULL
#' @param nfolds Number of foldid to use for cross validation. Default is 10,
#'     minimum is 3
#' @param foldid an optional vector of values between 1 and ‘nfold’
#' identifying what fold each observation is in. Default is NULL and
#' \code{cv.galasso} will automatically generate folds
#' @param maxit Maximum number of iterations to run. Default is 1000
#' @param eps Tolerance for convergence. Default is 1e-5
#' @return An object of type "cv.galasso" with 7 elements:
#' \describe{
#' \item{lambda}{Sequence of lambdas fit.}
#' \item{cvm}{Average cross validation error for each 'lambda'. For
#'            family = "gaussian", 'cvm' corresponds to mean squared error,
#'            and for binomial 'cvm' corresponds to deviance.}
#' \item{cvse}{Standard error of 'cvm'.}
#' \item{galasso.fit}{A 'galasso' object fit to the full data.}
#' \item{lambda.min}{The lambda value for the model with the minimum cross
#'                   validation error.}
#' \item{lambda.1se}{The lambda value for the  sparsest model within one
#'                   standard error of the minimum cross validation error.}
#' \item{df}{The number of nonzero coefficients for each value of lambda.}
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
#' pf       <- rep(1, 20)
#' adWeight <- rep(1, 20)

#' # Since 'Y' is a binary variable, we use 'family = "binomial"'
#' \donttest{
#' fit <- cv.galasso(x,y, pf, adWeight, family = "binomial")
#'
#' # Get selected variables from the 1 standard error rule
#' i <- which(fit$lambda == fit$lambda.1se)
#' coef <- fit$galasso.fit$beta[i,]
#' }
#' @references
#' TODO
#' @export
cv.galasso <- function(x, y, pf, adWeight, family = c("gaussian", "binomial"),
                       nlambda = 100, lambda.min.ratio = 1e-4, lambda = NULL,
                       nfolds = 10, foldid = NULL, maxit = 1000, eps = 1e-5)
{

    if (!is.numeric(nfolds) || length(nfolds) > 1)
        stop("'nfolds' should a be single number.")

    if (!is.null(foldid))
        if (!is.numeric(foldid) || length(foldid) != length(y[[1]]))
            stop("'nfolds' should a be single number.")

    fit <- galasso(x, y, pf, adWeight, family, nlambda, lambda.min.ratio,
                   lambda, maxit, eps)

    n <- length(y[[1]])
    if (!is.null(foldid)) {
        if (!is.numeric(foldid) || !is.vector(foldid) || length(foldid) != n)
            stop("'foldid' must be length n numeric vector.")
        nfolds <- max(foldid)
    } else {
        r     <- n %% nfolds
        q     <- (n - r) / nfolds
        folds <- c(rep(seq(nfolds), q), seq(r))
        folds <- sample(folds, n)
    }
    if (nfolds < 3)
        stop("'nfolds' must be bigger than 3.")

    m <- length(x)
    p <- ncol(x[[1]])
    x.scaled <- lapply(x, scale)

    cvm  <- matrix(0, nlambda, nfolds)
    cvse <- numeric(nlambda)
    for (j in seq(nfolds)) {
        x.test  <- lapply(x, function(.x) .x[folds == j, , drop = F])
        y.test  <- lapply(y, function(.y) .y[folds == j])

        y.train <- lapply(y, function(.y) .y[folds != j])
        x.train <- lapply(x.scaled, function(.x)
                              subset_scaled_matrix(.x, folds != j))

        cv.fit <- switch(match.arg(family),
            binomial = fit.galasso.binomial(x.train, y.train, fit$lambda, pf,
                                            adWeight, maxit, eps),
            gaussian = fit.galasso.gaussian(x.train, y.train, fit$lambda, pf,
                                            adWeight, maxit, eps)
        )

        cvm[, j] <- switch(match.arg(family),
            binomial = cv.galasso.err.binomial(cv.fit, x.test, y.test),
            gaussian = cv.galasso.err.gaussian(cv.fit, x.test, y.test)
        )

    }
    cvse <- apply(cvm, 1, stats::sd) / sqrt(nfolds)
    cvm  <- rowMeans(cvm)

    i <- which.min(cvm)
    lambda.min <- fit$lambda[i]
    j <- which((abs(cvm - cvm[i]) < cvse[i]))
    i <- which.min(fit$df[j])
    lambda.1se <- fit$lambda[j][i]

    structure(list(lambda = fit$lambda, cvm = cvm, cvse = cvse, galasso.fit =
                   fit, lambda.min = lambda.min, lambda.1se = lambda.1se, df =
                   fit$df), class = "cv.galasso")
}

# cv.err.gaussian calculates the cross validation error for the gaussian family
# via MSE
cv.galasso.err.gaussian <- function(cv.fit, x.test, y.test)
{
    m       <- length(x.test)
    nlambda <- length(cv.fit$lambda)

    cvm <- numeric(nlambda)
    mse <- rep(0, m)
    for (j in seq(nlambda)) {
        for (i in seq(m)) {
            res <- y.test[[i]] - x.test[[i]] %*% cv.fit$beta[-1, i, j]

            mse[i] <- mean((res - cv.fit$beta[1, i, j]) ^ 2)
            cvm[j] <- mean(mse)
        }
    }
    cvm
}

# cv.err.binomial calculates the cross validation error for the binomial family
# via deviance
cv.galasso.err.binomial <- function(cv.fit, x.test, y.test)
{
    m       <- length(x.test)
    nlambda <- length(cv.fit$lambda)

    cvm <- numeric(nlambda)
    dev <- rep(0, m)
    for (j in seq(nlambda)) {
        for (i in seq(m)) {
            eta <- x.test[[i]] %*% cv.fit$beta[-1, i, j] + cv.fit$beta[1, i, j]
            dev[i] <- -2 * mean(y.test[[i]] * eta - log(1 + exp(eta)))
            cvm[j] <- mean(dev)
        }
    }
    cvm
}
