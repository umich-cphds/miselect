#' Cross Validated Multiple Imputation Stacked Adaptive Elastic Net
#'
#' Does k-fold cross-validation for \code{saenet}, and returns optimal values
#' for lambda and alpha.
#'
#' \code{cv.saenet} works by stacking the multiply imputed data into a single
#' matrix and running a weighted adaptive elastic net on it. Simulations suggest
#' that the "stacked" objective function approaches tend to be more
#' computationally efficient and have better estimation and selection
#' properties.
#' @param x A length \code{m} list of \code{n * p} numeric matrices. No matrix
#'     should contain an intercept, or any missing values
#' @param y A length \code{m} list of length \code{n} numeric response vectors.
#'     No vector should contain missing values
#' @param pf Penalty factor of length \code{p}. Can be used to differentially
#'     penalize certain variables. 0 indicates to not penalize the covariate
#' @param adWeight Numeric vector of length p representing the adaptive weights
#'     for the L1 penalty
#' @param weights Numeric vector of length n containing the proportion observed
#'     (non-missing) for each row in the un-imputed data.
#' @param family The type of response. "gaussian" implies a continuous response
#'     and "binomial" implies a binary response. Default is "gaussian".
#' @param alpha Elastic net parameter. Can be a vector to cross validate over.
#'     Default is 1
#' @param nlambda Length of automatically generated 'lambda' sequence. If
#'     lambda' is non NULL, 'nlambda' is ignored. Default is 100
#' @param lambda.min.ratio Ratio that determines the minimum value of 'lambda'
#'     when automatically generating a 'lambda' sequence. If 'lambda' is not
#'     NULL, 'lambda.min.ratio' is ignored. Default is 1e-3
#' @param lambda Optional numeric vector of lambdas to fit. If NULL,
#'    \code{galasso} will automatically generate a lambda sequence based off
#'    of \code{nlambda} and code{lambda.min.ratio}. Default is NULL
#' @param nfolds Number of foldid to use for cross validation. Default is 5,
#'     minimum is 3
#' @param foldid an optional length \code{n} vector of values between 1 and
#     ‘nfold’ identifying what fold each observation is in. Default is NULL and
#'     \code{cv.galasso} will automatically generate folds
#' @param maxit Maximum number of iterations to run. Default is 1000
#' @param eps Tolerance for convergence. Default is 1e-5
#' @return An object of type "cv.saenet" with 9 elements:
#' \describe{
#' \item{lambda}{Sequence of lambdas fit.}
#' \item{cvm}{Average cross validation error for each lambda and alpha. For
#'            family = "gaussian", 'cvm' corresponds to mean squared error,
#'            and for binomial 'cvm' corresponds to deviance.}
#' \item{cvse}{Standard error of 'cvm'.}
#' \item{saenet.fit}{A 'saenet' object fit to the full data.}
#' \item{lambda.min}{The lambda value for the model with the minimum cross
#'                   validation error.}
#' \item{lambda.1se}{The lambda value for the  sparsest model within one
#'                   standard error of the minimum cross validation error.}
#' \item{alpha.min}{The alpha value for the model with the minimum cross
#'                   validation error.}
#' \item{alpha.1se}{The alpha value for the  sparsest model within one
#'                   standard error of the minimum cross validation error.}

#' \item{df}{The number of nonzero coefficients for each value of lambda and alpha.}
#' }
#' @references
#' TODO
#' @examples
#' library(miselect)
#' library(mice)
#'
#' set.seed(48109)
#'
#' # Using the mice defaults for sake of example only.
#' mids <- mice(miselect.df, m = 5, printFlag = FALSE)
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
#' weights  <- 1 - rowMeans(is.na(miselect.df))
#' pf       <- rep(1, 20)
#' adWeight <- rep(1, 20)
#'
#' # Since 'Y' is a binary variable, we use 'family = "binomial"'
#' \donttest{
#' fit <- cv.saenet(x, y, pf, adWeight, weights, family = "binomial")
#'
#' # By default 'coef' returns the betas for (lambda.min , alpha.min)
#' coef(fit)
#' }
#'
#' # You can also cross validate over alpha
#' \donttest{
#' fit <- cv.saenet(x, y, pf, adWeight, weights, family = "binomial",
#'                  alpha = c(.5, 1))
#' # Get selected variables from the 1 standard error rule
#' coef(fit, lambda = fit$lambda.1se, alpha = fit$alpha.1se)
#'
#' }
#' @export
cv.saenet <- function(x, y, pf, adWeight, weights, family =
                      c("gaussian", "binomial"), alpha = 1, nlambda = 100,
                      lambda.min.ratio = 1e-3, lambda = NULL, nfolds = 5,
                      foldid = NULL, maxit = 1000, eps = 1e-5)
{

    if (!is.numeric(nfolds) || length(nfolds) > 1)
        stop("'nfolds' should a be single number.")

    if (!is.null(foldid))
        if (!is.numeric(foldid) || length(foldid) != length(y[[1]]))
            stop("'nfolds' should a be single number.")

    fit <- saenet(x, y, pf, adWeight, weights, family, alpha, nlambda,
                  lambda.min.ratio, lambda, maxit, eps)

    n <- length(y[[1]])
    p <- ncol(x[[1]])
    m <- length(x)

    X <- do.call("rbind", x)
    Y <- do.call("c", y)

    weights <- rep(weights / m , m)

    if (!is.null(foldid)) {
        if (!is.numeric(foldid) || !is.vector(foldid) || length(foldid) != n)
            stop("'foldid' must be length n numeric vector.")
        nfolds <- max(foldid)
    } else {

        r     <- n %% nfolds
        q     <- (n - r) / nfolds
        foldid <- c(rep(seq(nfolds), q), seq(r))
        foldid <- sample(foldid, n)
        foldid <- rep(foldid, m)
    }
    if (nfolds < 3)
        stop("'nfolds' must be bigger than 3.")

    lambda  <- fit$lambda
    nlambda <- length(lambda)
    X.scaled <- scale(X, scale = apply(X, 2, function(.X) stats::sd(.X) * sqrt(m)))

    cvm  <- array(0, c(nlambda, length(alpha), nfolds))
    cvse <- matrix(nlambda, length(alpha))
    for (j in seq(nfolds)) {
        Y.train <- Y[foldid != j]
        X.train <- subset_scaled_matrix(X.scaled, foldid != j)
        w.train <- weights[foldid != j]

        X.test  <- X[foldid == j, , drop = F]
        Y.test  <- Y[foldid == j]
        w.test <- weights[foldid == j]

        cv.fit <- switch(match.arg(family),
            gaussian = fit.saenet.gaussian(X.train, Y.train, n, p, m, w.train,
                                           nlambda, lambda, alpha, pf, adWeight,
                                           maxit, eps),
            binomial = fit.saenet.binomial(X.train, Y.train, n, p, m, w.train,
                                           nlambda, lambda, alpha, pf, adWeight,
                                           maxit, eps)
        )

        cvm[,, j] <- cv.saenet.err(cv.fit, X.test, Y.test, w.test, m)
    }

    cvse <- apply(cvm, c(1, 2), stats::sd) / sqrt(nfolds)
    cvm <- apply(cvm, c(1, 2), mean)
    i <- which.min(cvm)

    row <- (i - 1) %% nlambda + 1

    col <- (i - 1) %/% nlambda + 1

    lambda.min <- fit$lambda[row]
    alpha.min  <- alpha[col]

    j <- abs(cvm - min(cvm)) < cvse[i]
    # Inf could cause NaN if df = 0
    i <- which.min(fit$df * ifelse(j, 1, 1e9))

    row <- (i - 1) %% nlambda + 1
    col <- (i - 1) %/% nlambda + 1

    lambda.1se <- fit$lambda[row]
    alpha.1se  <- alpha[col]

    structure(list(lambda = fit$lambda, alpha = alpha, cvm = cvm, cvse = cvse, saenet.fit =
                   fit, lambda.min = lambda.min, alpha.min = alpha.min,
                   lambda.1se = lambda.1se, alpha.1se = alpha.1se,
                   df = fit$df), class = "cv.saenet")
}


cv.saenet.err <- function(cv.fit, X.test, Y.test, w.test, m)
{
    nalpha <- length(cv.fit$alpha)
    nlambda <- length(cv.fit$lambda)
    cvm <- matrix(0, nlambda, nalpha)

    for (j in seq(nlambda)) {
        for (i in seq(nalpha)) {
            beta <- cv.fit$beta[j, i,]
            beta0 <- beta[1]
            beta  <- beta[-1]
            if ("saenet.gaussian" %in% class(cv.fit)) {
                mse <- m * mean((Y.test - X.test %*% beta - beta0) ^ 2 * w.test)
                cvm[j, i] <- mse
            }
            else {
                eta <- X.test %*% beta + beta0
                dev <-  w.test * (Y.test * eta - log(1 + exp(eta)))
                cvm[j, i] <- -2 * m * mean(dev)
            }
        }
    }
    cvm
}
