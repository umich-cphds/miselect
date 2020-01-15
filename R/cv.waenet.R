#' TODO
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
#' @param nfolds Number of foldid to use for cross validation. Default is 10
#' @param maxit Maximum number of iterations to run. Default is 1000
#' @param eps Tolerance for convergence. Default is 1e-5
#' @return TODO
#' @export
cv.waenet <- function(x, y, pf, adWeight, mids, alpha = 1, nlambda = 100,
                      lambda.min.ratio = 1e-3, lambda = NULL, nfolds = 10,
                      foldid = NULL, maxit = 1000, eps = 1e-5)
{

    if (!is.numeric(nfolds) || length(nfolds) > 1)
        stop("'nfolds' should a be single number.")

    if (!is.null(foldid))
        if (!is.numeric(foldid) || length(foldid) != length(y[[1]]))
            stop("'nfolds' should a be single number.")

    fit <- waenet(x, y, pf, adWeight, mids, alpha, nlambda, lambda.min.ratio,
                  lambda, maxit, eps)

    n <- length(y[[1]])
    p <- ncol(x[[1]])
    m <- length(x)

    X <- do.call("rbind", x)
    Y <- do.call("c", y)

    weight <- (1 - rowMeans(mids$where))
    weight <- rep(weight / m , m)

    if (!is.null(foldid)) {
        stop("Not implemented")
    } else {
        r     <- n %% nfolds
        q     <- (n - r) / nfolds
        foldid <- c(rep(seq(nfolds), q), seq(r))
        foldid <- sample(foldid, n)
        foldid <- rep(foldid, m)
    }

    lambda <- fit$lambda
    X.scaled <- scale(X, scale = apply(X, 2, function(.X) sd(.X) * sqrt(m)))
    cvm  <- numeric(nlambda)
    cvse <- numeric(nlambda)
    for (i in seq(nlambda)) {
        L2 <- lambda[i] * (1 - alpha) * pf
        L1 <- lambda[i] * alpha * adWeight * pf
        cv.dev <- numeric(nfolds)
        for (j in seq(nfolds)) {
            Y.train <- Y[foldid != j]
            X.train <- subset_scaled_matrix(X.scaled, foldid != j)
            w.train <- weight[foldid != j]

            X.test  <- X[foldid == j, , drop = F]
            Y.test  <- Y[foldid == j]
            w.test <- weight[foldid == j]

            cv.fit <- fit.waenet.binomial(X.train, Y.train, length(w.train) / m,
                                          p, m, w.train, L1, L2, maxit, eps)

            eta <- X.test %*% cv.fit$coef[-1] + cv.fit$coef[1]
            loglik <- mean(w.test * (Y.test * eta - log(1 + exp(eta))))
            cv.dev[j] <- -2 * m * loglik
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
                   fit$df), class = "cv.waenet")
}
