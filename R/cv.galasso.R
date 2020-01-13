#' @export
cv.galasso <- function(x, y, pf, adWeight, lambda = NULL, nlambda = 100,
                    lambda.min.ratio = 1e-3, nfolds = 10, foldid = NULL,
                    maxit = 1000, eps = 1e-5)
{

    if (!is.numeric(nfolds) || length(nfolds) > 1)
        stop("'nfolds' should a be single number.")

    if (!is.null(foldid))
        if (!is.numeric(foldid) || length(foldid) != length(y[[1]]))
            stop("'nfolds' should a be single number.")


    fit <- galasso(x, y, pf, adWeight, lambda, nlambda, lambda.min.ratio, maxit,
                   eps)

    n <- length(y[[1]])

    if (!is.null(foldid)) {
        stop("Not implemented")
    } else {
        r     <- n %% nfolds
        p     <- (n - r) / nfolds
        folds <- c(rep(1:nfolds, p), 1:r)
        folds <- sample(folds, n)
    }

    x.scaled <- lapply(x, scale)
    m <- length(x)
    cvm  <- numeric(nlambda)
    cvsd <- numeric(nlambda)
    for (i in seq(nlambda)) {
        L <- fit$lambda[i] * adWeight * pf
        cv.dev <- numeric(nfolds)
        for (j in seq(nfolds)) {
            x.test  <- lapply(x, function(.x) .x[folds == j, , drop = F])
            y.test  <- lapply(y, function(.y) .y[folds == j])

            y.train <- lapply(y, function(.y) .y[folds != j])
            x.train <- lapply(x.scaled, function(.x)
                              subset_scaled(.x, folds != j))


            cv.fit <- fit.galasso.binomial(x.train, y.train, L, maxit, eps)

            dev <- rep(0, m)
            for (k in seq(m)) {
                eta <- x.test[[k]] %*% cv.fit$coef[-1, k] + cv.fit$coef[1, k]
                dev[k] <- -2 * mean(y.test[[k]] * eta - log(1 + exp(eta)))
            }
            cv.dev[j] <- mean(dev)
        }
        cvm[i]  <- mean(cv.dev)
        cvsd[i] <- sqrt(sum((cv.dev - cvm[i]) ^ 2) / nfolds)
    }
    i <- which.min(cvm)
    lambda.min <- fit$lambda[i]
    j <- which((abs(cvm - cvm[i]) < cvsd[i]))
    i <- which.min(fit$df[j])
    lambda.1se <- fit$lambda[j][i]

    structure(list(lambda = fit$lambda, cvm = cvm, cvsd = cvsd, galasso.fit =
                   fit, lambda.min = lambda.min, lambda.1se = lambda.1se, df =
                   fit$df), class = "cv.galasso")
}

# Overcomes a technical limitation of R that drops almost all attributes from
# objects when subsetted using `[`
subset_scaled <- function(x, rows)
{
    .x <- x[rows, , drop = F]
    attr(.x, "scaled:center") <- attr(x, "scaled:center")
    attr(.x, "scaled:scale")  <- attr(x, "scaled:scale")
    .x
}
