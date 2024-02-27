#' Cross Validated Multiple Imputation Grouped Adaptive LASSO
#'
#' Does k-fold cross-validation for \code{galasso}, and returns an optimal value
#' for lambda.
#'
#' \code{cv.galasso} works by adding a group penalty to the aggregated objective
#' function to ensure selection consistency across imputations. Simulations
#' suggest that the "stacked" objective function approaches (i.e., \code{saenet})
#' tend to be more computationally efficient and have better estimation and
#' selection properties.
#' @param x A length \code{m} list of \code{n * p} numeric matrices. No matrix
#'     should contain an intercept, or any missing values
#' @param y A length \code{m} list of length \code{n} numeric response vectors.
#'     No vector should contain missing values
#' @param pf Penalty factor. Can be used to differentially penalize certain
#'     variables
#' @param adWeight Numeric vector of length p representing the adaptive weights
#'     for the L1 penalty
#' @param family The type of response. "gaussian" implies a continuous response
#'     and "binomial" implies a binary response. Default is "gaussian".
#' @param nlambda Length of automatically generated "lambda" sequence. If
#'     "lambda" is non NULL, "nlambda" is ignored. Default is 100
#' @param lambda.min.ratio Ratio that determines the minimum value of "lambda"
#'     when automatically generating a "lambda" sequence. If "lambda" is not
#'     NULL, "lambda.min.ratio" is ignored. Default is 1e-4
#' @param lambda Optional numeric vector of lambdas to fit. If NULL,
#'    \code{galasso} will automatically generate a lambda sequence based off
#'    of \code{nlambda} and code{lambda.min.ratio}. Default is NULL
#' @param nfolds Number of foldid to use for cross validation. Default is 5,
#'     minimum is 3
#' @param foldid an optional length \code{n} vector of values between 1 and
#     "nfold" identifying what fold each observation is in. Default is NULL and
#'     \code{cv.galasso} will automatically generate folds
#' @param maxit Maximum number of iterations to run. Default is 10000
#' @param eps Tolerance for convergence. Default is 1e-5
#' @returns An object of type "cv.galasso" with 7 elements:
#' \describe{
#' \item{call}{The call that generated the output.}
#' \item{lambda}{The sequence of lambdas fit.}
#' \item{cvm}{Average cross validation error for each "lambda". For
#'            family = "gaussian", "cvm" corresponds to mean squared error,
#'            and for binomial "cvm" corresponds to deviance.}
#' \item{cvse}{Standard error of "cvm".}
#' \item{galasso.fit}{A "galasso" object fit to the full data.}
#' \item{lambda.min}{The lambda value for the model with the minimum cross
#'                   validation error.}
#' \item{lambda.1se}{The lambda value for the  sparsest model within one
#'                   standard error of the minimum cross validation error.}
#' \item{df}{The number of nonzero coefficients for each value of lambda.}
#' }
#' @examples
#' \donttest{
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
#' pf       <- rep(1, 20)
#' adWeight <- rep(1, 20)
#'
#' fit <- cv.galasso(x, y, pf, adWeight)
#'
#' # By default 'coef' returns the betas for lambda.min.
#' coef(fit)
#' }
#' @references
#' Du, J., Boss, J., Han, P., Beesley, L. J., Kleinsasser, M., Goutman, S. A., ... 
#' & Mukherjee, B. (2022). Variable selection with multiply-imputed datasets: 
#' choosing between stacked and grouped methods. Journal of Computational and 
#' Graphical Statistics, 31(4), 1063-1075. <doi:10.1080/10618600.2022.2035739>
#' @export
cv.galasso <- function(x, y, pf, adWeight, family = c("gaussian", "binomial"),
                       nlambda = 100, lambda.min.ratio =
                         ifelse(isTRUE(all.equal(adWeight, rep(1, p))), 1e-3, 1e-6),
                       lambda = NULL, nfolds = 5, foldid = NULL, maxit = 1000,
                       eps = 1e-5)
{
  call <- match.call()
  
  if (!is.list(x))
    stop("'x' should be a list of numeric matrices.")
  if (any(sapply(x, function(.x) !is.matrix(.x) || !is.numeric(.x))))
    stop("Every 'x' should be a numeric matrix.")
  
  dim <- dim(x[[1]])
  n <- dim[1]
  p <- dim[2]
  m <- length(x)
  
  if (!is.numeric(nfolds) || length(nfolds) > 1)
    stop("'nfolds' should a be single number.")
  
  if (!is.null(foldid))
    if (!is.numeric(foldid) || length(foldid) != length(y[[1]]))
      stop("'nfolds' should a be single number.")
  
  fit <- galasso(x, y, pf, adWeight, family, nlambda, lambda.min.ratio,
                 lambda, maxit, eps)

  
  if (!is.null(foldid)) {
    if (!is.numeric(foldid) || !is.vector(foldid) || length(foldid) != n)
      stop("'foldid' must be length n numeric vector.")
    nfolds <- max(foldid)
  } else {
    r     <- n %% nfolds
    q     <- (n - r) / nfolds
    if(r == 0) {
      folds <- c(rep(seq(nfolds), q))
      folds <- sample(folds, n)
    } else {
      folds <- c(rep(seq(nfolds), q), seq(r))
      folds <- sample(folds, n)
    }
  }
  if (nfolds < 3)
    stop("'nfolds' must be bigger than 3.")
  
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
                     gaussian = fit.galasso.gaussian(x.train, y.train, fit$lambda, pf,
                                                     adWeight, maxit, eps),
                     binomial = fit.galasso.binomial(x.train, y.train, fit$lambda, pf,
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
  range = cvm[i] + cvse[i]
  id.all = which(cvm <= range)
  lambda.1se <- max(fit$lambda[id.all])
  
  structure(list(call = call, lambda = fit$lambda, cvm = cvm, cvse = cvse,
                 galasso.fit = fit, lambda.min = lambda.min, lambda.1se =
                 lambda.1se, df = fit$df), class = "cv.galasso")
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
      res <- y.test[[i]] - x.test[[i]] %*% cv.fit$coef[-1, i, j]
      
      mse[i] <- mean((res - cv.fit$coef[1, i, j]) ^ 2)
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
  for (j in seq(nlambda)) {
    dev <- rep(0, m)
    for (i in seq(m)) {
      eta <- x.test[[i]] %*% cv.fit$coef[-1, i, j] + cv.fit$coef[1, i, j]
      dev[i] <- -2 * mean(y.test[[i]] * eta - log(1 + exp(eta)))
    }
    cvm[j] <- mean(dev)
  }
  cvm
}

#' Print cv.galasso Objects
#'
#' \code{print.cv.galasso} print the fit and returns it invisibly.
#' @param x An object of type "cv.galasso" to print
#' @param ... Further arguments passed to or from other methods
#' @export
print.cv.galasso <- function(x, ...)
{
  nl <- length(x$lambda)
  
  out <- cbind(x$cvm, x$df)
  
  dimnames(out) <- list(paste0("l.", seq(nl)), c("cvm", "df"))
  cat("'cv.galasso' fit:\n")
  print(x$call)
  cat("Average cross validation error and df for each lambda\n")
  print(out)
  cat("lambda min:\n")
  cat(x$lambda.min, "\n", sep = "")
  cat("lambda 1 SE:\n")
  cat(x$lambda.1se, "\n", sep = "")
  invisible(x)
}