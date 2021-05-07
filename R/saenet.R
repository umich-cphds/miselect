#' Multiple Imputation Stacked Adaptive Elastic Net
#'
#' Fits an adaptive elastic net for multiply imputed data. The data is stacked
#' and is penalized that each imputation selects the same betas at each value
#' of lambda. "saenet" supports both continuous and binary responses.
#'
#' \code{saenet} works by stacking the multiply imputed data into a single
#' matrix and running a weighted adaptive elastic net on it. The objective
#' function is:
#' \deqn{ argmin_{\beta_j} -\frac{1}{n} \sum_{k=1}^{m} \sum_{i=1}^{n} o_i * L(\beta_j|Y_{ik},X_{ijk})}
#' \deqn{ + \lambda (\alpha \sum_{j=1}^{p} \hat{a}_j * pf_j |\beta_{j}|}
#' \deqn{+ (1 - \alpha)\sum_{j=1}^{p} pf_j * \beta_{j}^2)}
#' Where L is the log likelihood, \code{o = w / m}, \code{a} is the
#' adaptive weights, and \code{pf} is the penalty factor. Simulations suggest
#' that the "stacked" objective function approach (i.e., \code{saenet}) tends
#' to be more computationally efficient and have better estimation and selection
#' properties. However, the advantage of \code{galasso} is that it allows one
#' to look at the differences between coefficient estimates across imputations.
#' @param x A length \code{m} list of \code{n * p} numeric matrices. No matrix
#'     should contain an intercept, or any missing values
#' @param y A length \code{m} list of length \code{n} numeric response vectors.
#'     No vector should contain missing values
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
#' \item{coef}{nlambda x nalpha x p + 1 tensor representing the estimated betas
#'             at each value of lambda and alpha.}
#' \item{df}{Number of nonzero betas at each value of lambda and alpha.}
#' @examples
#' \donttest{
#' library(miselect)
#' library(mice)
#'
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
#' fit <- saenet(x, y, pf, adWeight, weights, family = "binomial")
#' }
#' @references
#' Variable selection with multiply-imputed datasets: choosing between stacked
#' and grouped methods. Jiacong Du, Jonathan Boss, Peisong Han, Lauren J Beesley,
#' Stephen A Goutman, Stuart Batterman, Eva L Feldman, and Bhramar Mukherjee. 2020.
#' arXiv:2003.07398
#' @export
saenet <- function(x, y, pf, adWeight, weights, family = c("gaussian", "binomial"),
                   alpha = 1, nlambda = 100, lambda.min.ratio =
                     ifelse(isTRUE(all.equal(adWeight, rep(1, p))), 1e-3, 1e-6),
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
  
  if (!is.numeric(pf) || !is.vector(pf) || any(pf < 0) || length(pf) != p)
  {
    stop("'pf' must be a non negative numeric vector of length p.")
  }
  
  if (!is.numeric(adWeight) || !is.vector(adWeight) || any(adWeight < 0) ||
      length(adWeight) != p)
  {
    stop("'adWeight' must be a non negative numeric vector of length p.")
  }
  
  if (!is.numeric(weights) || !is.vector(weights) || any(weights < 0) ||
      length(weights) != n)
  {
    stop("'weights' must be a non negative numeric vector of length n.")
  }
  
  weights <- rep(weights / m, m)
  
  if (!is.numeric(nlambda) || length(nlambda) > 1 || nlambda < 1)
    stop("'nlambda' should be an integer >= 1.")
  
  if (!is.numeric(lambda.min.ratio) || length(lambda.min.ratio) > 1  ||
      lambda.min.ratio < 0)
  {
    print(lambda.min.ratio)
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
    wY_X <- NULL
    if (match.arg(family) == "gaussian") {
      wY_X <- apply(X, 2, function(x) sum(x * Y * weights))
    }
    else {
      mu <- mean(Y)
      wY_X <- t(ifelse(Y == 1, 1 - mu, - mu) * weights) %*% X
    }
    l <- abs(wY_X / (n * adWeight * max(min(alpha), 0.01) * pf))
    lambda.max <- max(l[is.finite(l)])
    
    lambda <- exp(seq(log(lambda.max),
                      log(lambda.max * lambda.min.ratio),
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
  dimnames(fit$coef) <- list(NULL, NULL, c("(Intercept)", cn))
  
  return(fit)
}

fit.saenet.binomial <- function(X, Y, n, p, m, weights, nlambda, lambda, alpha,
                                pf, adWeight, maxit, eps)
{
  
  eta <- matrix(0, n * m, 1)
  pi  <- numeric(n * m)
  res <- numeric(n * m)
  
  WX2 <- apply(X, 2, function(x) x * x * weights)
  WX  <- apply(X, 2, function(x) x * weights)
  
  fit.model <- function(L1, L2, start)
  {
    beta0 = start[1]
    beta = start[-1]
    beta.old  <- beta - 1
    beta0.old <- beta0 - 1
    comp.set <- seq(p)
    
    it <- 0
    while (max(abs(beta - beta.old)) > eps && it < maxit) {
      it <- it + 1
      
      beta.old  <- beta
      beta0.old <- beta0
      
      eta <- X %*% beta + beta0
      eta  <- exp(eta) / (1 + exp(eta))
      hessian <- eta * (1 - eta)
      res     <- (Y - eta) / hessian
      
      # update beta0
      z2 <- sum(weights * hessian)
      z1 <- sum(weights * hessian * res)
      beta0 <- (z1 + beta0 * z2) / z2
      res <- res - (beta0 - beta0.old)
      
      # update beta
      for (j in comp.set) {
        z2 <- sum(hessian * WX2[, j])
        z1 <- sum(hessian * WX[, j] * res)
        z <- z1 + beta[j] * z2
        
        beta[j] <- S(z / n, L1[j]) / (z2 / n + 2 * L2[j])
        
        res <- res - X[, j] * (beta[j] - beta.old[j])
      }
    }
    
    mu <- attr(X, "scaled:center")
    sd <- attr(X, "scaled:scale")
    
    # transform back into scale
    coef <- rep(0, p + 1)
    coef[1]  <- beta0 - sum(mu / sd * beta) # intercept
    coef[-1] <- beta / sd
    
    eta <- X %*% beta + beta0

    list(coef = coef, beta = c(beta0, beta))
  }
  
  df   <- matrix(0, nlambda, length(alpha))
  beta <- array(0, c(nlambda, length(alpha), p + 1))
  coef <- array(0, c(nlambda, length(alpha), p + 1))
  start = rep(0, p + 1)
  for (j in seq(length(alpha))) {
    for (i in seq(nlambda)) {
      L1 <- lambda[i] * alpha[j] * adWeight * pf
      L2 <- lambda[i] * (1 - alpha[j]) * pf
      fit <- fit.model(L1, L2, start = start)
      start = fit$beta
      beta[i, j, ] <- fit$beta
      coef[i, j, ] <- fit$coef
      df[i, j]     <- sum(beta[i, j,] != 0) - 1
    }
  }
  
  structure(list(lambda = lambda, alpha = alpha, coef = coef, df = df), class = c("saenet.binomial", "saenet"))
}


fit.saenet.gaussian <- function(X, Y, n, p, m, weights, nlambda, lambda, alpha,
                                pf, adWeight, maxit, eps)
{
  res <- numeric(n * m)
  
  z2 <- apply(X, 2, function(x) mean(weights * x * x))
  WX <- apply(X, 2, function(x) weights * x)
  
  fit.model <- function(L1, L2, start)
  {
    beta  <- start[-1]
    beta0 <- start[1]
    beta.old  <- beta - 1
    beta0.old <- beta0 - 1
    comp.set <- seq(p)
    
    it <- 0
    while (max(abs(beta - beta.old)) > eps && it < maxit) {
      it <- it + 1
      beta.old  <- beta
      beta0.old <- beta0
      
      res <- Y - X %*% beta - beta0
      
      # soft threshold active betas
      for (j in comp.set) {
        z1 <- mean(WX[, j] * res)
        z  <- z1 + beta[j] * z2[j]
        
        beta[j] <- S(z, L1[j]) / (z2[j] + 2 * L2[j])
        
        # update residuals
        res <- res - X[, j] * (beta[j] - beta.old[j])
      }
    }
    
    # transform back into scale
    mu <- attr(X, "scaled:center")
    sd <- attr(X, "scaled:scale")
    
    coef <- rep(0, p + 1)
    
    coef[1]  <- beta0 - sum(mu / sd * beta) # intercept
    coef[-1] <- beta / sd
    
    list(coef = coef, beta = c(beta0, beta))
  }
  
  df   <- matrix(0, nlambda, length(alpha))
  beta <- coef <- array(0, c(nlambda, length(alpha), p + 1))
  
  start = rep(0, p + 1)
  for (j in seq(length(alpha))) {
    for (i in seq(nlambda)) {
      L1 <- lambda[i] * alpha[j] * adWeight * pf
      L2 <- lambda[i] * (1 - alpha[j]) * pf
      
      fit <- fit.model(L1, L2, start = start)
      beta[i, j, ] <- fit$beta
      coef[i, j, ] <- fit$coef
      df[i, j]     <- sum(beta[i, j, ] != 0) - 1
    }
  }
  
  structure(list(lambda = lambda, alpha = alpha, coef = coef, df = df), class = c("saenet.gaussian", "saenet"))
}
