#' Multiple Imputation Grouped Adaptive LASSO
#'
#' \code{galasso} fits an adaptive LASSO for multiply imputed data. "galasso"
#' supports both continuous and binary responses.
#'
#' \code{galasso} works by adding a group penalty to the aggregated objective
#' function to ensure selection consistency across imputations. The objective
#' function is:
#'
#' \deqn{argmin_{\beta_{jk}} - L(\beta_{jk}| X_{ijk}, Y_{ik})}
#' \deqn{+ \lambda * \Sigma_{j=1}^{p} \hat{a}_j * pf_j * \sqrt{\Sigma_{k=1}^{m} \beta_{jk}^2}}
#' Where L is the log likelihood,\code{a} is the adaptive weights, and
#' \code{pf} is the penalty factor. Simulations suggest that the "stacked"
#' objective function approach (i.e., \code{saenet}) tends to be more
#' computationally efficient and have better estimation and selection
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
#' @param family The type of response. "gaussian" implies a continuous response
#'     and "binomial" implies a binary response. Default is "gaussian".
#' @param lambda Optional numeric vector of lambdas to fit. If NULL,
#'    \code{galasso} will automatically generate a lambda sequence based off
#'    of \code{nlambda} and code{lambda.min.ratio}. Default is NULL
#' @param nlambda Length of automatically generated 'lambda' sequence. If
#'     lambda' is non NULL, 'nlambda' is ignored. Default is 100
#' @param lambda.min.ratio Ratio that determines the minimum value of 'lambda'
#'     when automatically generating a 'lambda' sequence. If 'lambda' is not
#'     NULL, 'lambda.min.ratio' is ignored. Default is 1e-4
#' @param maxit Maximum number of iterations to run. Default is 10000
#' @param eps Tolerance for convergence. Default is 1e-5
#' @return
#' An object with type "galasso" and subtype
#' "galasso.gaussian" or galasso.binomial", depending on which family was used.
#' Both subtypes have 4 elements:
#' \describe{
#' \item{lambda}{Sequence of lambda fit.}
#' 
#' \item{coef}{a list of length D containing the coefficient estimates from running 
#' galasso at each value of lambda. Each element in the list is a nlambda x (p+1) matrix.}
#' \item{df}{Number of nonzero betas at each value of lambda.}
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
#' pf       <- rep(1, 20)
#' adWeight <- rep(1, 20)
#'
#' fit <- galasso(x, y, pf, adWeight)
#' }
#' @references
#' Variable selection with multiply-imputed datasets: choosing between stacked
#' and grouped methods. Jiacong Du, Jonathan Boss, Peisong Han, Lauren J Beesley,
#' Stephen A Goutman, Stuart Batterman, Eva L Feldman, and Bhramar Mukherjee. 2020.
#' arXiv:2003.07398
#' @export
galasso <- function(x, y, pf, adWeight, family = c("gaussian", "binomial"),
                    nlambda = 100, lambda.min.ratio =
                      ifelse(isTRUE(all.equal(adWeight, rep(1, p))), 1e-3, 1e-6),
                    lambda = NULL, maxit = 10000, eps = 1e-5)
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
  
  if (!is.numeric(adWeight) || !is.vector(adWeight) ||
      length(adWeight) != p || any(adWeight < 0))
  {
    stop("'adWeight' should be a non negative vector of length p.")
  }
  
  if (!is.numeric(pf) || !is.vector(pf) || length(pf) != p || any(pf < 0))
    stop("'pf' should be a non negative vector of length p.")
  
  if (!is.numeric(nlambda) || length(nlambda) > 1 || nlambda < 1)
    stop("'nlambda' should be an integer >= 1.")
  
  if (!is.numeric(lambda.min.ratio) || length(lambda.min.ratio) > 1  ||
      lambda.min.ratio < 0)
  {
    stop("'lambda.min.ratio' should be an number >= 0.")
  }
  
  if (!is.numeric(maxit) || length(maxit) > 1 || maxit < 1)
    stop("'maxit' should be an integer >= 1.")
  
  if (!is.numeric(eps) || length(eps) > 1 || eps <= 0)
    stop("'eps' should be a postive number.")
  
  x <- lapply(x, function(x) scale(x))
  
  if (is.null(lambda)) {
    Y_X <- matrix(0, m, p)
    if (match.arg(family) == "gaussian") {
      for (i in seq(m))
        Y_X[i,] <- t(y[[i]]) %*% x[[i]]
    }
    else {
      for (i in seq(m)) {
        mu <- mean(y[[i]])
        Y_X[i,] <- t(ifelse(y[[i]] == 1, 1 - mu, - mu)) %*% x[[i]]
      }
    }
    l <- sqrt(apply(Y_X ^ 2, 2, sum)) / (n * adWeight * pf)
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
                gaussian = fit.galasso.gaussian(x, y, lambda, adWeight, pf, maxit, eps),
                binomial = fit.galasso.binomial(x, y, lambda, adWeight, pf, maxit, eps))
  
  cn <- colnames(x[[1]])
  if (is.null(cn))
    cn <- paste0("x", seq(p))
  dimnames(fit$coef)[[1]] <- c("(Intercept)", cn)
  fit$coef <- lapply(1:dim(fit$coef)[2], function(x) { t(fit$coef[,x,]) })
  
  return(fit)
}

fit.galasso.binomial <- function(x, y, lambda, adWeight, pf, maxit, eps)
{
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  m <- length(x)
  
  eta <- matrix(0, n, m)
  pi  <- matrix(0, n, m)
  res <- matrix(0, n, m)
  y.tilde <- y
  fit.model <- function(start, L) {
    beta  <- start[-1,]
    beta0 <- start[1,]
    
    it <- 0
    beta.old  <- beta - 1
    beta0.old <- beta0 - 1
    comp.set  <- seq(p)
    while (max(abs(beta - beta.old)) > eps && it < maxit) {
      it <- it + 1
      beta.old  <- beta
      beta0.old <- beta0
      for (i in seq(m)) {
        eta[, i] <- x[[i]] %*% beta[, i]
        eta[, i] <- eta[, i] + beta0[i]
        
        pi[, i] <- exp(eta[, i]) / (1 + exp(eta[, i]))
        
        y.tilde[[i]] <- eta[, i] + 4 * (y[[i]] - pi[, i])
        res[, i]     <- y.tilde[[i]] - eta[, i]
      }
      
      # update intercept, res
      for (i in seq(m)) {
        beta0[i] <- mean(y.tilde[[i]])
        res[, i] <- res[, i] - (beta0[i] - beta0.old[i])
      }
      
      # soft threshold each beta and update residuals
      for (j in comp.set) {
        z <- rep(0, m)
        for (i in seq(m))
          z[i] <- t(x[[i]][, j]) %*% res[, i] + n * beta[j, i]
        
        #soft threshold beta_j
        normz <- sqrt(sum((z/(4*n))^2))
        
        if(is.na(normz)) {
          beta = start[-1,]
          beta0 = start[1,]
          warning("galasso does not converge for small penalties")
          break
        } else {
          beta[j,] <- t(S(normz, L[j]) * z / (n*normz))
        }
        
        # update residuals with thresholded beta
        for (i in seq(m))
          res[, i] <- res[, i] - x[[i]][, j] * (beta[j, i] - beta.old[j, i])
      }
      
      if(is.na(normz)) {break} # to escape the while loop
    }
    if (max(abs(beta - beta.old)) >= eps)
      warning("galasso did not converge: delta = ", max(abs(beta - beta.old)))
    
    # tranform back into scale
    coef <- matrix(NA, p + 1, m)
    for (i in seq(m)) {
      mu <- attr(x[[i]], "scaled:center")
      sd <- attr(x[[i]], "scaled:scale")
      
      # intercept
      coef[1, i] <- beta0[i] - sum(mu / sd *  beta[, i])
      coef[2:(p + 1), i] <- beta[, i] / sd
    }
    
    list(beta = rbind(beta0, beta), coef = coef)
  } # end of fit.model function
  
  nlambda <- length(lambda)
  beta <- array(0, c(p + 1, m, nlambda))
  coef <- array(0, c(p + 1, m, nlambda))
  df   <- rep(0, nlambda)
  start  <- matrix(0, p + 1, m)
  for (i in seq(nlambda)) {
    
    # avoid optimization for small penalties once divergence happens
    if(i > 1 && all(fit$beta == start)) {
      start       <- fit$beta
      beta[,, i]  <- fit$beta
      coef[,, i]  <- fit$coef
      df[i]       <- sum(beta[-1, 1, i] != 0)
      next
    }
    
    if(i > 1) {start <- fit$beta}
    
    L <- lambda[i] * adWeight * pf
    fit <- fit.model(start, L)
    beta[,, i]  <- fit$beta
    coef[,, i]  <- fit$coef
    df[i]       <- sum(beta[-1, 1, i] != 0)
  }
  
  structure(list(lambda = lambda, coef = coef, df = df),
            class = c("galasso.binomial", "galasso"))
}

fit.galasso.gaussian <- function(x, y, lambda, adWeight, pf, maxit, eps)
{
  n <- nrow(x[[1]])
  p <- ncol(x[[1]])
  m <- length(x)
  
  res <- matrix(0, n, m)
  eta <- matrix(0, n, m)
  fit.model <- function(start, L)
  {
    beta  <- start[-1,]
    beta0 <- start[1,]
    it <- 0
    beta.old  <- beta - 1
    beta0.old <- beta0 - 1
    comp.set  <- seq(p)
    while (max(abs(beta - beta.old)) > eps && it < maxit) {
      it <- it + 1
      beta.old  <- beta
      beta0.old <- beta0
      
      for (i in seq(m)) {
        eta[, i] <- x[[i]] %*% beta[, i]
        # res[, i] <- y[[i]] - eta[, i] - beta0[i]
        res[, i] <- y[[i]] - eta[, i]
      }
      
      # update b
      for (j in comp.set) {
        # update z
        z <- rep(0, m)
        for (i in seq(m))
          z[i] <- t(x[[i]][, j]) %*% res[, i] + n * beta[j, i]
        
        # soft threshold beta_j
        normz    <- sqrt(sum(z ^ 2))
        if(is.na(normz)) {
          beta = start[-1,]
          beta0 = start[1,]
          warnings("galasso goes not converge for small penalties")
          break
        } else {
          beta[j,] <- t(S(normz / n, L[j]) * z / normz)
        }
        
        for (i in seq(m))
          res[, i] <- res[,i] - x[[i]][, j] * (beta[j, i] - beta.old[j, i])
      }
      if(is.na(normz)) {break}
    }
    if (max(abs(beta - beta.old)) >= eps)
      warning("galasso did not converge: delta = ", max(abs(beta - beta.old)))
    
    # tranform back into scale
    coef <- matrix(NA, p + 1, m)
    for (i in seq(m)) {
      mu <- attr(x[[i]], "scaled:center")
      sd <- attr(x[[i]], "scaled:scale")
      
      # intercept
      coef[1, i] <- beta0[i] - sum(mu / sd *  beta[, i])
      coef[2:(p + 1), i] <- beta[, i] / sd
    }
    
    list(beta = rbind(beta0, beta), coef = coef)
  }
  
  nlambda <- length(lambda)
  beta <- coef <- array(0, c(p + 1, m, nlambda))
  df   <- rep(0, nlambda)
  start  <- matrix(1, p + 1, m)
  start[1,] = 0
  for (i in seq(nlambda)) {
    
    if(i > 1 && all(fit$beta == start)) {
      start <- fit$beta
      beta[,,i] <- fit$beta
      coef[,,i] <- fit$coef
      df[i]    <- sum(beta[-1, 1, i] != 0)
      next
    }
    
    if(i > 1) {start <- fit$beta}
    
    L <- lambda[i] * adWeight * pf
    fit <- fit.model(start, L)
    beta[,,i] <- fit$beta
    coef[,,i] <- fit$coef
    df[i]    <- sum(beta[-1, 1, i] != 0)
  }
  structure(list(lambda = lambda, coef = coef, df = df),
            class = c("galasso.gaussian", "galasso"))
}