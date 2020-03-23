#' Extract Coefficients From a 'cv.saenet' Object
#' @param object A 'cv.saenet' fit
#' @param lambda Chosen value of lambda. Must be between 'min(lambda)' and
#'     'max(lambda)'. Default is 'lambda.min'
#' @param alpha Chosen value of alpha. Must be between 'min(alpha)' and
#'     'max(alpha)'. Default is 'alpha.min'
#' @param ... Additional unused arguments
#' @return A numeric vector containing the coefficients from running
#'     \code{saenet} on \code{lambda} and \code{alpha}.
#' @export
coef.cv.saenet <- function(object, lambda = object$lambda.min,
                           alpha = object$alpha.min, ...)
{
    if (class(object) != "cv.saenet")
        stop("'object' must have class 'cv.saenet'.")

    if (lambda < min(object$lambda) || lambda > max(object$lambda))
        stop("'lambda' must be between 'min(lambda)' and 'max(lambda)'.")
    if (alpha < min(object$alpha) || alpha > max(object$alpha))
        stop("'alpha' must be between 'min(alpha)' and 'max(alpha)'.")

    # Calculate lambda weights (lw) and alpha (wa) weights to attempt to
    # interpolate the coefficients if the user doesn't pick a lambda from the
    # lambda sequence.
    l <- abs(log(lambda) - log(object$lambda))
    l <- (l / (max(l) - min(l))) ^ -2

    lw <- l / sum(l)
    lw <- ifelse(is.nan(lw), 1, lw)

    a <- abs(alpha - object$alpha)
    a <- (a / (max(a) - min(a))) ^ -2

    aw <-  a / sum(a)
    aw <- ifelse(is.nan(aw), 1, aw)

    w <- lw %*% t(aw)
    apply(object$saenet.fit$beta, 3, function(x) sum(w * x))
}

#' Extract Coefficients From a 'saenet' Object
#'
#' \code{coef.galasso} averages the estimates across imputations to return a
#' single vector instead of a matrix.
#' @param object A 'cv.saenet' fit
#' @param lambda Chosen value of lambda. Must be between 'min(lambda)' and
#'     'max(lambda)'. Default is 'lambda.min'
#' @param alpha Chosen value of alpha. Must be between 'min(alpha)' and
#'     'max(alpha)'. Default is 'alpha.min'
#' @param ... Additional unused arguments
#' @return A numeric vector containing the coefficients from running
#'     \code{saenet} on \code{lambda} and \code{alpha}.
#' @export
coef.saenet <- function(object, lambda, alpha, ...)
{
    if (!("saenet" %in% class(object)))
        stop("'object' must have class 'saenet'.")

    if (lambda < min(object$lambda) || lambda > max(object$lambda))
        stop("'lambda' must be between 'min(lambda)' and 'max(lambda)'.")
    if (alpha < min(object$alpha) || alpha > max(object$alpha))
        stop("'alpha' must be between 'min(alpha)' and 'max(alpha)'.")

    # Calculate lambda weights (lw) and alpha (wa) weights to attempt to
    # interpolate the coefficients if the user doesn't pick a lambda from the
    # lambda sequence.
    l <- abs(log(lambda) - log(object$lambda))
    l <- (l / (max(l) - min(l))) ^ -2

    lw <- l / sum(l)
    lw <- ifelse(is.nan(lw), 1, lw)

    a <- abs(alpha - object$alpha)
    a <- (a / (max(a) - min(a))) ^ -2

    aw <-  a / sum(a)
    aw <- ifelse(is.nan(aw), 1, aw)

    w <- lw %*% t(aw)
    apply(object$beta, 3, function(x) sum(w * x))
}

#' Extract Coefficients From a 'cv.galasso' Object
#' @param object A 'cv.galasso' fit
#' @param lambda Chosen value of lambda. Must be between 'min(lambda)' and
#'     'max(lambda)'. Default is 'lambda.min'
#' @param ... Additional unused arguments
#' @return A numeric vector containing the coefficients from running
#'     \code{galasso} on \code{lambda}.
#' @export
coef.cv.galasso <- function(object, lambda = object$lambda.min, ...)
{
    if (class(object) != "cv.galasso")
        stop("'object' must have class 'cv.galasso'.")

    if (lambda < min(object$lambda) || lambda > max(object$lambda))
        stop("'lambda' must be between 'min(lambda)' and 'max(lambda)'.")

    # Calculate lambda weights to attempt to interpolate the coefficients if the
    # user doesn't pick a lambda from the lambda sequence.
    l <- abs(log(lambda) - log(object$lambda))
    l <- (l / (max(l) - min(l))) ^ -2

    w <- l / sum(l)
    w <- ifelse(is.nan(w), 1, w)

    apply(object$galasso.fit$beta, 2, function(x) sum(w * x))
}

#' Extract Coefficients From a 'galasso' Object
#' @param object A 'galasso' fit
#' @param lambda Chosen value of lambda. Must be between 'min(lambda)' and
#'     'max(lambda)'. Default is 'lambda.min'
#' @param ... Additional unused arguments
#' @return A numeric vector containing the coefficients from running
#'     \code{galasso} on \code{lambda}.
#' @export
coef.galasso <- function(object, lambda, ...)
{
    if (!("galasso" %in% class(object)))
        stop("'object' must have class 'galasso'.")

    if (lambda < min(object$lambda) || lambda > max(object$lambda))
        stop("'lambda' must be between 'min(lambda)' and 'max(lambda)'.")

    l <- abs(log(lambda) - log(object$lambda))
    l <- (l / (max(l) - min(l))) ^ -2

    w <- l / sum(l)
    w <- ifelse(is.nan(w), 1, w)

    apply(object$beta, 2, function(x) sum(w * x))
}
