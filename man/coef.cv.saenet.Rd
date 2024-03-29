% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R
\name{coef.cv.saenet}
\alias{coef.cv.saenet}
\title{Extract Coefficients From a "cv.saenet" Object}
\usage{
\method{coef}{cv.saenet}(object, lambda = object$lambda.min, alpha = object$alpha.min, ...)
}
\arguments{
\item{object}{A "cv.saenet" fit}

\item{lambda}{Chosen value of lambda. Must be between "min(lambda)" and
"max(lambda)". Default is "lambda.min"}

\item{alpha}{Chosen value of alpha. Must be between "min(alpha)" and
"max(alpha)". Default is "alpha.min"}

\item{...}{Additional unused arguments}
}
\value{
A numeric vector containing the coefficients from running
    \code{saenet} on \code{lambda} and \code{alpha}.
}
\description{
Extract Coefficients From a "cv.saenet" Object
}
