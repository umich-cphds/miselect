% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coef.R
\name{coef.cv.galasso}
\alias{coef.cv.galasso}
\title{Extract Coefficients From a "cv.galasso" Object}
\usage{
\method{coef}{cv.galasso}(object, lambda = object$lambda.min, ...)
}
\arguments{
\item{object}{A "cv.galasso" fit}

\item{lambda}{Chosen value of lambda. Must be between "min(lambda)" and
"max(lambda)". Default is "lambda.min"}

\item{...}{Additional unused arguments}
}
\value{
A list of numeric vectors containing the coefficients from running
    \code{galasso} on \code{lambda} for each imputation.
}
\description{
Extract Coefficients From a "cv.galasso" Object
}
