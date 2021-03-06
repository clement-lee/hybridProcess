% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{llik_nhpp_power}
\alias{llik_nhpp_power}
\title{Evaluate the log-likelihood under the power law process, equivalently a special case of hybrid process when the exponent parameter is 0}
\usage{
llik_nhpp_power(par, x)
}
\arguments{
\item{par}{vector of length 2, containing (in that order) the power parameter lambda, and the scale parameter.}

\item{x}{vector of cumulative/absolute event times.}
}
\description{
Evaluate the log-likelihood under the power law process, equivalently a special case of hybrid process when the exponent parameter is 0
}
