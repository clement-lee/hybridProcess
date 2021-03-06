% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{llik_nhpp_hybrid}
\alias{llik_nhpp_hybrid}
\title{Evaluate the log-likelihood under the hybrid process}
\usage{
llik_nhpp_hybrid(par, x)
}
\arguments{
\item{par}{vector of length 3, containing (in that order) the power parameter lambda, the exponent parameter theta, and the scale parameter.}

\item{x}{vector of cumulative/absolute event times.}
}
\description{
Evaluate the log-likelihood under the hybrid process
}
