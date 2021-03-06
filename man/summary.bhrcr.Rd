% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary.bhrcr}
\alias{summary.bhrcr}
\title{Summary Statistics for the Bayesian Clearance Estimator}
\usage{
\method{summary}{bhrcr}(object, ...)
}
\arguments{
\item{object}{an object of class \code{bhrcr}, given by \code{clearanceEstimatorBayes}.}

\item{...}{additional parameters.}
}
\value{
the input object is returned silently.
}
\description{
\code{summary.bhrcr} provides summary statistics for the effect of covariates on both log clearance rates 
and log half-lives.
}
\details{
This function provides mean, median, and credible intervals for gamma, which represents the effect of covariates 
on log clearance rates. It also provides those statistics for the effect of covariates on
log half-lives.
}
\examples{
data("posterior")
summary(posterior)
\dontshow{
data("pursat")
data("pursat_covariates")
data = pursat[pursat["id"] <= 80 & pursat["id"] > 70,]
covariates = pursat_covariates[71:80,]
out <- clearanceEstimatorBayes(data = data, covariates = covariates, outlier.detect = TRUE,
                              niteration = 3, burnin = 1, thin = 1)
summary(out)
}
\donttest{
data("pursat")
data("pursat_covariates")
out <- clearanceEstimatorBayes(data = pursat, covariates = pursat_covariates, 
                               niteration = 200, burnin = 50, thin = 10)
summary(out)
}
}
\author{
Colin B. Fogarty <cfogarty@mit.edu>, Saeed Sharifi-Malvajerdi <saeedsh@wharton.upenn.edu>, Feiyu Zhu <feiyuzhu@sas.upenn.edu>
}
