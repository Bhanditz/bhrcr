#' @title Summary Statistics for the Bayesian Clearance Estimator
#'
#' @description
#' \code{summary.bhrcr} provides summary statistics for the effect of covariates on both log clearance rates 
#' and log half-lives.
#' 
#' @param object an object of class \code{bhrcr}, given by \code{clearanceEstimatorBayes}.
#' @param ... additional parameters.
#' 
#' @return the input object is returned silently.
#' 
#' @details
#' This function provides mean, median, and credible intervals for gamma, which represents the effect of covariates 
#' on log clearance rates. It also provides those statistics for the effect of covariates on
#' log half-lives.
#'
#' @method summary bhrcr
#' @export
#'
#' @author Colin B. Fogarty <cfogarty@mit.edu>, Saeed Sharifi-Malvajerdi <saeedsh@wharton.upenn.edu>, Feiyu Zhu <feiyuzhu@sas.upenn.edu>
#' @examples
#' data("posterior")
#' summary(posterior)
#' \dontshow{
#' data("pursat")
#' data("pursat_covariates")
#' data = pursat[pursat["id"] <= 80 & pursat["id"] > 70,]
#' covariates = pursat_covariates[71:80,]
#' out <- clearanceEstimatorBayes(data = data, covariates = covariates, outlier.detect = TRUE,
#'                               niteration = 3, burnin = 1, thin = 1)
#' summary(out)
#' }
#' \donttest{
#' data("pursat")
#' data("pursat_covariates")
#' out <- clearanceEstimatorBayes(data = pursat, covariates = pursat_covariates, 
#'                                niteration = 200, burnin = 50, thin = 10)
#' summary(out)
#' }

summary.bhrcr <- function(object, ...) {
  
  # check input class
  if (!inherits(object, "bhrcr")) 
    stop("Object must be of class 'bhrcr'")
  
  cat("\n")
  print(object$CALL)
  cat("\n")
  
  cat("Posterior Estimates and Intervals for the Effect of Covariates on log Clearance Rates \n\n")
  
  table1 <- data.frame(round(object$gamma.mean,4), round(object$gamma.median,4), round(object$gamma.CI[,-2],4))
  colnames(table1) <- c("Mean ", "Median ", colnames(object$gamma.CI[,-2]))
  print(table1)
  
  cat("---\n")
  
  cat("Posterior Estimates and Intervals for the Effect of Covariates on log half-lives \n\n")
  table2 <- data.frame(round(object$halflifeslope.mean,4), round(object$halflifeslope.median,4), round(object$halflifeslope.CI[,-2],4))
  colnames(table2) <- c("Mean ", "Median ", colnames(object$halflifeslope.CI[,-2]))
  print(table2)
  
  cat("---\n")
  cat(paste("Detect Limit: ", object$detect.limit, ", Log Base: ", round(object$log.base,3), "\n"))
} 