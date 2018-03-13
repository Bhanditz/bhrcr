#' @title Print Function for the Bayesian Clearance Estimator
#'
#' @description
#' \code{print.bhrcr} prints the estimated effect of covariates on both log clearance rates and
#' log half-lives.
#' 
#' @param x an object of class \code{bhrcr}, given by \code{clearanceEstimatorBayes}.
#' @param ... additional parameters.
#' 
#' @return the input object is returned silently.
#' 
#' @details
#' This function prints the posterior mean value of gamma, which represents the effect of covariates 
#' on log clearance rates. It also prints the estimated impact of covariates on log half-lives.
#'
#' @method print bhrcr
#' @export
#' 
#' @author Colin B. Fogarty <cfogarty@mit.edu>, Saeed Sharifi-Malvajerdi <saeedsh@wharton.upenn.edu>, Feiyu Zhu <feiyuzhu@sas.upenn.edu>
#' @examples
#' data("posterior")
#' print(posterior)
#' \dontshow{
#' data("pursat")
#' data("pursat_covariates")
#' data = pursat[pursat["id"] <= 80 & pursat["id"] > 70,]
#' covariates = pursat_covariates[71:80,]
#' out <- clearanceEstimatorBayes(data = data, covariates = covariates, outlier.detect = TRUE,
#'                               niteration = 3, burnin = 1, thin = 1)
#' print(out)
#' }
#' \donttest{
#' data("pursat")
#' data("pursat_covariates")
#' out <- clearanceEstimatorBayes(data = pursat, covariates = pursat_covariates, 
#'                                niteration = 200, burnin = 50, thin = 10)
#' print(out)
#' }

print.bhrcr <- function(x, ...) {
  
  # check input class
  if (!inherits(x, "bhrcr")) 
    stop("Object must be of class 'bhrcr'")
  
  cat("\n")
  print(x$CALL)
  cat("\n")
  
  cat("Estimates (= Posterior Mean) for the Effect of Covariates on log Clearance Rates \n\n")
  
  table <- data.frame(round(x$gamma.mean,4) )
  colnames(table) <- c("Estimate")
  rownames(table) <- rownames(x$gamma.CI)
  print(table)
  
  cat("---\n")
  
  cat("Estimates (= Posterior Mean) for the Effect of Covariates on log half-lives \n\n")
  
  table <- data.frame(round(x$halflifeslope.mean,4) )
  colnames(table) <- c("Estimate")
  rownames(table) <- rownames(x$halflifeslope.CI)
  print(table)
  
  cat("---\n")
  cat(paste("Detect Limit: ", x$detect.limit, ", Log Base: ", round(x$log.base,3), "\n"))
} 