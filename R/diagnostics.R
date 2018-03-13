#' @title Diagnostics Function for MCMC
#'
#' @description
#' \code{diagnostics} provides diagnostic analysis for the MCMC process used in the main function 
#' \code{clearanceEstimatiorBayes}.
#'
#' @param object an object of class \code{bhrcr}, given by \code{clearanceEstimatorBayes}.
#' @param ... additional parameters.
#' 
#' @return the directory location under which all the output is saved.
#' 
#' @details
#' This function provides diagnostic analysis such as trace plots, ACF and PACF plots for some important parameters in the simulation process of Gibbs sampling.
#' With these diagnostic plots, we can be assured that we get the results after we have reached stationarity and have thinned sufficiently.
#'
#' @export
#' 
#' @author Colin B. Fogarty <cfogarty@mit.edu>, Saeed Sharifi-Malvajerdi <saeedsh@wharton.upenn.edu>, Feiyu Zhu <feiyuzhu@sas.upenn.edu>
#' @examples
#' \dontshow{
#' data("pursat")
#' data("pursat_covariates")
#' data = pursat[pursat["id"] <= 80 & pursat["id"] > 70,]
#' covariates = pursat_covariates[71:80,]
#' out <- clearanceEstimatorBayes(data = data, covariates = covariates, outlier.detect = TRUE,
#'                               niteration = 3, burnin = 1, thin = 1)
#' diagnostics(out)
#' }
#' \donttest{
#' data("posterior")
#' diagnostics(posterior)
#' }
#' \donttest{
#' data("pursat")
#' data("pursat_covariates")
#' out <- clearanceEstimatorBayes(data = pursat, covariates = pursat_covariates,
#'                                niteration = 200, burnin = 50, thin = 10)
#' diagnostics(out)
#' }

diagnostics <- function(object, ...) {
  
  # check input class
  if (!inherits(object, "bhrcr")) 
      stop("Object must be of class 'bhrcr'") 
    
  # save the current working directory
  dir.old <- getwd()
  # get current working directory
  dir <- getwd()
  # create a sub-directory for storing diagnostic plots
  if (!dir.exists(paste(dir, "/mcmcDiagnostics", sep=""))){
      dir.create(paste(dir, "/mcmcDiagnostics", sep=""))
  }
  setwd(paste(dir,"/mcmcDiagnostics", sep=""))
  
  # to initialize the sequential plot mechanism
  par(ask = FALSE)
  
  # get the customized inputs niteration and step size to thin the sample
  sampleSize <- length(object$var.error.post)
  nsim       <- length(object$var.epsilon.post)
  niteration <- nsim - object$burnin
  step       <- niteration/sampleSize
  
  # get names of covariates 
  names <- rownames(object$gamma.CI)
  
  # get the x-axis (simulation sequence) for the diagnostic plots
  # x is the whole sequence of simulations including the burn-in period
  x <- seq(1, nsim)
  # x.thin is the sequence after burn-in and thinning
  x.thin <- seq(1, sampleSize)
  
  # trace plots for the whole simulation
  for (i in 1:nrow(object$gamma.post)) {
    filename = paste('trace_of_gamma_', names[i], sep = '') 
    pdf(file = sprintf("%s.pdf", filename))
    plot(x, object$gamma.post[i,], 'l', main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
    dev.off()
  }
  
  filename = paste('trace_of_variance_of_epsilon')
  pdf(file = sprintf("%s.pdf", filename))
  plot(x, object$var.epsilon.post, 'l', main=filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = paste('trace_of_pi_lag')
  pdf(file = sprintf("%s.pdf", filename))
  plot(x, object$p.lag, 'l', main=filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = paste('trace_of_pi_tail')
  pdf(file = sprintf("%s.pdf", filename))
  plot(x, object$p.tail, 'l', main=filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  # autocorrelation diagnostics (for the whole samples)
  for (i in 1:nrow(object$gamma.post)) {
    filename = paste('ACF of gamma_', names[i], sep = '')
    pdf(file = sprintf("%s.pdf", filename))
    acf(object$gamma.post[i,], lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
    dev.off()
    
    filename = paste('PACF of gamma_', names[i], sep = '')
    pdf(file = sprintf("%s.pdf", filename))
    pacf(object$gamma.post[i,], lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
    dev.off()
  }
  
  filename = 'ACF_of_variance_of_epsilon'
  pdf(file = sprintf("%s.pdf", filename))
  acf(object$var.epsilon.post, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'PACF_of_variance_of_epsilon'
  pdf(file = sprintf("%s.pdf", filename))
  pacf(object$var.epsilon.post, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'ACF_of_pi_lag'
  pdf(file = sprintf("%s.pdf", filename))
  acf(object$p.lag, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'PACF_of_pi_lag'
  pdf(file = sprintf("%s.pdf", filename))
  pacf(object$p.lag, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'ACF_of_pi_tail'
  pdf(file = sprintf("%s.pdf", filename))
  acf(object$p.tail, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'PACF_of_pi_tail'
  pdf(file = sprintf("%s.pdf", filename))
  pacf(object$p.tail, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  # trace plots for the samples after burn-in and thinning
  for (i in 1:nrow(object$gamma.post.thin)) {
    filename = paste('trace_of_gamma_', names[i], ' after_thinning', sep = '')
    pdf(file = sprintf("%s.pdf", filename))
    plot(x.thin, object$gamma.post.thin[i,], 'l', main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
    dev.off()
  }
  
  filename = 'trace_of_variance_of_epslion_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  plot(x.thin, object$var.error.post, 'l', main=filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'trace_of_pi_lag_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  plot(x.thin, object$p.lag.thin, 'l', main=filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'trace_of_pi_tail_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  plot(x.thin, object$p.tail.thin, 'l', main=filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  # autocorrelation diagnostics (thinned Gibbs samples)
  for (i in 1:nrow(object$gamma.post.thin)) {
    filename = paste('ACF_of_gamma_', names[i], ' after_thinning', sep = '')
    pdf(file = sprintf("%s.pdf", filename))
    acf(object$gamma.post.thin[i,], lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
    dev.off()
    
    filename = paste('PACF_of_gamma_', names[i], ' after_thinning', sep = '')
    pdf(file = sprintf("%s.pdf", filename))
    pacf(object$gamma.post.thin[i,], lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
    dev.off()
  }
  
  filename = 'ACF_of_variance_of_epsilon_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  acf(object$var.error.post, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'PACF_of_variance_of_epsilon_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  pacf(object$var.error.post, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'ACF_of_pi_lag_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  acf(object$p.lag.thin, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'PACF_of_pi_lag_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  pacf(object$p.lag.thin, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'ACF_of_pi_tail_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  acf(object$p.tail.thin, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()
  
  filename = 'PACF_of_pi_tail_after_thinning'
  pdf(file = sprintf("%s.pdf", filename))
  pacf(object$p.tail.thin, lag.max = step, main = filename, sub=NULL, xlab = 'Simulations', ylab = 'Values')
  dev.off()

  # restore the original working directory
  print("all diagnostic plots are saved under ./mcmcDiagnostics")
  setwd(dir.old)
} 