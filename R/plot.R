#' @title Bayesian Clearance Estimator Plotting
#'
#' @description
#' \code{plot.bhrcr} plots the posterior results from \code{clearanceEstimatorBayes}.
#' 
#' @param x output given by \code{clearanceEstimatorBayes}
#' @param plot.post indicator of whether or not the posterior samples should be plotted
#' @param id.plot patients' IDs
#' @param thin an optional vector showing which posterior samples to be plotted
#' @param ... additional arguments passed to the \code{plot.bhrcr} function
#'
#' @return the directory location under which all the plots are saved.
#'
#' @details
#' This function plots clearance profile of each individual along with their fitted Bayesian model, 
#' Flegg's PCE estimates, and posterior samples.
#'
#' @import graphics
#' @import grDevices
#' 
#' @method plot bhrcr
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
#' plot(out)
#' }
#' \donttest{
#' data("posterior")
#' plot(posterior)
#' }
#' \donttest{
#' data("pursat")
#' data("pursat_covariates")
#' out <- clearanceEstimatorBayes(data = pursat, covariates = pursat_covariates, 
#'                                niteration = 200, burnin = 50, thin = 10)
#' plot(out)
#' }


plot.bhrcr <- function(x, plot.post = T, id.plot = NULL, thin = NULL, ...) {

# check input class
if (!inherits(x, "bhrcr")) stop("Object must be of class 'bhrcr'")

# save the current working directory
dir.old <- getwd()
# get current working directory
dir <- getwd()   
# create a sub-directory for storing diagnostic plots
if (!dir.exists(paste(dir, "/plots", sep=""))){
    dir.create(paste(dir, "/plots", sep=""))
}
setwd(paste(dir,"/plots", sep=""))

log.base = x$log.base
detect.limit = x$detect.limit
begin = 1

nobs=length(x$clearance.median)
end = length(x$var.error.post)
if(is.null(thin))
{
thin = seq.int(1,length(x$var.error.post), by = 1)
}
index = x$index

if(is.null(id.plot))
{
	id.plot = 1:nobs
}
for(j in id.plot) {
    filename = paste('patient_', j, '.pdf', sep = '') 
    pdf(file = filename)
    
	lc = log(x$counts[index[[j]]], base = log.base)
	m1 = names(sort(-table(x$lag.post[j,thin])))[1]
	m2 = names(sort(-table(x$lag2.post[j,thin])))[1]
	xt = x$t.overall[index[[j]]]
	if(length(xt) > 1)
	{
		plot(x$t.overall[index[[j]]], log(x$counts.current[index[[j]]], base = log.base), main = paste("Patient id:",j), type = "n", ylab = "Log Parasite Count", xlab = "Time (Hours)")
	} else { 
		ylim = c(min(x$intercept.post[j,]) - max(x$clearance.post[j,])*(xt[1] + 12), max(x$intercept.post[j,]))
		
		plot(c(xt,xt+12),c(0,12), main = paste("Patient id:",j), type = "n", ylab = "Log Parasite Count", xlab = "Time (Hours)", ylim = ylim)
		
	}
	# Posterior Sample
	posterior.sample = matrix(0,nrow = min(20,length(thin)), ncol = length(x$t.overall[index[[j]]]))
	if(plot.post == T)
	{
	i=1
	for(k in (sample(thin, min(20,length(thin)), replace = F)))
	{
		
	changelag = x$changelag.post[j,k]
	changetail = x$changetail.post[j,k]
	t = x$t.overall[index[[j]]]
	tlag.plot = c(t[t < changelag], changelag)
	tdecay.plot = c(changelag, t[t > changelag & t < changetail], changetail)
	ttail.plot = c(changetail, t[t>changetail])
	tlag = rep(changelag, length(tlag.plot))
	tdecay = tdecay.plot
	if(length(t) == 1)
	{
		tdecay = c(t, t + 12)
	}
	tdecay.plot = tdecay
	ttail = rep(changetail, length(ttail.plot))
	llag = x$intercept.post[j,k] -x$clearance.post[j,k]*tlag
	ldecay = x$intercept.post[j,k] -x$clearance.post[j,k]*tdecay
	ltail = x$intercept.post[j,k] -x$clearance.post[j,k]*ttail
	lines(tlag.plot, llag, col = "grey")
	lines(tdecay.plot, ldecay, col = "grey")
	lines(ttail.plot, ltail, col = "grey")
	t2 = c(tlag.plot, tdecay.plot, ttail.plot)
	llag.new = llag
	ltail.new = ltail
	if (length(tlag.plot) > 1) llag.new = llag[-length(llag)]
	if (length(ttail.plot) > 1) ltail.new = ltail[-1]
	ldecay.new = ldecay[-c(1,length(ldecay))]
	posterior.sample[i,] = c(llag.new,ldecay.new,ltail.new)
	i = i+1
	}
	}
	#######################New Median Curve#############################
	times = x$t.overall[index[[j]]]
	med.val = rep(0,length(times))
	for (t in 1:length(times)){
	  med.val[t] = median(posterior.sample[,t])
	}
	lines(times, med.val,col = "red", lwd = 2)
	####################################################################
	# Mean Curve
	changelag = mean(x$changelag.post[j,])
	changetail = mean(x$changetail.post[j,])
	t = x$t.overall[index[[j]]]
	tlag.plot = c(t[t < changelag], changelag)
	tdecay.plot = c(changelag, t[t > changelag & t < changetail], changetail)
	ttail.plot = c(changetail, t[t>changetail])
	tlag = rep(changelag, length(tlag.plot))
	tdecay = tdecay.plot
	if(length(t) == 1)
	{
		tdecay = c(t, t + 12)
	}
	tdecay.plot = tdecay
	ttail = rep(changetail, length(ttail.plot))
	llag = mean(x$intercept.post[j,begin:end]) + mean(-x$clearance.post[j,begin:end])*tlag
	ldecay = mean(x$intercept.post[j,begin:end]) + mean(-x$clearance.post[j,begin:end])*tdecay
	ltail = mean(x$intercept.post[j,begin:end]) + mean(-x$clearance.post[j,begin:end])*ttail
	lines(tlag.plot, llag, lwd = 2)
	lines(tdecay.plot, ldecay, lwd = 2)
	lines(ttail.plot, ltail, lwd = 2)
	
	# Median Curve
	changelag = median(x$changelag.post[j,])
	changetail = median(x$changetail.post[j,])
	t = x$t.overall[index[[j]]]
	tlag.plot = c(t[t < changelag], changelag)
	tdecay.plot = c(changelag, t[t > changelag & t < changetail], changetail)
	ttail.plot = c(changetail, t[t>changetail])
	tlag = rep(changelag, length(tlag.plot))
	tdecay = tdecay.plot
	if(length(t) == 1)
	{
		tdecay = c(t, t + 12)
	}
	tdecay.plot = tdecay
	ttail = rep(changetail, length(ttail.plot))
	llag = median(x$intercept.post[j,begin:end]) + median(-x$clearance.post[j,begin:end])*tlag
	ldecay = median(x$intercept.post[j,begin:end]) + median(-x$clearance.post[j,begin:end])*tdecay
	ltail = median(x$intercept.post[j,begin:end]) + median(-x$clearance.post[j,begin:end])*ttail
	lines(tlag.plot, llag, lwd = 2, col = "blue")
	lines(tdecay.plot, ldecay, lwd = 2, col = "blue")
	lines(ttail.plot, ltail, lwd = 2, col = "blue")
	
	#hist(beta.post[j, thin])
	# Flegg Estimates
	lines(xt, x$predicted.flegg[index[[j]]], col = "purple", lwd = 2)
	
	#abline(median(x$intercept.post[j,begin:end]), median(-x$clearance.post[j,begin:end]), lwd = 2)
	#abline(mean(x$intercept.post[j,begin:end]), mean(-x$clearance.post[j,begin:end]), col = "blue", lwd = 2)
	# Observations
	points(x$t.overall[index[[j]]][lc>= log(detect.limit, base = log.base)], lc[lc >= log(detect.limit, base = log.base)], pch = 16)
	# Censored Observations
	points(x$t.overall[index[[j]]][lc < log(detect.limit, base = log.base)], rep(log(detect.limit, base = log.base), times = sum(lc < log(detect.limit, base = log.base))), col = "green", pch = 17)
	
	legend("bottomleft", c("Flegg et al Estimate", "Bayes Estimate - Median", "Posterior Median", "Bayes Estimate - Mean", "Posterior Sample", "Censored Observation"), 
	       lwd = c(2,2,2,2,2, NA), pch = c(NA, NA, NA, NA, NA, 17), col=c("purple", "blue","red", "black", "grey", "green"), cex = .8)
	#print(outlier.freq[[j]])
	
	dev.off()
}
# restore the original working directory
setwd(dir.old)
print("all plots are saved under ./plots")
}