#' @title Bayesian Hierarchical Regression on Clearance Rates
#'
#' @description
#' \code{clearanceEstimatorBayes} estimates the parasite clearance rates by using a Bayesian hierarchical model. 
#' Moreoever, it provides regression analysis of clearance rates on given covariates.
#' 
#' @param data a data frame containing the profiles of patients.
#' This data frame must contain \code{id}, \code{time}, and \code{count} columns, 
#' in that order. The first column represents the IDs of patients. 
#' The second and third columns contain parasite measurements (per microliter) in different times.
#' @param covariates an optional data frame containing individual level covaraites. This argument may be \code{NULL}, 
#' in which case estimation of clearance rates is of primary interest.
#' @param seed a user-sepcified number used to initialize a pseudorandom number generator. 
#' The default value is set to be 1234 for reproducibility. If \code{seed = NULL}, then its value will be 
#' automatically obtained from the system clock.
#' @param detect.limit detection limit of the parasite density in blood (parasites per microliter)
#' @param outlier.detect indicator of whether or not to use Flegg's outlier detection method.
#' \code{outlier.detect = TRUE} is recommended.
#' @param conf.level required confidence level for reporting credible intervals
#' @param niteration total number of simulations after the burn-in period
#' @param burnin length of the burn-in priod in the MCMC used in \code{clearanceEstimatorBayes}
#' @param thin step size of the thinning process in the MCMC used in \code{clearanceEstimatorBayes}
#' @param filename the name of the csv file used to store some output elements. This file contains
#' \code{id}, \code{clearance.mean}, \code{lag.median}, and \code{tail.median}.
#'  
#' @import stats
#' @import MCMCpack
#' @import MASS
#' @import msm
#' @import mvtnorm
#' @import survival
#' @import AER
#' @import Cairo
#' @import utils
#' 
#' @export
#' 
#' @return The function \link{summary} (i.e., \link{summary.bhrcr}) can be used to obtain a summary of the results.
#' \code{clearanceEstimatorBayes} returns an object of class "bhrcr" which is a list containing:
#' \item{CALL}{function call}
#' \item{clearance.post}{posterior distributions of clearance rates}
#' \item{clearance.mean}{mean values of the posterior distributions of clearance rates}
#' \item{clearance.median}{median values of the posterior distributions of clearance rates}
#' \item{intercept.post}{posterior distributions of the intercepts (alpha_i's) in the model}
#' \item{gamma.post}{posterior distribution of gamma}
#' \item{gamma.post.thin}{thinned posterior sample of gamma}
#' \item{gamma.mean}{mean values of the posterior distribution of gamma}
#' \item{gamma.median}{median values of the posterior distribution of gamma}
#' \item{gamma.CI}{Credible intervals for gamma}
#' \item{halflifeslope.post}{posterior distribution for the effect of covariates on log half-lives}
#' \item{halflifeslope.mean}{mean values of the posterior distribution for the effect of covariates on log half-lives}
#' \item{halflifeslope.median}{median values of the posterior distribution for the effect of covariates on log half-lives}
#' \item{halflifeslope.CI}{Credible intervals for the effect of covariates on log half-lives}
#' \item{predicted.pce}{PCE estimates}
#' \item{eta.post}{posterior distribution of eta}
#' \item{changelag.post}{posterior distributions of changetime between lag and decay phases}
#' \item{changetail.post}{posterior distributions of changetime between decay and tail phases}
#' \item{lag.median}{median values of the posterior distributions of changetime between lag and decay phases}
#' \item{tail.median}{median values of the posterior distributions of changetime between decay and tail phases}
#' \item{var.epsilon.post}{posterior variance of epsilon after simulation}
#' \item{var.error.post}{thinned posterior sample of variance of epsilon}
#' \item{var.alpha.post}{posterior distribution of variance of alpha}
#' \item{var.beta.post}{posterior distribution of variance of beta}
#' \item{index}{a list containing each patient's indices in the data}
#' \item{counts}{Original parasite counts of all patients}
#' \item{counts.current}{Parasite counts of all patients after sampling censored measurements}
#' \item{t.overall}{measurement times of all patients}
#' \item{p.lag}{posterior value of the priori probability of there being a lag phase after simulation}
#' \item{p.lag.thin}{thinned posterior sample of the priori probability of there being a lag phase}
#' \item{p.tail}{posterior value of the priori probability of there being a tail phase after simulation}
#' \item{p.tail.thin}{thinned posterior sample of the priori probability of there being a tail phase}
#' \item{var1.post}{posterior distribution of c^2}
#' \item{var2.post}{posterior distribution of d^2}
#' \item{mu1.post}{posterior distribution of a}
#' \item{mu2.post}{posterior distribution of b}
#' \item{detect.limit}{the detection limit of parasitemia}
#' \item{lag.post}{posterior distributions of index of changetime between lag and decay phases}
#' \item{lag2.post}{posterior distributions of index of changetime between decay and tail phases}
#' \item{theta.post}{posterior distributions of log-parasite-count's mean in lag phase}
#' \item{theta2.post}{posterior distributions of log-parasite-count's mean in tail phase}
#' \item{burnin}{length of the burn-in period}
#' 
#' @details
#' This function estimates parasite clearance rates, along with the effect of covariates on them, 
#' by using the Bayesian hierarchical model which was introduced in Fogarty et al. (2015). 
#' A change point model is used on the log of the parasite densities to account for three potential phases: 
#' (1) a constant phase (the lag phase); (2) a phase with a linear decrease (decay phase); 
#' (3) another constant phase (the tail phase). Hence the estimation of the parasite clearance rate is only 
#' based on observations within the decay phase. The Bayesian approach allows us to treat the delineation between
#' lag, decay, and tail phases within an individual's clearance profile as themselves being random variables,
#' thus taking into account the additional uncertainty of boundaries between phases.
#' Details are in Fogarty et al. (2015).
#' 
#' @author Colin B. Fogarty <cfogarty@mit.edu>, Saeed Sharifi-Malvajerdi <saeedsh@wharton.upenn.edu>, Feiyu Zhu <feiyuzhu@sas.upenn.edu>
#' 
#' @references Flegg, J. A., Guerin, P. J., White, N. J., & Stepniewska, K. (2011). 
#' Standardizing the measurement of parasite clearance in falciparum malaria: the parasite clearance estimator. 
#' Malaria journal, 10(1), 339.
#' @references Fogarty, C. B., Fay, M. P., Flegg, J. A., Stepniewska, K., Fairhurst, R. M., & Small, D. S. (2015). 
#' Bayesian hierarchical regression on clearance rates in the presence of "lag" and "tail" phases 
#' with an application to malaria parasites. Biometrics, 71(3), 751-759.
#' 
#' @examples
#' \dontshow{
#' data("pursat")
#' data("pursat_covariates")
#' data = pursat[pursat["id"] <= 80 & pursat["id"] > 70,]
#' covariates = pursat_covariates[71:80,]
#' out <- clearanceEstimatorBayes(data = data, covariates = covariates, outlier.detect = TRUE,
#'                               niteration = 3, burnin = 1, thin = 1)
#' }
#' \donttest{
#' data("pursat")
#' data("pursat_covariates")
#' out <- clearanceEstimatorBayes(data = pursat, covariates = pursat_covariates, outlier.detect = TRUE,
#'                                niteration = 200, burnin = 50, thin = 10)
#' }
#' 


clearanceEstimatorBayes = function(data, covariates = NULL, seed = 1234, detect.limit = 40, outlier.detect = TRUE, conf.level = .95, niteration = 100000, burnin = 500, thin = 50, filename = "output.csv"){

base.log = exp(1)

if (!identical(
  dimnames(data)[[2]][1:3],
  c("id","time","count"))){
  stop("first three columns of data must be named 'id', 'time', and 'count' in that order, names are case sensitive")
}

if(!is.null(seed)) set.seed(seed)

# show users this message ealier than the actual start of MCMC
print(paste("Progress starts. It may take a while..."))

# save the current working directory
dir.old <- getwd()
# get current working directory
dir <- getwd()

if(is.null(data$Predicted)) {
  par(ask = FALSE)
  #cat("\n Calculating WWARN PCE Estimates... \n")
  print("Calculating WWARN PCE Estimates...")
  if (!dir.exists(paste(dir, "/PceEstimates", sep="")))
    {
      dir.create(paste(dir, "/PceEstimates", sep=""))
    }
    setwd(paste(dir, "/PceEstimates", sep=""))
    sink("PceEstimates.txt")
    suppressWarnings(
        suppressMessages(calculateParasiteClearance(Name = "pceestimates", minpara = detect.limit, data1 = data))
    )
    sink()
	data = read.csv("pceestimates_Estimates.csv")
	data = as.data.frame(data)
	
	# If users want to use Flegg's method to detect outliers (outlier.detect == TRUE), then
	# after running the PCE method, we detect and exclude points with outlier status 1 and 2
	# outlier status 3 appears in a tail phase, we don't remove that
	if (outlier.detect == TRUE) {
	  data = data[(data$Outlier!=1 & data$Outlier!=2), ]
	}
	
	# If users don't want to detect outliers in this case (i.e. data doesn't have the Predicted column)
	# a warning message should be reported
	if (outlier.detect == FALSE) {
	  warning("You may forget to detect outliers.", 
	          "\n If outliers before the tail phase were removed from the input data by yourself, please ignore this warning.",
	          "\n Otherwise, please re-run the function by setting outlier.detect = TRUE and the outliers will be automatically removed for you.")
	}
}
# restore the original working directory
setwd(dir.old)

# If data contains predicted PCE (generated by users or by our program) but the type 1 or type 2 outliers are not eliminated, 
# then the program ends and produces an error message.
if(!is.null(data$Precited) & (is.element(1, data$Outlier) | is.element(2, data$Outlier))) {
  stop("There are still outliers before the tail phase, please remove them and then predict Flegg's PCE.",
       "\nYou can choose to delete the old predicted PCE column and then run the program with outlier.detect = TRUE,",
       "\nwhich will automatically clean the data set and generate a predicated PCE column by Flegg's method.")
}

#cat("\n Conducting Bayesian Analysis... \n")
print("Conducting Bayesian Analysis...")
id = data$id
unique.id = unique(id)
t.overall = data[,2]
if(is.null(t.overall))
{
  t.overall = data$time
}
counts = data[,3]
predicted.flegg = data$Predicted/log(base.log)
outlier.detect = F

if(!is.null(covariates)){
  d = c()
  covariates = as.data.frame(covariates)
  for(i in 1:ncol(covariates))
  {
    if(all(covariates[,i] == 1))
    {
      d = c(d, i)
    }
  }
  if(!is.null(d)){
    covariates = covariates[,-c(d)]
    
  }
  covariates = as.data.frame(covariates)
  names.cov = c(colnames(covariates))
  b = runif(nrow(covariates))
  X.cov = lm(b ~., data = as.data.frame(covariates), x = T)$x
  X.cov = as.matrix(X.cov)
  
} 


if(is.null(covariates)){ 
  X.cov = t(t(rep(1, length(unique.id))))
}

nobs = length(unique(id))

lag2.current = rep(0, nobs)

nsim = burnin + niteration

index = vector("list", nobs)
alpha.current <- rep(1, nobs)
beta.current <- rep(-.05, nobs)
outlier = vector("list",nobs)
maximum.t = rep(0, nobs)
for(i in 1:nobs)
{
  
  temp = counts[id==unique.id[i]]
  count.temp = temp
  count.temp[count.temp == 0] = detect.limit
  index[[i]] = which(id == unique.id[i]) #indices for the ith individual
  maximum.t[i] = max(t.overall[index[[i]]]) #maximum t observed for ith individual
  outlier[[i]] = rep(0,length(index[[i]]))
  lag2.current[i] = length(count.temp)
  vec = lm(log(count.temp, base = base.log)~t.overall[index[[i]]])$coef
  alpha.current[i] <- vec[1]
  beta.current[i] <- vec[2]
  if(beta.current[i] >= 0 || is.na(beta.current[i]))
  {
    beta.current[i] = -.05/log(base.log, 10)
  }
  if(alpha.current[i] <= 0 || is.na(alpha.current[i]))
  {
    alpha.current[i] = 5/log(base.log, 10)
  }
  
  
}
prior.lag = 1 
prior.tail = 1 
prior.nolag = 1 
prior.notail = 1
change.current = rep(0, nobs)
change2.current = maximum.t
lag.current = rep(0, nobs)
outlier.post = rep(list(outlier),nsim)
theta.current = alpha.current
theta2.current = rep(0, nobs)
gamma.current = rep(0, ncol(X.cov))
eta.current = rep(0, ncol(X.cov))
alpha.pop.current = rep(0, nobs)
beta.pop.current = rep(0, nobs)
var.epsilon.current = 1/(log(base.log, 10))^2
p.lag = .25
p.tail = .25
c.out = 2

var.epsilonout.current = 5/(log(base.log, 10))^2
var.theta.current = 1
var.theta2.current = 1
var.beta.current = 1
var.alpha.current = 1
mu1.current =  log(3) 
mu2.current = log(1)
var1.current = 1
var2.current = 1
mu1.post = rep(0, nsim)
mu2.post = rep(0, nsim)
var1.post = rep(0, nsim)
var2.post = rep(0, nsim)
probs = rep(0, length(t.overall))
var.epsilon.post = rep(1, nsim)
var.epsilonout.post = rep(1, nsim)
var.theta.post = rep(1, nsim)
var.theta2.post = rep(1, nsim)
alpha.post = matrix(0, nobs, nsim)
beta.post = matrix(0, nobs, nsim)
beta.pop.post = matrix(0, nobs, nsim)
alpha.pop.post = matrix(0, nobs, nsim)
gamma.post = matrix(0, ncol(X.cov), nsim)
halflifeslope.post = matrix(0, ncol(X.cov), nsim)
eta.post = matrix(0, ncol(X.cov), nsim)
theta.post = matrix(0, nobs, nsim)
theta2.post = matrix(0, nobs, nsim)
var.alpha.post = rep(1, nsim)
c.out.post = rep(2, nsim)
var.beta.post = rep(1, nsim)
p.lag.post = rep(0, nsim)
p.tail.post = rep(0, nsim)
shape.post = rep(0,nsim)
rate.post = rep(0, nsim)
lag.post = matrix(0, nobs, nsim)
lag2.post = matrix(0, nobs, nsim)
change.post = matrix(0, nobs, nsim)
change2.post = matrix(0, nobs, nsim)
counts.current = counts
X = cbind(rep(1, length(t.overall)), t.overall)
prob.out = .2
same = rep(0, nobs)
prob.out.post = rep(0,nsim)
tripX = solve(t(X.cov)%*%X.cov)%*%t(X.cov)

for(i in 1:nsim)
{
  SSE = 0
  sum.delta = 0
  sum.delta2 = 0
  n.total = 0
  n.delta = 0
  n.delta2 = 0
  SSEout = 0
  sum.deltaout = 0
  sum.delta2out = 0
  s.delta2out= rep(0, nobs)
  n.totalout = 0
  n.deltaout = 0
  n.delta2out = 0
  nindivid.out = rep(0, nobs)
  for(j in 1:nobs)
  {	
    ind = index[[j]]
    out = outlier[[j]]
    cc = counts[ind]
    t = t.overall[ind]
    lt = length(t)
    #################### sample Censored Observations ######################
    s = which(cc < detect.limit)
    if(length(s) > 0)
    {
      for(k in s)
      {
        if(lag2.current[j] < k)
        {
          new = base.log^(rtnorm(1, theta2.current[j], sqrt(var.epsilon.current)*(1-out[k]) + sqrt(var.epsilonout.current)*out[k], upper = log(detect.limit,base=base.log)))
        }
        else
        {
          new = base.log^(rtnorm(1, alpha.current[j] + beta.current[j]*t[k], sqrt(var.epsilon.current)*(1-out[k]) + sqrt(var.epsilonout.current)*out[k], upper = log(detect.limit,base = base.log)))
        }
        counts.current[ind[k]] = new
      }
    }
    #########################################################################
    
    lc = unlist(log(counts.current[index[[j]]], base = base.log))
    beta.old = -beta.current[j]
    alpha.old = alpha.current[j]
    
    t.ind = ((lag.current[j]+1):lag2.current[j])
    
    ########## Defining the z_ij values ############################
    if(lag.current[j] != lag2.current[j])
    {
      tnew = c(rep(change.current[j], lag.current[j]), t[t.ind], rep(change2.current[j], (length(t)-lag2.current[j])))
    }else{
      tnew = c(rep(change.current[j], lag.current[j]), rep(change2.current[j], (length(t)-lag2.current[j])))
    }
    ################################################################
    outnew = out
    Xt = cbind(rep(1, length(t)), -tnew)
    
    m.a = exp(alpha.pop.current[j] + .5*var.alpha.current)
    m.b = exp(beta.pop.current[j] + .5*var.beta.current)
    v.a = (exp(var.alpha.current) - 1)*exp(2*alpha.pop.current[j] + var.alpha.current)
    v.b = (exp(var.beta.current) - 1)*exp(2*beta.pop.current[j] + var.beta.current)
    y.new = unlist(c(lc, m.a, m.b))
    X.new = rbind(Xt, diag(2))
    V.new = var.epsilon.current*diag((lt+2))
    if(lt != 1)
    {
      diag(V.new[1:lt,1:lt]) = var.epsilon.current*(1-outnew) + var.epsilonout.current*outnew
    }else {
      V.new[1:lt,1:lt] = var.epsilon.current*(1-outnew) + var.epsilonout.current*outnew
    }	   
    Vinv = diag(1/diag(V.new[1:lt, 1:lt])) 
    if(lt == 1)
    {
      Vinv = 1/V.new[1,1]
    }
    diag(V.new)[(lt+1):(lt+2)] = c(v.a, v.b)
    X3 = solve(t(X.new)%*%solve(V.new)%*%X.new)%*%t(X.new)%*%solve(V.new)
    Var.b = 1.3*solve(t(X.new)%*%solve(V.new)%*%X.new)
    mean = X3%*%y.new
    t.list = -t.overall[index[[j]]]
    vec = c(0,0)
    
    vec.old = c(alpha.old, beta.old)
    vec[2] = rtnorm(1, mean[2], sqrt(Var.b[2,2]), lower = 0)
    condvar = Var.b[1,1] - Var.b[2,1]^2/Var.b[2,2]
    condmean = mean[1] + Var.b[2,1]/Var.b[2,2]*(vec[2] - mean[2])
    vec[1] = rtnorm(1, condmean, sqrt(condvar), lower = 0)
    condmean.old = mean[1] + Var.b[2,1]/Var.b[2,2]*(vec.old[2] - mean[2])
    beta.new = vec[2]
    alpha.new = vec[1]
    lratio = -.5*t(lc - Xt%*%vec)%*%Vinv%*%(lc-Xt%*%vec)-log(beta.new) - log(alpha.new) - (log(beta.new)-beta.pop.current[j])^2/(2*var.beta.current) -(log(alpha.new)-alpha.pop.current[j])^2/(2*var.alpha.current) - (-.5*t(lc - Xt%*%vec.old)%*%Vinv%*%(lc-Xt%*%vec.old) -log(beta.old) - log(alpha.old) - (log(beta.old)-beta.pop.current[j])^2/(2*var.beta.current) - (log(alpha.old)-alpha.pop.current[j])^2/(2*var.alpha.current)) + dtnorm(beta.old, mean[2], sqrt(Var.b[2,2]), lower = 0, log = T) + dtnorm(alpha.old, condmean.old, sqrt(condvar), lower = 0,log = T) - dtnorm(beta.new, mean[2], sqrt(Var.b[2,2]), lower = 0, log = T) - dtnorm(alpha.new, condmean, sqrt(condvar), lower = 0, log = T)
    # lratio 
    vec[2] = -beta.old
    vec[1] = alpha.old
    if(log(runif(1)) <= lratio)
    {
      beta.current[j] = -beta.new
      alpha.current[j] = alpha.new
      theta2.current[j] = alpha.current[j]+ beta.current[j]*change2.current[j]
      theta.current[j] = alpha.current[j] + beta.current[j]*change.current[j]
    }
    
    if(i > 2)
    {
      if(beta.current[j] == beta.post[j,(i-1)])
      {
        same[j] = same[j] + 1
      }else{same[j] = 0}
    }				
    
    lc = unlist(log(counts.current[index[[j]]], base = base.log))
    t = t.overall[index[[j]]]
    int = rep(0, lag2.current[j]+1)
    for(k in 0:(lag2.current[j]))
    {
      
      if(k == 0)
      {
        int[(k+1)] = densityintegrate.LN(0, change2.current[j], 0, lag2.current[j], alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu1.current, var1.current)*(1-p.lag)
      }else if(k == length(t)){
        int[(k+1)] = densityintegrate.LN(t[length(t)], change2.current[j], k, lag2.current[j], alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu1.current, var1.current)*(p.lag)
      }else{
        int[(k+1)] = integrate(densityintegrate.LN, lower = t[k], upper = min(t[(k+1)], change2.current[j]), change2.current[j], k, lag2.current[j], alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu1.current, var1.current)$value*(p.lag)
      }	
    }
    
    
    probs = int/sum(int)
    lag.current[j] = sample(0:(lag2.current[j]),size = 1, prob = probs)
    choice = lag.current[j] + 1
    
    k = lag.current[j]
    if(k == 0)
    {
      change.current[j] = 0
    }else if(k==length(index[[j]])){change.current[j] = t[length(t)]
    
    }else{	
      logmaxval = optimize(densitymax.LN, c(t[k], min(t[(k+1)], change2.current[j])), change2.current[j], k, lag2.current[j], alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu1.current, var1.current, maximum = TRUE)$objective
      stop = FALSE
      while(stop == FALSE)
      {
        prop = runif(1, t[k], min(t[(k+1)], change2.current[j]))
        logdensity.value = densitymax.LN(prop, change2.current[j], k, lag2.current[j], alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out,mu1.current, var1.current)
        ratio = exp(logdensity.value- logmaxval)
        if(runif(1) <= ratio){
          change.current[j] = prop
          stop = TRUE
        }
      }
    }
    
    
    
    
    N = length((lag.current[j]):length(index[[j]]))
    int = rep(0, N)
    for(k in (lag.current[j]):(length(index[[j]])))
    {
      
      if(k == length(index[[j]]))
      {
        int[(k - (lag.current[j]-1))] = densityintegrate2.LN(maximum.t[j], change.current[j], lag.current[j], k, alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu2.current, var2.current)*(1-p.tail)
      }else if(k == 0){
        int[(k - (lag.current[j]-1))] = densityintegrate2.LN(0, change.current[j], lag.current[j], k, alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu2.current, var2.current)*(p.tail)
      }else{
        
        int[(k - (lag.current[j]-1))] = integrate(densityintegrate2.LN, lower = max(t[k], change.current[j]), upper = t[(k+1)], change.current[j], lag.current[j], k, alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu2.current, var2.current)$value*(p.tail)
      }
      
    }
    prior = c(rep((1-p.tail), (N-1)), p.tail)
    probs2 = int/sum(int)
    if(length(probs2) > 1)
    {
      lag2.current[j] = sample((lag.current[j]):(length(index[[j]])), size = 1, prob = probs2)
    }
    if(length(probs2)==0)
    {
      lag2.current[j] = length(index[[j]])
    }
    
    k = lag2.current[j]
    
    if(k == length(index[[j]]))
    {
      change2.current[j] = t[length(t)]
    }else if(k==0){
      change2.current[j] = 0
    }else{
      logmaxval2 = optimize(densitymax2.LN, c(max(t[k], change.current[j]), t[k+1]), change.current[j], lag.current[j], k, alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu2.current, var2.current, maximum = TRUE)$objective
      
      
      stop = FALSE
      while(stop == FALSE)
      {
        prop = runif(1, max(t[k], change.current[j]), t[(k+1)])
        logdensity.value = densitymax2.LN(prop, change.current[j], lag.current[j], k, alpha.current[j], beta.current[j], t, lc, var.epsilon.current, var.epsilonout.current, out, mu2.current, var2.current)
        ratio = exp(logdensity.value- logmaxval2)
        if(runif(1) <= ratio)
        {
          change2.current[j] = prop
          stop = TRUE
        }
      }
    }
    
    
    
    
    if(outlier.detect == T)
    {
      outlier[[j]] = identify.outlier(alpha.current[j], beta.current[j],theta.current[j], theta2.current[j], var.epsilon.current, var.epsilonout.current, var.epsilon.current, var.epsilonout.current, var.epsilon.current, var.epsilonout.current,t, lc, lag.current[j], lag2.current[j], prob.out)
    }
    out.current = outlier[[j]]
    
    ind = index[[j]]
    t.ind = (lag.current[j]+1):(lag2.current[j])	
    if(lag.current[j] != lag2.current[j])
    {
      tnew = c(rep(change.current[j], lag.current[j]), t[t.ind], rep(change2.current[j], (length(t)-lag2.current[j])))
    }else{
      tnew = c(rep(change.current[j], lag.current[j]), rep(change2.current[j], (length(t)-lag2.current[j])))
    }
    
    
    fit = alpha.current[j] + beta.current[j]*tnew
    
    SSE = SSE + sum((1-out.current)*(lc - fit)^2)
    SSEout = SSEout + sum((out.current)*(lc - fit)^2)
    n.total = n.total + sum(1-out.current)
    n.totalout = n.totalout + sum(out.current)
  }
  
  shape = .5*(n.total + n.delta + n.delta2 +n.totalout + n.deltaout + n.delta2out) 
  rate = .5*(SSE + sum.delta + sum.delta2) + .5*(SSEout + sum.deltaout + sum.delta2out)/c.out
  var.epsilon.current = 1/rgamma(1, shape, rate)
  
  
  
  
  if(outlier.detect == T)
  {
    var.epsilonout.current = var.epsilon.current
    shapeout = .5*(n.totalout + n.deltaout + n.delta2out) + 2.5
    rateout = .5*(SSEout + sum.deltaout + sum.delta2out)/var.epsilon.current + 5
    e = 1
    while(var.epsilonout.current <= var.epsilon.current & e <= 40)
    {
      c.out = 1/rgamma(1, shapeout, rateout)
      var.epsilonout.current = c.out*var.epsilon.current
      if(e == 40)
      {
        var.epsilonout.current = var.epsilon.current
      }
      e = e + 1
    }
    if(var.epsilonout.current > 1e3)
    {
      var.epsilonout.current = 1e3
    }
  }	
  gamma.current = mvrnorm(1, tripX%*%log(-beta.current), var.beta.current*solve(t(X.cov)%*%X.cov))
  eta.current = mvrnorm(1, tripX%*%log(alpha.current), var.alpha.current*solve(t(X.cov)%*%X.cov))
  theta.pop.current = rnorm(1,mean(theta.current), sqrt(var.theta.current/nobs))
  
  beta.pop.current = X.cov%*%gamma.current
  alpha.pop.current = X.cov%*%eta.current
  shape.beta = .5*nobs - 1/2
  rate.beta = .5*sum((log(-beta.current) - beta.pop.current)^2)
  var.beta.current = 1/rgamma(1, shape.beta, rate.beta)
  shape.alpha = .5*nobs - 1/2
  rate.alpha = .5*sum((log(alpha.current) - alpha.pop.current)^2)
  var.alpha.current = 1/rgamma(1, shape.alpha, rate.alpha)
  
  
  lcurrent = log(counts.current, base = base.log)
  
  c1.temp = change.current[change.current > 0]
  c2.temp = change2.current[change.current > 0]
  if(length(c1.temp) >= 0)
  {
    
    v.post = 1/(length(c1.temp)/var1.current + 1/.25)
    m.post = v.post*(sum(log(c1.temp))/var1.current + log(6)/.25) 		
    
    
    mu1.prop = rnorm(1, m.post, sqrt(v.post))
    
    lratio = sum(pnorm(log(c2.temp)), mu1.current, sqrt(var1.current), log.p = T) -  sum(pnorm(log(c2.temp)), mu1.prop, sqrt(var1.current), log.p = T)
    if(log(runif(1)) < lratio)
    {
      mu1.current = mu1.prop
    }
    var1.prop = 1/rgamma(1, length(c1.temp)/2+1, .5*sum((log(c1.temp) - mu1.current)^2) + 1 )
    lratio = sum(pnorm(log(c2.temp)), mu1.current, sqrt(var1.current), log.p = T) -  sum(pnorm(log(c2.temp)), mu1.current, sqrt(var1.prop), log.p = T)
    if(log(runif(1)) < lratio)
    {
      var1.current= var1.prop
    }
    
  }
  c2.temp = change2.current[change2.current < maximum.t]
  c1.temp =  change.current[change2.current < maximum.t]
  maxt.temp = maximum.t[change2.current < maximum.t]
  if(length(c2.temp) >= 0)
  {
    v.post = 1/(length(c2.temp)/var2.current + 1/.25)
    m.post = v.post*(sum(log(maxt.temp - c2.temp))/var2.current + log(6)/.25)
    
    mu2.prop = rnorm(1, m.post, sqrt(v.post))
    lratio = sum(pnorm(log(maxt.temp - c1.temp)), mu2.current, sqrt(var2.current), log.p = T) -  sum(pnorm(log(maxt.temp - c1.temp)), mu2.prop, sqrt(var2.current), log.p = T)
    if(log(runif(1)) < lratio)
    {
      mu2.current = mu2.prop
    }
    var2.prop = 1/rgamma(1, length(c1.temp)/2 + 1, .5*sum((log(maxt.temp - c2.temp) - mu2.current)^2) + 1 )
    lratio = sum(pnorm(log(maxt.temp - c1.temp)), mu2.current, sqrt(var2.current), log.p = T) -  sum(pnorm(log(maxt.temp - c1.temp)), mu2.current, sqrt(var2.prop), log.p = T)
    if(log(runif(1)) < lratio)
    {
      var2.current= var2.prop
    }
    
  }
  
  
  
  if(outlier.detect == T)
  {	
    prob.out = rbeta(1,length(t.overall)*.005 +n.deltaout + n.delta2out + n.totalout, length(t.overall)*.095 + n.delta + n.delta2 + n.total)
  }
  
  nlag = sum(lag.current!=0)
  ntail = sum(change2.current != maximum.t)
  p.lag = rbeta(1, nlag+ prior.lag, nobs - nlag + prior.nolag)
  p.tail = rbeta(1, ntail+prior.tail, nobs - ntail + prior.notail)
  beta.post[,i] = beta.current
  alpha.post[,i] = alpha.current
  theta.post[,i] = theta.current
  theta2.post[,i] = theta2.current
  beta.pop.post[,i] = beta.pop.current
  alpha.pop.post[,i] = alpha.pop.current
  gamma.post[,i] = gamma.current
  eta.post[,i] = eta.current	
  var.epsilon.post[i] = var.epsilon.current
  var.epsilonout.post[i] = var.epsilonout.current
  c.out.post[i] = c.out
  p.lag.post[i] = p.lag
  p.tail.post[i] = p.tail
  var.beta.post[i] = var.beta.current
  var.alpha.post[i] = var.alpha.current
  var.theta.post[i] = var.theta.current
  var.theta2.post[i] = var.theta2.current
  lag.post[,i] = lag.current
  lag2.post[,i] = lag2.current
  change.post[,i] = change.current
  change2.post[,i] = change2.current
  prob.out.post[i] = prob.out
  outlier.post[[i]] = outlier
  var1.post[i] = var1.current
  var2.post[i] = var2.current
  mu1.post[i] = mu1.current
  mu2.post[i] = mu2.current
  
  #if(i == 1)
      # print(paste("Progress starts. It may take a while..."))
  
  if(i %% 50 == 0)
  {
    print(paste("Progress:", i, "out of", nsim))
  }
  
}


thin = seq.int(burnin+1, nsim, thin)

halflifeslope.post = as.vector( solve(t(X.cov)%*%X.cov)%*%t(X.cov)%*%rep(log(log(2)), nrow(X.cov)))- gamma.post

clearance.median = -apply(beta.post[,thin], 1, median)
clearance.mean = -apply(beta.post[,thin], 1, mean)
gamma.median = apply(gamma.post[,thin], 1, median)
gamma.mean = apply(gamma.post[,thin], 1, mean)
halflifeslope.median = apply(halflifeslope.post[,thin], 1, median)
halflifeslope.mean = apply(halflifeslope.post[,thin], 1, mean)

outlier.freq = outlier
for(i in 1:nobs)
{
  temp = matrix(0, length(thin),length(index[[i]]))
  for(j in 1:length(thin))
  {
    temp[j,] = outlier.post[[thin[j]]][[i]]
  }
  outlier.freq[[i]] = colMeans(temp)
}
lag.median = apply(change.post[,thin], 1, median)
tail.median = apply(change2.post[,thin], 1, median)

a = (1-conf.level)/2
A = 100*a

if(!is.null(covariates))
{
  gamma.CI = apply(gamma.post[,thin], 1, quantile, c(a, .5, 1-a))
  gamma.CI = t(gamma.CI)
  rownames(gamma.CI) = colnames(X.cov)
  colnames(gamma.CI) = c(paste("CI ", A, "%", sep = ""), "Median", paste("CI ", 100-A, "%", sep = ""))
  
  halflifeslope.CI = apply(halflifeslope.post[,thin], 1, quantile, c(a, .5, 1-a))
  halflifeslope.CI = t(halflifeslope.CI)
  rownames(halflifeslope.CI) = colnames(X.cov)
  colnames(halflifeslope.CI) = c(paste("CI ", A, "%", sep = ""), "Median", paste("CI ", 100-A, "%", sep = ""))
}else{
  
  
  gamma.CI = quantile(gamma.post, c(a, .5, 1-a))
  names(gamma.CI) = c(paste("CI ", A, "%", sep = ""), "Median", paste("CI ", 100-A, "%", sep = ""))
  
  halflifeslope.CI = quantile(halflifeslope.post, c(a, .5, 1-a))
  names(halflifeslope.CI) = c(paste("CI ", A, "%", sep = ""), "Median", paste("CI ", 100-A, "%", sep = ""))
  
}

write.csv(cbind(unique.id, clearance.mean, lag.median, tail.median), filename)

rownames(change.post)  <- unique.id
rownames(change2.post) <- unique.id
names(clearance.mean) <- unique.id
names(clearance.median) <- unique.id
clearance.post = -beta.post[,thin]
dimnames(clearance.post) <- list(unique.id, NULL)
dimnames(alpha.post) <- list(unique.id, NULL)

rownames(lag.post) <- unique.id
rownames(lag2.post) <- unique.id
rownames(theta.post) <- unique.id
rownames(theta2.post) <- unique.id
names(lag.median) <- unique.id
names(tail.median) <- unique.id

## We define a class called "bhrcr" for the output, following the rule of S3 method in R.
results <- list(CALL = match.call(), clearance.median = clearance.median, clearance.mean = clearance.mean, clearance.post = clearance.post, 
                intercept.post = alpha.post[,thin], gamma.post = gamma.post, gamma.post.thin = gamma.post[,thin], gamma.median = gamma.median, gamma.mean = gamma.mean, 
                gamma.CI = gamma.CI, halflifeslope.post = halflifeslope.post[,thin], halflifeslope.median = halflifeslope.median,
                halflifeslope.mean = halflifeslope.mean, halflifeslope.CI = halflifeslope.CI,
                eta.post = eta.post[,thin], var.epsilon.post = var.epsilon.post, var.error.post = var.epsilon.post[thin], changelag.post = change.post[,thin], 
                changetail.post = change2.post[,thin], index = index, lag.post = lag.post[,thin], 
                lag2.post = lag2.post[,thin], counts.current = counts.current, t.overall = t.overall, 
                theta.post = theta.post[,thin], theta2.post = theta2.post[,thin], predicted.pce = predicted.flegg, 
                counts = counts, var.beta.post = var.beta.post[thin], var.alpha.post = var.alpha.post[thin], 
                p.lag = p.lag.post, p.lag.thin = p.lag.post[thin], p.tail = p.tail.post, p.tail.thin = p.tail.post[thin], detect.limit = detect.limit,
                var1.post = var1.post, var2.post = var2.post, mu1.post = mu1.post, mu2.post = mu2.post, 
                lag.median = lag.median, tail.median = tail.median, burnin = burnin, id = unique.id)
class(results) <- c("bhrcr", "list")
return(results)
}