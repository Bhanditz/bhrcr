make.probs = function(alpha, beta, theta, var.epsilon, var.epsilonout, var.delta, var.deltaout, t, lag2, change, change2, y, out, rho1, prior = NULL)
{
	lt = length(change)
	if(is.null(prior))
	{
		prior = rep(1,lt)
	}
	prior = prior/sum(prior)
	p = rep(0, lt)
	for(i in 1:(lt))
	{
		tnew = c(rep(change[i], (i-1)), t[(i):lag2], rep(change2, (length(t)-lag2)))
		p[i] = prod(dnorm(y, alpha[i] + beta[i]*tnew, sqrt(var.epsilon)*(1-out)+sqrt(var.epsilonout)*out))*prior[i]*(change[i] + .05)^(-rho1)
	 
	}

	p = p/sum(p)
}
make.probs2 = function(alpha, beta, theta2, var.epsilon, var.epsilonout, var.delta2, var.delta2out, t, lag, change, change2, y, out, rho2, prior = NULL)
{
	ind = (lag+1):length(t)
	tnew = t[ind]

	p = rep(0,(length(ind)-1))
	lt = length(t)
	if(is.null(prior))
	{
		prior = rep(1,lt)
	}
	prior = prior/sum(prior)
	for(i in (lag+2):lt) 
	{
		jj = i - (lag + 1)
		tdesign = c(rep(change, lag), t[(lag+1):i], rep(change2[jj], length(t)-i))
		p[jj] = prod(dnorm(y, alpha[jj] + beta[jj]*tdesign, sqrt(var.epsilon)*(1-out)+sqrt(var.epsilonout)*out))*prior[jj]
	}
		#if(all(p == 0))
		#{
		#	p = runif(length(ind))
		#}
	p=p/sum(p)
	p
}

identify.outlier = function(alpha, beta,theta, theta2, var.epsilon, var.epsilonout, var.delta, var.deltaout, var.delta2, var.delta2out,t, y, lag, lag2, prob.out)
{
	outlier = rep(0, length(t))
	for(i in 1:length(t))
	{
		p = c(1,0)
		 if(i <= lag)
		 {
			 p[1] = dnorm(y[i], theta, sqrt(var.delta))*(1-prob.out)
			 p[2] = dnorm(y[i], theta, sqrt(var.deltaout))*prob.out
		 }
		 if(i > lag & i <= lag2)
		 {
			 p[1] = dnorm(y[i], alpha + beta*t[i], sqrt(var.epsilon))*(1-prob.out)
			 p[2] = dnorm(y[i], alpha + beta*t[i], sqrt(var.epsilonout))*prob.out
		 }
		 if(i > lag2)
		 {
			 p[1] = dnorm(y[i], theta2, sqrt(var.delta2))*(1-prob.out)
			 p[2] = dnorm(y[i], theta2, sqrt(var.delta2out))*prob.out

		 }
		p = p/sum(p)
		outlier[i] = (runif(1) < p[2])
	}
	outlier
}
densitymax = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, rho1, var.change=NULL)
			{
				
			
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T)) 
				if(change1 > 0 & !is.null(var.change))
				{
					s = s + dtnorm(change1, 0, sqrt(var.change), lower = t[lag1], upper = t[(lag1 + 1)], log = T)
				}
				if(change1 >0 & is.null(var.change))
				{
					s = s #+ dunif(change1, t[lag1], t[(lag1+1)], log = T)
				}
				s
			}

densitymax2 = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, rho2)
			{
				
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				
				
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T))
				if(change2 != t[length(t)])
				{
					s = s #+ dunif(change2, t[lag2], t[(lag2 + 1)], log = T)
				}
				s
			}

densityintegrate = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, var.change = NULL)
			{
				value = rep(0, length(change1))
				for(i in 1:length(change1))
				{
					
						tnew = c(rep(change1[i], lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
					
					p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
					if(change1[i] > 0 & !is.null(var.change))
					{
						p = p*dtnorm(change1[i], 0, sqrt(var.change), lower = 0, upper = t[(lag2 - 1)])
					}
					if(change1[i] > 0 & is.null(var.change))
					{
						p = p*dunif(change1[i], 0, change2)
					}
					value[i] = p
				}
				value
			}

densityintegrate2 = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out)
			{
				value = rep(0, length(change2))
				for(i in 1:length(change2))
				{
					
						tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2[i], (length(t)-lag2)))
					

				 	
				 	p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
				 	if(change2[i] != t[length(t)])
					{
					p = p*dunif(change2[i], change1, t[length(t)])
					}
					value[i] = p
				}
				value
			}




			
densitymax.lt = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, rho1, var.change=NULL)
			{
				
				if(lag1 != lag2){
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				}else{tnew = c(rep(change1, lag1), rep(change2, (length(t)-lag2)))}
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T)) 
				if(change1 > 0 & !is.null(var.change))
				{
					s = s + dtnorm(change1, 0, sqrt(var.change), lower = t[lag1], upper = min(t[(lag1+1)], change2), log = T)
				}
				if(change1 >0 & is.null(var.change))
				{
					s = s #+ dunif(change1, t[lag1], min(t[(lag1+1)], change2), log = T)
				}
				s
			}

densitymax2.lt = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, rho2, var.change = NULL)
			{
				if(lag1 != lag2){
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				}else{tnew = c(rep(change1, lag1), rep(change2, (length(t)-lag2)))}
				
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T))
				if(change2 != t[length(t)] & is.null(var.change))
				{
					s = s #+ dunif(change2, max(t[lag2], change1) , t[(lag2+1)], log = T)
				}
				if(change2 != t[length(t)] & !is.null(var.change))
				{
					s = s + dtnorm(change2, max(t), sqrt(var.change), lower = max(t[lag2], change1) , upper = t[(lag2+1)], log = T)
				}
					

				s
			}

densityintegrate.lt = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, var.change = NULL)
			{
				value = rep(0, length(change1))
				for(i in 1:length(change1))
				{
					if(lag1 != lag2){
						tnew = c(rep(change1[i], lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
					}else{tnew = c(rep(change1[i], lag1), rep(change2, (length(t)-lag2)))}
					p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
					if(change1[i] > 0 & !is.null(var.change))
					{
						p = p*dtnorm(change1[i], 0, sqrt(var.change), lower = 0, upper = change2)
					}
					if(change1[i] > 0 & is.null(var.change))
					{
						p = p*dunif(change1[i], 0, change2)
					}
					value[i] = p
				}
				value
			}

densityintegrate2.lt = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, var.change = NULL)
			{
				value = rep(0, length(change2))
				for(i in 1:length(change2))
				{
					if(lag1 != lag2){
						tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2[i], (length(t)-lag2)))
					}else{tnew = c(rep(change1, lag1), rep(change2[i], (length(t)-lag2)))}

				 	
				 	p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
				 	if(change2[i] != t[length(t)] & is.null(var.change))
					{
					p = p*dunif(change2[i], change1, t[length(t)])
					}
					if(change2[i] != t[length(t)] & !is.null(var.change))
					{
					p = p*dtnorm(change2[i], max(t), sqrt(var.change), lower = change1 , upper = max(t))
					}
					value[i] = p
				}
				value
			}
			
identify.outlier.lt = function(alpha, beta,theta, theta2, var.epsilon, var.epsilonout, var.delta, var.deltaout, var.delta2, var.delta2out,t, y, lag, lag2, prob.out)
{
	outlier = rep(0, length(t))
	for(i in 1:length(t))
	{
		 p = c(1,0)
		 if(i <= lag)
		 {
			 p[1] = dnorm(y[i], theta, sqrt(var.delta))*(1-prob.out)
			 p[2] = dnorm(y[i], theta, sqrt(var.deltaout))*prob.out
		 }
		 if(i > lag & i <= lag2)
		 {
			 p[1] = dnorm(y[i], alpha + beta*t[i], sqrt(var.epsilon))*(1-prob.out)
			 p[2] = dnorm(y[i], alpha + beta*t[i], sqrt(var.epsilonout))*prob.out
		 }
		 if(i > lag2)
		 {
			 p[1] = dnorm(y[i], theta2, sqrt(var.delta2))*(1-prob.out)
			 p[2] = dnorm(y[i], theta2, sqrt(var.delta2out))*prob.out

		 }
		p = p/sum(p)
		outlier[i] = (runif(1) < p[2])
	}
	outlier
}
densitymax.W = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, sh, sc)
			{
				
				if(lag1 != lag2){
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				}else{tnew = c(rep(change1, lag1), rep(change2, (length(t)-lag2)))}
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T)) 
				
				if(change1 >0 )
				{
					s = s + dweibull(change1, sh, sc, log = T)
				}
				s
			}

densitymax2.W = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, sh, sc)
			{
				max.t = max(t)
				if(lag1 != lag2){
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				}else{tnew = c(rep(change1, lag1), rep(change2, (length(t)-lag2)))}
				
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T))
				if(change2 != t[length(t)] )
				{
					s = s + dweibull((max.t - change2), sh, sc, log = T)
					#s = s + dweibull((change2-change1), sh, sc, log = T)
				}
					
				s
			}

densityintegrate.W = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, sh, sc)
			{
				value = rep(0, length(change1))
				for(i in 1:length(change1))
				{
					if(lag1 != lag2){
						tnew = c(rep(change1[i], lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
					}else{tnew = c(rep(change1[i], lag1), rep(change2, (length(t)-lag2)))}
					p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
					if(change1[i] > 0)
					{
						p = p*dweibull(change1[i], sh, sc)/pweibull(change2, sh, sc)
					}
					value[i] = p
				}
				value
			}

densityintegrate2.W = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, sh, sc)
			{
				max.t = max(t)
				value = rep(0, length(change2))
				for(i in 1:length(change2))
				{
					if(lag1 != lag2){
						tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2[i], (length(t)-lag2)))
					}else{tnew = c(rep(change1, lag1), rep(change2[i], (length(t)-lag2)))}

				 	
				 	p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
				 	if(change2[i] != t[length(t)])
					{
					p = p*dweibull((max.t - change2[i]), sh, sc)/pweibull((max.t-change1), sh, sc)
					#p = p*dweibull((change2[i]-change1), sh, sc)

					}
					
					value[i] = p
				}
				value
			}
			
			
moments.weibull=function(x,sh.prev)
{
	mle.weib = function(sh,x){
		
		mx = min(x)
		1/sh - sum(x^sh*log(x) - mx^sh*log(mx))/sum(x^sh-mx^sh) + mean(log(x))
	}
	
	diff = 10
	mx = min(x)
	while(diff > .001)
	{
	mx = min(x)
	# sh = sh.prev
	# f = 1/sh - sum(x^sh*log(x) - mx^sh*log(mx))/sum(x^sh-mx^sh) + mean(log(x))
	# fx = -1/sh^2 - (sum(x^sh*log(x)^2 - mx^sh*log(mx)^2)*sum(x^sh-mx^sh) - sum(x^sh*log(x) - mx^sh*log(mx))*sum(x^sh*log(x)-mx^sh*log(mx)))/sum(x^sh-mx^sh)^2
	# sh = (sh - f/fx)
	# diff = abs(sh - sh.prev)
	# sh.prev = sh
	
	 sh = sh.prev
	 sh = 1/(sum(x^sh*log(x) - mx^sh*log(mx))/sum(x^sh-mx^sh) - mean(log(x)))
	 diff = abs(sh - sh.prev)
	 sh.prev = sh
	}
	

	
	
	#sh = uniroot(mle.weib, interval = c(.1, 10), x)$root
	mx = min(x)
	sc = mean(x^sh - mx^sh)^(1/sh)

		
	dxx = sum((sh^2*(x/sc)^sh * log(x/sc)^2 + 1)/sh^2)
	dyy = -sum((sh - (sh)*(sh+1)*(x/sc)^sh)/sc^2)
	dxy = -sum(((x/sc)^sh *(sh*log(x/sc) + 1) - 1)/sc)
	v = solve(cbind(c(dxx,dxy), c(dxy, dyy)))
	return(list(m = c(sh, sc), v = v ))
}

# shape.weibull=function(x,sh.prev, sc)
# {
	# mle.weib = function(sh,x){
		
		# mx = min(x)
		# 1/sh - sum(x^sh*log(x) - mx^sh*log(mx))/sum(x^sh-mx^sh) + mean(log(x))
	# }
	
	# diff = 10
	# iter = 0
	# while(diff > .001)
	# {
	# # mx = min(x)
	# # sh = sh.prev
	# # f = 1/sh - sum(x^sh*log(x) - mx^sh*log(mx))/sum(x^sh-mx^sh) + mean(log(x))
	# # fx = -1/sh^2 - (sum(x^sh*log(x)^2 - mx^sh*log(mx)^2)*sum(x^sh-mx^sh) - sum(x^sh*log(x) - mx^sh*log(mx))*sum(x^sh*log(x)-mx^sh*log(mx)))/sum(x^sh-mx^sh)^2
	# # sh = (sh - f/fx)
	# # diff = abs(sh - sh.prev)
	# # sh.prev = sh
		# mx = min(x)
	# sh = sh.prev
	# fx = length(x)/sh + sum(log(x/sc)) - sum(sh*log(x/sc)*(x/sc)^sh)
	# fxx = -sum((sh^2*(x/sc)^sh * log(x/sc)^2 + 1)/sh^2)
	# sh = (sh - .5*fx/fxx)
	# diff = abs(sh - sh.prev)
	# sh.prev = sh
	# sc = (mean(x^sh))^(1/sh)
	# iter = iter + 1
	# print(iter)
	# }

	
	
	# #sh = uniroot(mle.weib, interval = c(.1, 10), x)$root
	# mx = min(x)
	# sc = mean(x^sh)^(1/sh)
	
	# lsh = log(sh)

	# dxx = sum((sh^2*(x/sc)^sh * log(x/sc)^2 + 1)/sh^2)
	
	# lxx = dxx*exp(2*log(sh))
	# return(list(m = sh, s = sqrt(1/dxx), lm = lsh, ls = sqrt(1/lxx)))
# }

scale.weibull=function(x,sh)
{
	mx = min(x)
	sc = mean(x^sh - mx^sh)^(1/sh)
	dyy = -sum((sh - (sh)*(sh+1)*(x/sc)^sh)/sc^2)
	return(list(m = sc, s = sqrt(1/dyy)))
}

shape.weibull <- function (x=NULL,n=NULL){
if(is.null(x))x <- c(7,12.1,22.8,23.1,25.7,26.7,29.0,29.9,39.5,41.9)
r <- length(x)
if(is.null(n)){n<-r}else{if(r>n||r<2){
return("x must have length r with: 2 <= r <= n")}}
xs <- sort(x)
# MPF 2015/02/25: survival attached by Depends
#if(!exists("survreg"))library(survival)
#tests whether survival package is loaded, if not, 
# then it loads survival
if(r<n){
statusx <- c(rep(1,r),rep(0,n-r))
dat.weibull <- data.frame(c(xs,rep(xs[r],n-r)),statusx)
14}else{statusx <- rep(1,n)
dat.weibull <- data.frame(xs,statusx)}
names(dat.weibull)<-c("time","status")
out.weibull <- survreg(Surv(time,status)~1,dist="weibull",data=dat.weibull)
alpha.hat <- exp(out.weibull$coef)
beta.hat <- 1/out.weibull$scale
parms <- c(alpha.hat,beta.hat)
sh = beta.hat
sc = alpha.hat
lsh = log(sh)

dxx = sum((sh^2*(x/sc)^sh * log(x/sc)^2 + 1)/sh^2)
	
lxx = dxx*exp(2*log(sh))
return(list(m = sh, s = sqrt(1/dxx), lm = lsh, ls = sqrt(1/lxx)))
}


grid.weibull = function(x=NULL,ub,size = 100, n=NULL)
{
if(is.null(x))x <- c(7,12.1,22.8,23.1,25.7,26.7,29.0,29.9,39.5,41.9)
r <- length(x)
if(is.null(n)){n<-r}else{if(r>n||r<2){
return("x must have length r with: 2 <= r <= n")}}
xs <- sort(x)
# MPF 2015/02/25, survival attached by Depends
#if(!exists("survreg"))library(survival)
#tests whether survival package is loaded, if not, then it loads # survival
if(r<n){
statusx <- c(rep(1,r),rep(0,n-r))
dat.weibull <- data.frame(c(xs,rep(xs[r],n-r)),statusx)
14}else{statusx <- rep(1,n)
dat.weibull <- data.frame(xs,statusx)}
names(dat.weibull)<-c("time","status")
out.weibull <- survreg(Surv(time,status)~1,dist="weibull",data=dat.weibull)
alpha.hat <- exp(out.weibull$coef)
beta.hat <- 1/out.weibull$scale
parms <- c(alpha.hat,beta.hat)
sh = beta.hat
sc = alpha.hat
m.sh = log(sh)
m.sc = log(sc)

dxx = sum((sh^2*(x/sc)^sh * log(x/sc)^2 + 1)/sh^2)

	
lxx = dxx*exp(2*log(sh))
dyy = -sum((sh - (sh)*(sh+1)*(x/sc)^sh)/sc^2)

lyy = dyy*exp(2*log(sc))
s.sh = 1/sqrt(lxx)
s.sc = 1/sqrt(lyy)

m.scpr = log(2)
s.scpr = log(10)/3
m.shpr = log(2)
s.shpr = log(10)/3

s.scpo = sqrt(1/(1/s.sc^2 + 1/s.scpr^2)) 
s.shpo = sqrt(1/(1/s.sh^2 + 1/s.shpr^2)) 
m.scpo = s.scpo^2*(m.sc/s.sc^2 + m.scpr/s.scpr^2)
m.shpo = s.shpo^2*(m.sh/s.sh^2 + m.shpr/s.shpr^2)

range.sc = m.sc+ c(-3,3)*s.sc
range.sc = c(max(range.sc[1], log(.1)), min(range.sc[2], log(50)))
range.sh = m.sh + c(-3,3)*s.sh
range.sh = c(max(range.sh[1], log(.5)), min(range.sh[2], log(5)))

grid = matrix(0, size, size)
grid.sc = seq(range.sc[1], range.sc[2], length.out = 100)
grid.sh = seq(range.sh[1], range.sh[2], length.out = 100)
sh.marg = rep(0, size)
sc.cond = grid
for(i in 1:size)
{
	for(j in 1:size)
	{
		sc = exp(grid.sc[j])
		sh  = exp(grid.sh[i])
		grid[i,j] = sum(dweibull(x, sh, sc, log = T))- sum(pweibull(ub, sh, sc, log.p = T)) + log(1/sc) #dnorm(log(sc), m.scpr, s.scpr, log = T) + dnorm(log(sh), m.shpr, s.shpr, log = T)+log(1/sc) + log(1/sh)
	}
		
}
M = max(grid)
grid = exp(grid - M)
grid = grid/sum(grid)
for(i in 1:size)
{
	sh.marg[i] = sum(grid[i,])
	sc.cond[i,] = grid[i,]/sh.marg[i]

}

sh.samp = sample(grid.sh, 1, prob = sh.marg)
ii = which(grid.sh == sh.samp)
sc.samp = sample(grid.sc, 1, prob = sc.cond[ii,])

sh.samp = exp(sh.samp)
sc.samp = exp(sc.samp)

return(list(shape = sh.samp, scale = sc.samp))
}


grid.exponential = function(x=NULL,ub,size = 200, n=NULL)
{

sh = 1
sc = mean(x)
m.sh = log(sh)
m.sc = log(sc)

dxx = sum((sh^2*(x/sc)^sh * log(x/sc)^2 + 1)/sh^2)

	
lxx = dxx*exp(2*log(sh))
dyy = -sum((sh - (sh)*(sh+1)*(x/sc)^sh)/sc^2)

lyy = dyy*exp(2*log(sc))
s.sh = 1/sqrt(lxx)
s.sc = 1/sqrt(lyy)


range.sc = m.sc+ c(-3,3)*s.sc
range.sc = c(max(range.sc[1], log(.1)), min(range.sc[2], log(50)))
#range.sh = m.sh + c(-3,3)*s.sh
#range.sh = c(max(range.sh[1], log(.5)), min(range.sh[2], log(5)))

grid = rep(0, size)
grid.sc = seq(range.sc[1], range.sc[2], length.out = size)

sc.cond = grid
for(j in 1:size)
{
	
		sc = exp(grid.sc[j])
		sh  = 1
		grid[j] = sum(dweibull(x, sh, sc, log = T))- sum(pweibull(ub, sh, sc, log.p = T)) + log(1/sc) #dnorm(log(sc), m.scpr, s.scpr, log = T) + dnorm(log(sh), m.shpr, s.shpr, log = T)+log(1/sc) + log(1/sh)
	
		
}
M = max(grid)
grid = exp(grid - M)
grid = grid/sum(grid)


sc.samp = sample(grid.sc, 1, prob = grid)

sh.samp = 1
sc.samp = exp(sc.samp)

return(list(shape = sh.samp, scale = sc.samp))
}






densitymax.LN = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, mu, sigma2)
			{
				
				if(lag1 != lag2){
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				}else{tnew = c(rep(change1, lag1), rep(change2, (length(t)-lag2)))}
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T)) 
				
				if(change1 >0 )
				{
					s = s + dnorm(log(change1), mu, sqrt(sigma2), log = T) + log(1/change1)
				}
				s
			}

densitymax2.LN = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, mu, sigma2)
			{
				max.t = max(t)
				if(lag1 != lag2){
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				}else{tnew = c(rep(change1, lag1), rep(change2, (length(t)-lag2)))}
				
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T))
				if(change2 != t[length(t)] )
				{
					s = s + dnorm(log((max.t - change2)), mu, sqrt(sigma2), log = T) + log(1/(max.t - change2)) 					
				}
					
				s
			}

densityintegrate.LN = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, mu, sigma2)
			{
				value = rep(0, length(change1))
				for(i in 1:length(change1))
				{
					if(lag1 != lag2){
						tnew = c(rep(change1[i], lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
					}else{tnew = c(rep(change1[i], lag1), rep(change2, (length(t)-lag2)))}
					p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
					if(change1[i] > 0)
					{
						p = p*dnorm(log(change1[i]), mu, sqrt(sigma2))*(1/change1[i]) /pnorm(log(change2),mu ,sqrt(sigma2))					}
					value[i] = p
				}
				value
			}

densityintegrate2.LN = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, mu, sigma2)
			{
				max.t = max(t)
				value = rep(0, length(change2))
				for(i in 1:length(change2))
				{
					if(lag1 != lag2){
						tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2[i], (length(t)-lag2)))
					}else{tnew = c(rep(change1, lag1), rep(change2[i], (length(t)-lag2)))}

				 	
				 	p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
				 	if(change2[i] != t[length(t)])
					{
					p = p* dnorm(log(max.t - change2[i]), mu, sqrt(sigma2))*(1/(max.t - change2[i])) 	/pnorm(log(max.t-change1),mu , sqrt(sigma2))
					

					}
					
					value[i] = p
				}
				value
			}
			
			

densitymax.N = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, mu, sigma2)
			{
				
				if(lag1 != lag2){
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				}else{tnew = c(rep(change1, lag1), rep(change2, (length(t)-lag2)))}
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T)) 
				
				if(change1 >0 )
				{
					s = s + dnorm(change1, mu, sqrt(sigma2), log = T)
				}
				s
			}

densitymax2.N = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, mu, sigma2)
			{
				max.t = max(t)
				if(lag1 != lag2){
					tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
				}else{tnew = c(rep(change1, lag1), rep(change2, (length(t)-lag2)))}
				
				s = sum(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out, log = T))
				if(change2 != t[length(t)] )
				{
					s = s + dnorm(max.t - change2, mu, sqrt(sigma2), log = T) 					
				}
					
				s
			}

densityintegrate.N = function(change1, change2, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, mu, sigma2)
			{
				value = rep(0, length(change1))
				for(i in 1:length(change1))
				{
					if(lag1 != lag2){
						tnew = c(rep(change1[i], lag1), t[(lag1+1):lag2], rep(change2, (length(t)-lag2)))
					}else{tnew = c(rep(change1[i], lag1), rep(change2, (length(t)-lag2)))}
					p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
					if(change1[i] > 0)
					{
						p = p*dnorm(change1[i], mu, sqrt(sigma2))				}
					value[i] = p
				}
				value
			}

densityintegrate2.N = function(change2, change1, lag1, lag2, alpha, beta, t, y, var.epsilon, var.epsilonout, out, mu, sigma2)
			{
				max.t = max(t)
				value = rep(0, length(change2))
				for(i in 1:length(change2))
				{
					if(lag1 != lag2){
						tnew = c(rep(change1, lag1), t[(lag1+1):lag2], rep(change2[i], (length(t)-lag2)))
					}else{tnew = c(rep(change1, lag1), rep(change2[i], (length(t)-lag2)))}

				 	
				 	p = prod(dnorm(y, alpha + beta*tnew, sqrt(var.epsilon)*(1-out) + sqrt(var.epsilonout)*out))
				 	if(change2[i] != t[length(t)])
					{
					p = p* dnorm(max.t - change2[i], mu, sqrt(sigma2))
					

					}
					
					value[i] = p
				}
				value
			}
			





