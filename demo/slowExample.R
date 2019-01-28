library(bhrcr)
data(pursat)
data(pursat_covariates)
data(posterior)
# The following takes a long time, so we calculated it ahead of time 
# and saved it as a data set.
# posterior <- clearanceEstimatorBayes(data = pursat, covariates = pursat_covariates, seed = 1234, 
#                                       detect.limit = 15, niteration = 25500, burnin = 500, thin = 100)

summary(posterior)
print(posterior)

# diagnostic plots of the Bayesian analysis
diagnostics(posterior)

# plot patients' profiles
plot(posterior)