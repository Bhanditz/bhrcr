library(bhrcr)
data(pursat)
data(pursat_covariates)
# Try it out for smaller number of iterations.
# For reproducibility, set seed = 1234.
out <- clearanceEstimatorBayes(data = pursat,covariates=pursat_covariates, seed=1234, 
                               detect.limit = 15, burnin=50, niteration=100, thin=10)

summary(out)
print(out)

# diagnostic plots of the Bayesian analysis
diagnostics(out)

# plot patients' profiles
plot(out)