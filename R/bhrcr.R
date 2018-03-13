#' Bayesian Hierarchical Regression on Clearance Rates
#' 
#' "bhrcr" provides tools for calculating, analyzing, and visualizing parasite clearance rates 
#' in the presence of "lag" and "tail" phases through the use of a Bayesian hierarchical linear model.
#' The main function for the Bayesian hierarchical linear model is \link{clearanceEstimatorBayes}.
#' Also the function \link{calculatePCE} performs the method presented in Flegg et al (2011).
#' The hierarchical approach enables us to appropriately incorporate the uncertainty in both estimating 
#' clearance rates in patients and assessing the potential impact of covariates on these rates into the 
#' posterior intervals generated for the parameters associated with each covariate. Furthermore, it permits 
#' users to incorporate information about individuals for whom there exists only one observation time before 
#' censoring, which alleviates a systematic bias affecting inference when these individuals are excluded.
#' The detailed model and simulation study are presented in the paper "Bayesian Hierarchical Regression on 
#' Clearance Rates in the Presence of Lag and Tail Phases with an Application to Malaria Parasites" by Fogarty et al. (2015).
#' 
#' @references Flegg, J. A., Guerin, P. J., White, N. J., & Stepniewska, K. (2011). 
#' Standardizing the measurement of parasite clearance in falciparum malaria: the parasite clearance estimator. 
#' Malaria journal, 10(1), 339.
#' @references Fogarty, C. B., Fay, M. P., Flegg, J. A., Stepniewska, K., Fairhurst, R. M., & Small, D. S. (2015). 
#' Bayesian hierarchical regression on clearance rates in the presence of "lag" and "tail" phases 
#' with an application to malaria parasites. Biometrics, 71(3), 751-759.
#' 
#' @docType package
#' @name bhrcr-package
#' @aliases bhrcr
NULL