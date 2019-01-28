#' @title Description of the Dataset \code{pursat_covariates.rda}
#' 
#' @description Individual level covariates of patients in the \code{pursat} dataset.
#' 
#' @keywords datasets
#' @name pursat_covariates
#' @usage data("pursat_covariates")
#' @format A data frame with following variables
#' \describe{
#'   \item{Sex}{A factor variable with two levels \code{F} and \code{M}}
#'   \item{agegroup}{21+ (21 years of age or older), or 21- (younger than 21 years)}
#'   \item{vvkv}{whether or not an individual was from  Veal Veng or Kranvanh}
#'   \item{HbE}{the number of alleles of Hemoglobin E variant}
#'   \item{athal}{the number of alleles of alpha-thalassaemia variant}
#'   \item{g6pd}{the number of alleles of G6PD deficient variant}
#'   \item{lnPf0}{Log initial parasite density}
#'   \item{year2010}{\code{TRUE} if 2010, \code{FALSE} if 2009}
#'   \item{group}{1 if group 1, 0 if group 2}
#' }
NULL