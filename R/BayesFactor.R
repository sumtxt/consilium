#' Computes Bayes Factor of Two Partial M-Probit Models
#'
#' \code{BayesFactor} computes the Bayes factor given two partial m-probit posteriors for which a marginal likelihood has been computed.
#'
#' @param m1 (required) partial m-probit model 1 estimated with chib95=TRUE
#' @param m2 (required) partial m-probit model 2 estimated with chib95=TRUE
#' @param log report BayesFactor on log-scale? 
#'
#' @return A scalar estimate of the Bayes factor. 
#' 
#' @references 
#'	Marbach, Moritz. 2016. `Analyzing Decision Records from Committees.'' Working Paper.
#' 
#' 	Chib, Siddhartha. 1995. Marginal Likelihood from the Gibbs Output. Journal of the American Statistical Association 90(432), 1313–1321.
#' 
#' @export
BayesFactor <- function(m1,m2, log=FALSE){
	if ( is.null( attr(m1, "logmarglik") ) ) stop("Model 1 does not contain a marginal likelihood estimate. Re-estimate with chib95=TRUE.")
	if ( is.null( attr(m2, "logmarglik") ) ) stop("Model 2 does not contain a marginal likelihood estimate. Re-estimate with chib95=TRUE.")
	logBF <- attr(m1, "logmarglik") - attr(m2, "logmarglik")
	logBF <- as.numeric(logBF)
	if (log==TRUE) return(logBF)
	else return(exp(logBF))
	}