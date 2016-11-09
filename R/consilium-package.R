#'	
#' Estimates multivariate probit models with partial
#' observability (partial m-probits) using a double-augmented Gibbs sampler and maximum likelihood. 
#' When using this package, please cite the paper and the package, see \code{citation("consilium")}. 
#' 
#' Feedback is very welcome! 
#'
#' 
#' @details
#' \tabular{ll}{
#'		Package: \tab consilium\cr
#'		Type: \tab Package\cr
#'		Version: \tab 0.0.0.9000\cr
#'		Date: \tab 2016-11-01\cr
#'		License: \tab  GPL-3\cr
#'		}
#'
#' @name consilium-package
#' 
#' @docType package
#' @aliases consilium
#' @title Estimating Multivariate Probit Models with Partial Observability
#' @author Moritz Marbach \email{marbach@ipz.uzh.ch}
#' @references 
#'	Marbach, Moritz. 2016. `Analyzing Decision Records from Committees.'' Working Paper.
#'
#' 
#' @import Rcpp 
#' @import RcppArmadillo
#' @import coda
#' @import Formula
#' @import poibin
#' @importFrom Matrix rankMatrix
#' @importFrom plyr ddply dlply
#' @importFrom texreg createTexreg screenreg
#' @importFrom combinat combn
#' @importFrom mnormt dmnorm
#' @importFrom utils methods
#' @importFrom methods new
#' @importFrom stats model.frame model.matrix optim pbinom pnorm rnorm runif terms update
#' 
#' @useDynLib consilium 
NULL

