#' Calculates the Rank of the Aggregate Design Matrix 
#'
#' \code{RankTest} calcualtes the rank of the aggregate design matrix of a dataset
#'
#'  
#' @param formula (required) a formula object of the form \code{ y | b ~ x1 + ... xK}. All \code{x} and \code{b} must be 
#' 		without missing values. Intercept is included automatically. For each proposal, the \code{b} vector must be identical 
#' 		across actors. \code{y |} is optional supplying known vote choices for a proposal.
#' 		Only used if all votes for a proposal are none-missing.  
#' @param data (required) a data.frame object that contains all data used in the \code{formula} and the variables specified
#' 		 in \code{prpslid} and \code{gntid}. 
#' @param prpslid (required) the name for the consecutive numbered (\eqn{[1, J]}) integer variable that identifies 
#' 		each proposal uniquely in the data. 
#' @param gntid (required) the name for the consecutive numbered (\eqn{[1, M]}) integer variable that identifies 
#' 		each actor in the data.
#' @param getAggData if \code{TRUE} the aggregate design matrix is returned 
#' 
#' @details 
#' This function is useful to evaluate a necessary condition for parameteric identification of the partial m-probit model. 
#' 
#' @examples 
#'  \dontrun{
#'
#' 	require(plyr)
#' 
#'  set.seed(10)
#'  J <- 250 	# proposal 
#'  I <- 10 	# actors
#'  R <- 6		# majority threshold
#'
#'  # Simualte roll call voting record 
#'  beta <- c(0,0.4)
#'  X <- data.frame(x0=1,x1=runif(J*I,-2,2))
#'  y <- rbinom(J*I, 1, pnorm(as.matrix(X) %*% beta))
#'
#'  # Bundle data with IDs
#'  data <- data.frame(gntid=sort(rep(seq(1,I), J)), 
#'  		prpslid=rep(seq(1,J), I), 
#'  		y, X)
#'  
#'  # Generate decision record 
#'  data <- ddply(data, "prpslid" ,function(x) { 
#'  		x$y.agg <- as.numeric(sum(x$y) >= R)
#'  		return(x)
#'  		})
#'
#'   RankTest(formula=y.agg ~ x1, data=data, prpslid="prpslid", gntid="gntid")  
#' }
#' 
#' @references 
#'	Marbach, Moritz. 2016. `Analyzing Decision Records from Committees.'' Working Paper.
#' 
#' 
#' @export
RankTest <- function(formula,data, prpslid, gntid, getAggData=FALSE, ...){

	# Data sanity check 
	if (  is.numeric(data[,gntid]) != TRUE | is.numeric(data[,prpslid]) != TRUE ) stop("Supply prpslid and gntid as numeric.") 
	data$m <- unclass(as.factor(data[,gntid]))
	data$j <- unclass(as.factor(data[,prpslid]))

	# Sorting observations -- !! VERY IMPORTANT !! --
	data <- data[order(data$j, data$m),]
	constJ = max(data$j)
	constM = max(data$m)

	f <- Formula(formula)
	f <- update(f, . ~ . - 1)
	mf <- model.frame(f,data, na.action="na.pass")
	X <- model.matrix(f, mf)

	X_ <- cbind(x0=1,X, j=data$j)
	X_ <- ddply(as.data.frame(X_), "j", function(x) colMeans(x) )
	X_$j <- NULL 

	res <- rankMatrix(X_, ...)

	cat("Required (at least): ", ncol(X_), "\n" )
	cat("Actual: ", as.numeric(res), "\n")

	if ( getAggData == TRUE) return(X_)

	}