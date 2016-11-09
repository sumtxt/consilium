#' Estimates Bayesian Partial M-Probit via Maximum Likelihood.
#'
#' \code{MaxLikPMP} obtains MLE estimates for a partial m-probit (pmp). 
#'
#' @param formula (required) a formula object of the form \code{b ~ x1 + ... xK}. All \code{x} and \code{b} must be 
#' 		without missing values. Intercept is included automatically. For each proposal, the \code{b} vector must be identical 
#' 		across actors. 
#' @param data (required) a data.frame object that contains all data used in the \code{formula} and the variables specified
#' 		 in \code{prpslid} and \code{gntid}. 
#' @param prpslid (required) the name for the consecutive numbered (\eqn{[1, J]}) integer variable that identifies 
#' 		each proposal uniquely in the data. 
#' @param gntid (required) the name for the consecutive numbered (\eqn{[1, M]}) integer variable that identifies 
#' 		each actor in the data. 
#' @param I integer. The number of actors. Only required if \code{max(gntid)==1},
#' @param V integer. The number of actors that have have to agree to pass a proposal. The function 
#'		assumes that the first V actors in the data are the once with veto power. That is, after sorting
#'		the data for each proposal by \code{gntid} the first \code{V} entries are used as they 
#'		belong to actors with veto power. 
#' @param R (required) integer. the number of actors that have to agree to pass a proposal excluding the 
#'		actors that have veto power. 
#' @param q/pweights integer matrix of dimension \code{M x 1}. The vote weights sorted by gntid. 
#' @param q/pvote scalar. The threshold for the vote weights. 
#' @param betastart a vector of starting value vectors for \eqn{\beta}. It is advised to supply starting values. 
#' @param whichloglik if set to \code{"comb"} the function will not make any attempt to find a more 
#' 		efficient representation of the likelihood function (see notes below)
#' @param verbose logical value. Use \code{TRUE} for progress report, \code{FALSE} for no output. 
#' @param method string value. Which of \code{optim()}'s optimization methods to use? 
#' @param control passed on to \code{optim()}.
#' @param behahat a matrix of size \code{K x H} of known \code{K} coefficient values for the model 
#' 		for which the \code{H} likelihood values are computed
#' 
#' 
#' 
#' @details 
#' Let \eqn{x_{ij}} be a vector of length \eqn{K} that collects all observed covariates for one of \eqn{M} 
#' actors and a proposal \eqn{j}. Let \eqn{y_{ij}} be the unobserved vote choice of actor \eqn{i} for proposal \eqn{j}. 
#' Let \eqn{b_j} be the observed vector of length \eqn{J} that collects the observed decisions made by the \eqn{M} actors
#'  according to a q-voting rule with threshold \eqn{R}. The model takes the following form: 
#'  
#' \deqn{Prob(b_j=0) = Prob( \sum_{i=1}^M y_{ij} < R) }
#' \deqn{Prob(y_{ij}=0) = \phi(x_{ij}\beta)}
#'
#' where \eqn{\phi()} is the standard univariate normal distribution and \eqn{\beta} is the parameter vector of interest.
#' 
#' 
#' Different to \code{BayesPMP(.)} this function can not handle proposal-specific voting rules, partially observed 
#' voting records or varying intercepts.  
#' 
#' 
#' Unless \code{betahat} is supplied, the function optimizes the implied likelihood using \code{optim()}. 
#' If \code{betahat} is supplied it computes the likelihood given using \code{betahat} as coefficients. 
#' 
#' To make the computation as efficient as possible, it attempts to take advantage of two special cases of the partial m-probit likelihood. 
#' 
#' a) If \code{max(gntid)==1} and \code{p/qweights==NULL} the function assumes that all 
#' covariate values are the same across actors for each proposal. In this case the likelihood 
#' is equivalent to a Bernoulli likelihood with a special (voting-rule depended) link function that Marbach (2016) refers to as 'B-link'. 
#' 
#' b) If \code{p/qweights==NULL}, the likelihood is the product of Poisson’s Binomial distribution functions 
#' which can be efficiently calculated using the discrete Fourier transform of the distribution's characteristic 
#' function (Hong 2013). 
#' 
#' 
#' If none of the special cases applies the likelihood is computed by summing across the probabilities of all 
#' potential voting coalition implied by the voting rule. Even for moderately sized committees 
#' this computation is very time- and memory-intensive. To force the function to not attempt to detect a special 
#' case set \code{whichloglik} to \code{"comb"}.
#' 
#'  
#' @return \code{maxlikpmp} object with the following slots: 
#' \itemize{
#' 	\item \code{n} the number of observation for the RHS of the \code{formula}
#' 	\item \code{coef} estimated coefficients
#' 	\item \code{coefnames} the corresponding variable names for the coefficients 
#' 	\item \code{hessian} the estimated hessian
#' 	\item \code{form} the formula used 
#' 	\item \code{code} the \code{optim()} convergence code 
#' 	\item \code{iter} the number of iterations from \code{optim()}
#' 	\item \code{I} the number of actors 
#' 	\item \code{R} the majority threshold 
#' 	\item \code{V} the number of actors with veto power
#' }
#'
#' The standard R functions \code{vcov}, \code{coef}, \code{logLik} and \code{summary} can be used on any \code{maxlikpmp} object.
#' 
#' 
#' @seealso \code{\link{optim}}
#'
#' @examples 
#'  \dontrun{
#'  # Example 1: Using Simulated Data # 
#'  ###########################
#'
#'  require(plyr)
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
#'  # Estimate partial m-probit 
#'  m1 <- MaxLikPMP(formula=y.agg ~ x1, prpslid="prpslid", gntid="gntid", data=data, I=I, R=R)
#'
#'  summary(m1)
#' 
#'  # Plot log-likelihood function
#'  require(lattice)
#'
#'  beta0 <- beta1 <- seq(-2,2,length=30)
#'  beta <- expand.grid(beta0,beta1)
#' 
#'  logliks <- MaxLikPMP(formula=y.agg ~ x1, prpslid="prpslid", gntid="gntid",
#'				     data=data, I=I, R=R, betahat=as.matrix(beta))
#'
#'  logliks <- data.frame(loglik = logliks * (-1), beta0=beta[,1], beta1=beta[,2])
#'
#'  wireframe(loglik ~ beta0 + beta1, data=logliks, 
#'		    shade=TRUE, main="Log-Likelihood", scales=list(arrows=FALSE) , zlim=c(-500,0) )
#'
#'
#'  # Example 2: Using Simulated Data # 
#'  ###########################
#'
#'  No variation across actors for a proposal. Uses B-Link. 
#'
#'  require(plyr)
#'
#'  set.seed(10)
#'  J <- 250 	# proposal 
#'  I <- 10 	# actors
#'  R <- 6		# majority threshold
#'
#'  # Simualte recorded votes
#'  beta <- c(0,0.4)
#'  X <- data.frame(x0=1,x1=rep(runif(J,-2,2), I) ) # changed compared to Ex.1
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
#'  # Estimate partial m-probit 
#'  m2 <- MaxLikPMP(formula=y.agg ~ x1, prpslid="prpslid", gntid="gntid", data=data[data$gntid==1,], I=I, R=R)
#'
#'  summary(m2)
#'
#'  }
#'
#'
#' @references 
#'	Marbach, Moritz. 2016. `Analyzing Decision Records from Committees.'' Working Paper.
#' 
#'  Hong, Yili. 2013. On Computing the Distribution Function for the Poisson Binomial Distribution. Computational Statistics & Data Analysis 59, 41–51.
#' 
#' @export
MaxLikPMP <- function(formula,data, prpslid, gntid, R, I=NULL, V=NULL, qweights = NULL, pweights = NULL, qvote =NULL, pvote =NULL, 
	betastart=NULL,verbose=TRUE, method="BFGS", control = list(), whichloglik=NULL, betahat = NULL){

	I_ = I
	V_= V
	R_= R
	qweights_ = qweights
	pweights_ = pweights
	qvote_ =qvote
	pvote_ =pvote

	# Data sanity check 
	if (  is.numeric(data[,gntid]) != TRUE | is.numeric(data[,prpslid]) != TRUE ) stop("Supply prpslid and gntid as numeric.") 
	data$m <- unclass(as.factor(data[,gntid]))
	data$j <- unclass(as.factor(data[,prpslid]))

	# Sorting observations -- !! VERY IMPORTANT !! --
	data <- data[order(data$j, data$m),]
	constJ = max(data$j)
	constM = max(data$m)

	f <- Formula(formula)
	mf <- model.frame(f,data, na.action="na.pass")

	if( length(f)[1] == 2 ){
		stop("MaxLikPMP can not handle partially observed roll call voting records.")
	} else {
		yagg <- as.integer(model.part(f, mf, lhs = 1, drop=TRUE))
		if ( sum(is.na(yagg)) > 0 ) stop(" Second part of RHS can not have 
			missing values if no information is supplied on the first part.")
		}	

	X <- model.matrix(f, mf)	
	constK = ncol(X)
	if( is.null(betastart) ) betastart <- rnorm(constK,0,0.01)	
	
	if ( length(R) > 1 | length(qvote) > 1 | length(pvote) > 1  ) stop("MaxLik requires decision-invariant voting rules.")		
	if ( is.null(V) ) V=0L

	# Choose Optimization method #
	use = ""
	if ( is.null(qvote) & is.null(pvote) & is.null(qweights) & is.null(pweights) ) use = "poibin"
	if (  constM == 1 & is.null(qweights) & is.null(pweights)  ) use = "blink"
	if (  use == "" ) use = "comb"
	if ( !is.null(whichloglik) ) use = whichloglik
	

	# USE BLINK LIKELIHOOD #
	###################

	if ( use == "blink" ){

		if( verbose==TRUE ) { cat("Calculate B-Link likelihood ... \n") }
		if ( !is.null(betahat) ) {
			if ( is.matrix(betahat) & ncol(betahat) == constK ) {
				values <- apply(betahat, 1, function(x) {
					cat(".")
					return(blinklikpmp(x, y=yagg, X=X, I=I, R=R,V=V))
					}) 

				attr(values, "lik") <- use
				return(values)
				}
			}

		if( verbose==TRUE ) cat("Optimizing B-Link liklihood ...  \n")
		if ( is.null(I) ) stop("B-Link requires supplying I.")
		out.optim <- optim(betastart, blinklikpmp, y=yagg, X=X, I=I, R=R,V=V,
			method=method, hessian=TRUE, control=control)
		}


	# PREPARE FOR PB / COMB LIKELIHOOD #
	#############################

	if ( use == "poibin" | use == "comb"){
		
		# Build default vote weights #
		if ( is.null(qweights) )  qweights <- matrix(0L, constM, 1)
		if ( is.null(pweights) )  pweights <- matrix(0L, constM, 1)
		
		if ( is.null(qvote) )  qvote <- 0L
		if ( is.null(pvote) )  pvote <- 0L

		# Reshape the data
		rdata <- as.data.frame(cbind(X, j=data$j))

		X_0 <- dlply(rdata[yagg==0,], "j", function(x){ 
			x$j <- NULL
			return(as.matrix(x)) 
			})

		X_1 <- dlply(rdata[yagg==1,], "j", function(x){ 
			x$j <- NULL
			return(as.matrix(x)) 
			})
		}

	# USE PB LIKELIHOOD  #
	#################

	if ( use == "poibin" ){

		# Valuation 
		if ( !is.null(betahat) ) {
			if ( is.matrix(betahat) & ncol(betahat) == constK ) {
				if( verbose==TRUE ) { cat("Calculate partial m-probit likelihood via poibin ... \n") }
				values <- apply(betahat, 1, function(x){
				 	cat(".")
				 	return(poibinlikpmp(x, X_0=X_0,X_1=X_1, R=R, V=V, I=constM, verbose=FALSE))
				 	})
				attr(values, "lik") <- use
				return(values)
				}
			}

		# Optimization
		if( verbose==TRUE ) { cat("Optimizing partial m-probit likelihood via poibin ... \n") }
		out.optim <- optim(betastart, poibinlikpmp, X_0=X_0,X_1=X_1, R=R, V=V, I=constM,
				verbose=verbose, method=method, hessian=TRUE, control=control)
		} 

	# USE COMBINATORIAL LIKELIHOOD #
	##########################

	if ( use == "comb" ){
		if( verbose==TRUE ) cat("Calculating V+/V- sets ... \n")
	
		# Calculate possible coalitions 
		cols <- mapply(function(x) t(combn(seq(1,constM),x, fun=constr,a=constM)), seq(1,constM), SIMPLIFY=FALSE)
		cols <- do.call(rbind, cols)
		
		# Separate into accept / reject coalitions depending on voting rule	
		st <- apply(cols, 1, function(x) check_vprofile(x,V=V, R=R, I=constM, 
			qvote=qvote, pvote=pvote, qweights=qweights, pweights=pweights) )
		st <- as.logical(st)

		yescols <- matrix(cols[st,], ncol=constM)
		nocols <- matrix(cols[!st,], ncol=constM)
		nocols <- rbind(nocols, rep(FALSE, constM))
		
		if( verbose==TRUE ) {
			cat(" |V-|: ", nrow(nocols), "\n")
			cat(" |V+|: ", nrow(yescols), "\n")
			}

		if ( !is.null(betahat) ) {
			if ( is.matrix(betahat) & ncol(betahat) == constK ) {
				if( verbose==TRUE ) { cat("Calculate partial m-probit likelihood ... \n") }
				values <- apply(betahat, 1, function(x) {
					cat(".")
					return(comblikpmp(x, X_0=X_0,X_1=X_1, nocols=nocols, yescols=yescols, verbose=FALSE))
					}) 
				attr(values, "lik") <- use
				return(values)
				}
			}

		if( verbose==TRUE ) {
			cat("Optimizing partial m-probit likelihood... \n")
			}
		
		out.optim <- optim(betastart, comblikpmp, X_0=X_0,X_1=X_1, nocols=nocols, yescols=yescols, 
			verbose=verbose, method=method, hessian=TRUE, control=control)
		}

	
	# EXPORT #
	########

	if( verbose==TRUE ) cat("Convergence: ", out.optim$convergence, 
				 "; iter: ", as.vector(out.optim$counts)[1], "\n" )


	res <- new("maxlikpmp", 
		loglik = -out.optim$value,
		n = nrow(X),
		coef = out.optim$par,
		coefnames = colnames(X), 
		hessian = solve(out.optim$hessian),
		form = formula,
		code = paste(out.optim$convergence),
		iter =as.vector(out.optim$counts)[1], 
		lik = use, 
		I = I_, 
		R = R_, 
		V = V_, 
		qweights = qweights_, 
		pweights = pweights_,
		qvote =qvote_,
		pvote =pvote_
		)

	return(res)
	}