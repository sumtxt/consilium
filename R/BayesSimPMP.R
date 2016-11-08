#' Simulates Posterior Predicted Probabilities for Proposal Passage from a Partial M-Probit Posterior 
#'
#' \code{BayesSimPMP} simulates posterior predicted probabilities for an estimated partial m-probit using Monte Carlo. 
#'
#' 
#' @param posterior (required) a posterior object from \code{BayesPMP()}. Can have multiple chains. 
#' @param data (required) a \code{data.frame} object that contains all RHS variables used in the \code{formula} to obtain 
#' 		 posterior using \code{BayesPMP()}. 
#' @param prpslid (required) the name for the consecutive numbered (\eqn{[1, J]}) integer variable that identifies 
#' 		each proposal uniquely in the data. 
#' @param gntid (required) the name for the consecutive numbered (\eqn{[1, M]}) integer variable that identifies 
#' 		each voting member in the data.
#' @param groups the name for the integer variable in \code{data} that identifies groups of proposals for which a varying intercept should be estimated.   
#' @param R (required) the number of members that have to agree to pass a proposal excluding the 
#'		members that have veto power. Can be a vector of length \code{n_distinct(prpslid)} or a single number. 
#' @param V the number of veto members that have to agree to pass a proposal. Can 
#'		can be a vector of length \code{n_distinct(prpslid)} or a single number. The function assumes that 
#'		the first V members in the data are the once with veto power. That is, after sorting the data for each  
#'		proposal by \code{gntid} the first \code{V} entries are used as they belong to members with 
#'		with veto power. 
#' @param q/pweights integer matrix of dimension \code{M x T} where \code{T} is the number of distinct voting weights. The vote weights sorted by gntid. 
#' @param q/pvote scalar or integer vector of length \code{J}. The threshold for the vote weights. 
#' @param q/pslct integer vector of length \code{J}. Each entry refers to the applicable column in the \code{q/pweights} for a particular proposal \code{j}. 
#' 
#' 
#' @param NSIM integer. The number of Monte Carlo simulations to run.
#' @param verbose integer number: use \code{0} for no output; use any positive natural number to report progress at each 
#' 		\code{verbose}-th iteration.
#' @param seed integer number for the seeding value.
#' 
#' @param postSamp number of posterior samples to be used in the calculations. 
#' @param postMean only calculate the predicted probabilities for the posterior mean? 
#'   
#'   
#' @details 
#' To obtain the posterior density for the predicted probability that proposal \code{j} passes, the function uses a 
#' a Monte Carlo algorithm. Conditional on the covariates for all voting members with respect to proposal \code{j} and a single draw from 
#' the posterior density of the coefficients, it simulates \code{NSIM} vote profiles. The number of 
#' adoption decisions implied by these vote profiles relative to \code{NSIM} gives an estimate about the probability for the passage of the 
#' \code{j^{th}} proposal. Iterating over all posterior draws gives the posterior density for the predicted probability that proposal \code{j} passes. 
#' 
#' The function runs in \code{C++}.
#' 
#' The function can not (yet) obtain posterior predicted probabilities if the posterior includes a varying intercept. 
#'
#' @return \code{matrix} object of size \code{max(prpslid)} times \code{niter(posterior) * nchain(posterior)}. 
#'
#' @examples 
#'  \dontrun{
#'  # Example 1: Using Simulated Data # 
#'  ###########################
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
#'  m1 <- BayesPMP(formula=y.agg ~ x1,R=R, prpslid="prpslid", gntid="gntid", data=data)
#'	
#'  # Simulate posterior predicted probability for first proposal (prpslid==1).
#'  m1sim <- BayesSimPMP(m1, R=R, prpslid="prpslid", gntid="gntid",data=data[data$prpslid==1,])
#'  hist(m1sim, breaks=100)
#'  }
#'
#' @references 
#'	Marbach, Moritz. 2016. `Analyzing Decision Records from Committees.'' Working Paper.
#' 
#' 
#' 
#' @export
BayesSimPMP <- function(posterior, data, prpslid, gntid, R, V=NULL, qweights = NULL, pweights = NULL, 
						qvote =NULL, pvote =NULL, pslct=NULL, qslct=NULL, method = "mc",  
						NSIM=500, verbose=500, seed=NULL, postSamp=NULL, postMean = FALSE){

	if ( attr(posterior, "type") != "consilium" ) stop("Posterior is not from BayesPMP().") 
	if ( attr(posterior, "groups") != 1 ) stop("Simulations with varying intercept not (yet) implemented.") 

	# Data sanity check 
	if (  is.numeric(data[,gntid]) != TRUE | is.numeric(data[,prpslid]) != TRUE ) stop("Supply prpslid and gntid as numeric.") 
	data$m <- unclass(as.factor(data[,gntid]))
	data$j <- unclass(as.factor(data[,prpslid]))

	# Sorting observations -- !! VERY IMPORTANT !! --
	data <- data[order(data$j, data$m),]
	constJ = max(data$j)
	maxM = max(data$m)

	formula <- attr(posterior, "call")
	ff <- formula( Formula(formula), rhs=1, lhs=0)
	X <- model.matrix(ff, data)	
	constK = ncol(X)
	X <- as.data.frame(X)
	X$j <- data$j

	Ms <- ddply(data, "j", function(x) nrow(x) )$V1
	if ( is.null(V) ) { Vs = rep(0,constJ) } else {
		if ( length(V) == 1 ) Vs = rep(V, constJ) else Vs <- V
		}
	if ( length(R)==1 ) Rs = rep(R,constJ) else Rs <- R

	if( !is.null(seed) ) set.seed(seed)
	
	
	posterior <- as.matrix(posterior)
	if ( !is.null(postSamp) ) {
		posterior <- posterior[sample(seq(1,nrow(posterior)),postSamp),]
		}
	if ( postMean == TRUE ) { 
		posteriormu <- colMeans(posterior)
		posterior <- matrix(posteriormu, 1, ncol(posterior) )
		}
	

	if ( is.null(qvote) & is.null(pvote) & is.null(qweights) & is.null(pweights) & method != "mc" ) { 

		if ( is.null(V) ) V=0L

		require(poibin)
		getPP <- function(data, betadraw,R,V,I) {
			pp <- pnorm(as.matrix(data) %*% betadraw)
			return( poibinprob(pp,R,V,I) )
			}

		ress <- ddply(X, "j", function(x) {
			x$j <- NULL
			out <- apply(posterior, 1, function(draw) getPP(data=as.matrix(x),betadraw=draw,R=R,V=V,I=maxM)  )
			return(out)
			})
		ress$j <- NULL 

		} else { 

		# Build default vote weights. 
		if ( is.null(pslct) )  pslct <- rep(1L, constJ)
		if ( is.null(qslct) )  qslct <- rep(1L, constJ)
		if ( is.null(qweights) )  qweights <- matrix(0L, maxM, 1)
		if ( is.null(pweights) )  pweights <- matrix(0L, maxM, 1)

		if ( is.null(qvote) )  qvote <- rep(0L, constJ)
		if ( length(qvote)==1 ) qvote = rep(qvote,constJ) 

		if ( is.null(pvote) )  pvote <- rep(0L, constJ)
		if ( length(pvote)==1 ) pvote = rep(pvote,constJ) 

		NSAMP <- as.integer(nrow(posterior))

		ress <- consiliumppc(posterior=posterior, X=as.matrix(X), 
					Rs=Rs, Vs=Vs, Is=Ms, qvote=qvote, pvote=pvote, 
					qweights=qweights, pweights=pweights, pslct=pslct-1, qslct=qslct-1,
					J=constJ, K=constK, 
					NSAMP=NSAMP, SIM=as.integer(NSIM), 
					VERBOSE=as.integer(verbose) )

		}

	rownames(ress) <- seq(1,constJ)

	return(ress)

	}
