#' Estimates Bayesian Partial M-Probit via Gibbs Sampling.
#'
#' \code{BayesPMP} obtains a set of draws from the posterior distribution of a (mixed) partial m-probit. 
#'
#'
#'
#' @param formula (required) a formula object of the form \code{ y | b ~ x1 + ... xK}. All \code{x} and \code{b} must be 
#' 		without missing values. For each proposal, the \code{b} vector must be identical 
#' 		across members. \code{y |} is optional supplying known vote choices for a proposal.
#' 		Only used if all votes for a subset of proposal are none-missing.  
#' @param data (required) a \code{data.frame} object that contains all data used in \code{formula} and the variables specified
#' 		 in \code{prpslid} and \code{gntid}. 
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
#' @param b0 the prior mean of length \eqn{K+1}. If \code{b0=NULL} the function assumes \code{b0=0}.
#' @param invB0 the prior precision matrix of dimension \eqn{(K+1) \times (K+1)}. If \code{invB0=NULL}
#' 		the function assumes \code{invB0=solve(diag(K+1)*1000)} where \code{K} is the number of coefficients implied by \code{formula}. 
#' @param e0/f01  the parameters for the inverse-gamma prior density if a varying intercept is included (\code{groups} is not \code{NULL})
#' @param chib95 calculate the marginal likelihood using the method of Chib (1995)?
#' @param betastart the list of starting value vectors for \eqn{\beta}. 
#' @param burnin the integer number of samples used as burn-in.
#' @param ngibbs the integer number of samples from the posterior distribution.
#' @param thin the n-th posterior draw to be recorded. Must be a number that yields a 
#' 		positive integer when \code{ngibbs/thin}. 
#' @param verbose integer number: use \code{0} for no output; use any positive natural number to report progress at each 
#' 		\code{verbose}-th iteration.
#' @param chains integer number for the number of chains.
#' @param seed integer number for the seeding value.
#' @param adapt if \code{TRUE} (default) the Gibbs sampler uses the more efficient algorithm 2 in Marbach (2016) otherwise algorithm 1.
#' @param rate schedule for the adaption if \code{adapt=TRUE}.
#' @param step \eqn{\epsilon}-parameter if \code{adapt=TRUE}. 
#' @param monitor if \code{adapt=TRUE} prints adaptation steps for each draw.  
#'  
#'   
#'   
#' @details 
#' Let \eqn{x_{ij}} be a vector of length \eqn{K} that collects all observed covariates for one of \eqn{M} 
#' members and a proposal \eqn{j}. Let \eqn{y_{ij}} be the unobserved vote choice of member \eqn{i} for proposal \eqn{j}. 
#' Let \eqn{b_j} be the observed vector of length \eqn{J} that collects the observed decisions made by the \eqn{M} members
#'  according to a q-rule with threshold \eqn{R}. The model takes the following form: 
#'  
#' \deqn{Prob(b_j=0) = Prob( \sum_{i=1}^M y_{ij} < R) }
#' \deqn{Prob(y_{ij}=0) = \phi(x_{ij}\beta)}
#'
#' where \eqn{\phi()} is the standard univariate normal distribution and \eqn{\beta} is the parameter vector of interest.
#' The prior density for \eqn{\beta} is a multivariate normal with user-specified prior mean vector 
#' and precision matrix. The function here simulates draws from the posterior density of \eqn{\beta}. 
#' 
#' Notice, that \code{BayesPMP}  can handle other voting rules than a q-rules with proposal-invariant threshold \eqn{R}.
#' The threshold and the number of members with veto powers are allowed to vary across proposals. Weighted voting rules are
#' also supported. 
#' 
#' If a varying intercept is included, it is assumed to be drawn from a normal density with a variance that has an inverse-gamma prior density. 
#' 
#' See the references for the derivation of the posterior and a description of the Gibbs sampler. The default (\code{adapt=TRUE}) uses 
#' algorithm 2 in Marbach (2016). The function runs in \code{C++}.
#' 
#' The marginal likelihood (\code{chib95=TRUE}) can not (yet) be calculated if a) a varying intercept is included (\code{groups} is \code{NULL}), 
#' b) the voting rule varies across decisions or c) with a partially observed voting record.
#' 
#' See \code{coda} documentation for help analyzing posterior samples. It is important to assess the convergence of the chains. 
#'
#'  
#' @return \code{coda} object of \code{ngibbs/thin} draws from the posterior distribution of the coefficients 
#'  	including the intercept. The \code{coda} object has 3 non-standard attributes: \code{type} which is set 
#' 	to \code{consilium}, \code{formula} which is  set to \code{formula} used in the function syntax 
#' 	and \code{timecode} recording the \code{Sys.time()} when the function finished. If \code{chib95=TRUE} the marginal 
#'  likelihood estimate is also included. 
#'
#' #' @seealso \code{\link{summary.mcmc}}, \code{\link{plot.mcmc}} and other coda-functions.
#'
#' @examples 
#'  \dontrun{
#'  # Example 1: q-rule # 
#'  ###########################
#'
#'  set.seed(10)
#'  J <- 250 	# proposal 
#'  I <- 10 	# members
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
#'  # Generate selected individual voting record 
#'  data$y.slct <- ifelse(data$prpslid <= 10, data$y, NA)
#'
#'  # Estimate partial m-probit with observed votes 
#'  m2 <- BayesPMP(formula=y.slct | y.agg ~ x1,R=R, prpslid="prpslid", gntid="gntid", data=data)
#'  
#'  summary(m1)
#'  summary(m2)
#'
#'  # Example 2: Weighted q-rule 
#'  #############################
#'
#'  rm(list=ls())
#'  set.seed(10)
#'  J <- 250 	# proposal 
#'  I <- 10 	# members
#'  R <- 6		# majority threshold
#'  voteweights <- matrix(0L, I, 1)
#'  voteweights[1:3,1] <- 10
#'  qvote <- 20
#'  qslct <- rep(1, J)
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
#'  		x$y.agg <- as.numeric( sum(x$y) >= R & (x$y %*% voteweights) >= qvote )
#'  		return(x)
#'  		})
#'  
#'  # Estimate partial m-probit 
#'  m2 <- BayesPMP(formula=y.agg ~ x1,R=R, prpslid="prpslid", 
#' 		gntid="gntid", data=data, qweights=voteweights, qvote=qvote, qslct=qslct)
#'  summary(m2)
#'  }
#'
#'
#' @references 
#'	Marbach, Moritz. 2016. `Analyzing Decision Records from Committees.'' Working Paper.
#' 
#' 	Chib, Siddhartha. 1995. Marginal Likelihood from the Gibbs Output. Journal of the American Statistical Association 90(432), 1313–1321.
#' 
#' 
#' 
#' 
#' @export
BayesPMP <- function(formula, data, prpslid, gntid, groups=NULL, R, V=NULL, qweights = NULL, pweights = NULL, qvote =NULL, pvote =NULL, pslct=NULL, qslct=NULL,
								  b0=NULL, invB0=NULL, e0=0.1, f0=0.1, chib95=FALSE, betastart=NULL, omegastart=NULL, 
								  burnin=500, ngibbs=2000, verbose=500, thin=1, chains=2,
								  seed=42, rate=200, step=0.05, adapt=TRUE, monitor = FALSE){

	# Data sanity check 
	if (  is.numeric(data[,gntid]) != TRUE | is.numeric(data[,prpslid]) != TRUE ) stop("Supply prpslid and gntid as numeric.") 
	data$m <- unclass(as.factor(data[,gntid]))
	data$j <- unclass(as.factor(data[,prpslid]))

	# Sorting observations -- !! VERY IMPORTANT !! --
	data <- data[order(data$j, data$m),]
	constJ = max(data$j)
	maxM = max(data$m)

	# Grouping 
	if ( !is.null(groups) ){
		alpha_idx <- as.factor(data[,groups])
		alpha_labs <- paste("(a_", levels(alpha_idx), ")",sep="")
		alpha_idx <- unclass(alpha_idx)
		G <- max(alpha_idx) 
		nG <- table(alpha_idx, useNA="no")
		if (chib95==TRUE) stop("Chib95=TRUE only works with groups=NULL.") 
		}
	else{
		G <- 0 
		nG <- 1
		alpha_idx <- 1
		}

	f <- Formula(formula)
	mf <- model.frame(f,data, na.action="na.pass")

	if( length(f)[1] == 2 ){
		yind <- as.integer(model.part(f, mf, lhs = 1, drop=TRUE))
		yagg <- as.integer(model.part(f, mf, lhs = 2, drop=TRUE))
		if ( sum(is.na(yind) & is.na(yagg)) ) stop("There are missing values in 
			first and second part of the RHS simultaneously.")
		yind_name <- attr(terms(f,rhs=0), "term.labels")[1]
		miss <- ddply(data, "j", function(x) as.integer(sum(is.na(x[,yind_name]))==nrow(x))  )$V1
	} else {
		yagg <- as.integer(model.part(f, mf, lhs = 1, drop=TRUE))
		if ( sum(is.na(yagg)) > 0 ) stop(" Second part of RHS can not have missing 
			values if no information is supplied on the first part.")
		yind <- NA
		miss <- rep(1, constJ)
		}	

	X <- model.matrix(f, mf)
	constK = ncol(X)

	
	# Check Gibbs Parameter #
	###################

	if( !is.null(seed) ) set.seed(seed)

	if ( as.integer(burnin) < 0) stop("Check: 0 < burnin < Inf!")
	if ( as.integer(verbose) < 0) stop("Check: 0 <= verbose < Inf!")
	if ( as.integer(ngibbs) < 0) stop("Check: 0 < ngibbs < Inf!")
	if ( as.integer(thin) < 1  ) stop("Check: 1 < thin < burnin!")
	if ( !is.naturalnumber(ngibbs/thin) ) stop("Check: !is.integer(ngibbs/thin)")

	# Make prior #
	###########

	if( is.null(b0) ) b0 <- rep(0,constK)
	if( is.null(invB0) ) invB0 <- solve(diag(constK)*100)

	# Starting Values #
	###############

	if ( is.null(omegastart) ){
		omegastart <- runif(chains,0.5,1)	
		}
	if ( !is.null(omegastart) & length(omegastart) != chains ) stop("Supply a list of omega startvalues of length = to #chains.")

	if ( is.null(betastart) ){
		betastart <- list()
		for(i in 1:chains) betastart[[i]] <- rnorm(constK,0,0.01)	
		}
	if ( !is.null(betastart) & length(betastart) != chains ) stop("Supply a list of beta startvalues of length = to #chains.")
	
	# Build Voting Rules #
	################

	Ms <- ddply(data, "j", function(x) nrow(x) )$V1

	V_ <- V
	if ( is.null(V) ) { Vs = rep(0,constJ) } else {
		if ( length(V) == 1 ) Vs = rep(V, constJ) else Vs <- V
		}
	R_ <- R
	if ( length(R)==1 ) Rs = rep(R,constJ) else Rs <- R

	qweights_ <- qweights
	pweights_ <- pweights
	if ( is.null(qweights) )  qweights <- matrix(0L, maxM, 1)
	if ( is.null(pweights) )  pweights <- matrix(0L, maxM, 1)
	
	qvote_ <- qvote
	if ( is.null(qvote) )  qvote <- rep(0L, constJ)
	if ( length(qvote)==1 ) qvote = rep(qvote,constJ) 

	pvote_ <- pvote
	if ( is.null(pvote) )  pvote <- rep(0L, constJ)
	if ( length(pvote)==1 ) pvote = rep(pvote,constJ) 

	if ( is.null(pslct) )  pslct <- rep(1L, constJ)
	if ( is.null(qslct) )  qslct <- rep(1L, constJ)

	if ( 	length(f)[1] == 2 & 
	 	is.null(pslct) & 
	 	is.null(qslct) & 
	 	length(unique(Ms)) & 
	 	length(R_) != 1 &
	 	length(unique(V)) & chib95==TRUE ) {
	 	stop("Chib95=TRUE only works with decision-invariant majority voting 
	 		rules and without partial observed recorded votes.") 
	 	}

	 poibin = is.null(qvote_) & is.null(pvote_) & is.null(qweights_) & is.null(pweights_) & is.null(V_)

	# BEGIN SAMPLING #
	###############

	if (verbose!=0) cat("Found max. ", max(Ms), "actors, ", 
		constJ, " proposals, ", constK-1, "covariates.\n")

	# Run
	ress <- list()

	for(i in 1:chains){
		if (verbose!=0) cat("\n== Chain ", i, " ==\n")
		ress[[i]] <- consiliumcA(X=X, z=yagg, y=yind , K=constK, Is=Ms, J=constJ, R=Rs, V=Vs, 
				qweights = qweights, pweights=pweights, qvote = qvote, pvote = pvote, pslct=pslct-1, qslct=qslct-1,
				miss=miss, NGIBBS=as.integer(ngibbs), NTHIN=as.integer(thin), NBURN=as.integer(burnin), 
				VERBOSE=as.integer(verbose), 
				adapt=as.integer(adapt), rate=as.integer(rate), step=as.double(step), VERBOSE2=as.integer(monitor), 
				b0=b0, invB0=invB0, betastart=betastart[[i]], omegastart=omegastart[i], G=G, 
				alpha_idx=alpha_idx-1, nG=nG, e0=e0, f0=f0)
		if (verbose!=0) cat("\n")
		}

	if ( G > 0 ){
		betachains <- lapply(ress, function(x) {
			tmp <- cbind(t(x$beta),x$omega,t(x$alpha))
			o <- as.mcmc(tmp) 
			colnames(o) <- c(colnames(X),"omega", alpha_labs)
			return(o)
			})
	} else {
		betachains <- lapply(ress, function(x) {
			o <- as.mcmc(t(x$beta)) 
			colnames(o) <- colnames(X)
			return(o)
			})
	}

	betachains <- as.mcmc.list( betachains )

	attr(betachains, "type") <- "consilium"
	attr(betachains, "call") <- formula
	attr(betachains, "groups") <- nG
	attr(betachains, "timecode") <- Sys.time()

	if (chib95 ==TRUE ){

		if (verbose!=0)  cat("\n Calculate Marginal Likelihood ... \n")
		#
		bmuchains <- lapply(ress, function(x) {
			o <- as.mcmc(t(x$bmu)) 
			colnames(o) <- colnames(X)
			return(o)
		})

		bmuchains <- as.mcmc.list( bmuchains )
		bmuchains <- as.matrix(bmuchains)

		# 
		bsd <- solve(invB0 + t(X) %*% X) 

		#
		betahat <- apply(as.matrix(betachains), 2, mean)

		#
		loglik <-  (-1) * MaxLikPMP(formula=formula, prpslid=prpslid, gntid=gntid, data=data, 
			 I=unique(Ms), V=V_, R=R_, qweights = qweights_, pweights=pweights_, 
			 qvote = qvote_, pvote = pvote_,
			 betahat=t(as.matrix(betahat)) , verbose=FALSE )
		logprior <- log(dmnorm(betahat, mean=b0, varcov=solve(invB0) )) 
		logplugin <- apply(bmuchains, 1, function(x) dmnorm(betahat, mean=x, varcov=bsd, log=FALSE ) )

		logprior <- as.numeric(logprior)
		logplugin <- as.numeric(logplugin)
		loglik <- as.numeric(loglik)

		attr(betachains, "logmarglik") <- logprior + loglik  - log(mean(logplugin))
		attr(betachains, "logprior") <- logprior
		attr(betachains, "loglik") <- loglik
		attr(betachains, "logplugin") <- -log(mean(logplugin))
	}


	return(betachains)

	}
