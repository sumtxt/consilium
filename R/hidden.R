is.naturalnumber <- function(x, tol = .Machine$double.eps^0.5)  x > tol & abs(x - round(x)) < tol

constr <- function(x,a){
	v <- rep(0,a)
	v[x] <- 1
	return(as.logical(v))
	}

#' @export
setClassUnion("numericOrNULL", c("numeric","NULL"))

#' @export
setClass("maxlikpmp", slots = c(
	loglik = "numeric",
	n = "numeric",  
	coef = "numeric", 
	coefnames = "character",
	hessian = "matrix", 
	form = "formula",
	code = "character", 
	iter = "numeric", 
	lik = "character",
	I = "numericOrNULL", 
	R = "numeric", 
	V = "numericOrNULL",
	qweights = "numericOrNULL",
	pweights = "numericOrNULL", 
	qvote = "numericOrNULL", 
	pvote = "numericOrNULL"
	))

#' @export
setMethod("vcov", "maxlikpmp", function(object) object@hessian )

#' @export
setMethod("coef", "maxlikpmp", function(object) object@coef )

#' @export
setMethod("logLik", "maxlikpmp", function(object) object@loglik )

#' @export
setMethod("summary", "maxlikpmp", function(object) { 
	require(texreg)
	lik <- logLik(object)
	n <- object@n
	k <- length(coef(object))
	aic <-  2*k - 2*(lik)
	bic <-  log(n)*k - 2*(lik)
	gof <- c(aic, bic, lik, n)
	gof.names <- c("AIC", "BIC", "Log Likelihood", "Num. obs.")
	decimal.places <- c(TRUE, TRUE, TRUE, FALSE) 
	ses <- sqrt(diag(vcov(object)))
	zvalue <- coef(object)/ses
	pvalue <-  2*pnorm(-abs(zvalue))
	res <- createTexreg(
		coef.names = object@coefnames,
		coef = coef(object),
		se = ses,
		pvalues = pvalue,
		gof.names = gof.names,
		gof  = gof,
		gof.decimal = decimal.places 
		)
	screenreg(res)
	})


# These two functions are necessary to calculate the log-likelihood based on 
# an exact set of possible coalitions (the V^+/V^- sets)
logliks <- function(data, beta0, coalitions){
	out <- sapply(data,  function(x) {
		pro <- pnorm(x %*% beta0) 
		neg <- 1-pro
		probs <- apply(coalitions, 1, function(coalition) {
			prod(pro[coalition==TRUE]) * prod(neg[coalition==FALSE])
			})
		return(log(sum(probs)) )
		})
	return(out)
	}

comblikpmp <- function(beta0, X_1, X_0, nocols, yescols,verbose){
	
	logliks_0 <- logliks(X_0,beta0,nocols)
	logliks_1 <- logliks(X_1,beta0,yescols)

	logl <- sum(c(logliks_0, logliks_1))

	if( verbose==TRUE ) {
		cat(paste(" Log: ", as.character(round(logl,3)),  
			"\tBeta: ", paste(as.character(round(beta0,3)), sep=" ", collapse=" | "),"\n"))
		}
	return(-logl)
	}


poibinprob <- function(x,R,V,I) {
    if (V == 0) return( 1-ppoibin( (R-1), x ) )
    if (V == 1) return( x[V] * (1-ppoibin( (R-1), x[2:I] ))  )
    if (V > 1 & V < I) return( prod(x[1:V]) * (1-ppoibin( (R-1), x[(V+1):I])) )
    if ( V == I) return( prod(x) )
    }

poibinlikpmp <- function(beta0, X_1, X_0, R,V,I, verbose=TRUE){

	p1 <- sapply(X_1, function(x) {
		pp <- pnorm(x %*% beta0) 
		return( poibinprob(pp,R,V,I) )
		} )
	
	p0 <- sapply(X_0, function(x) {
			pp <- pnorm(x %*% beta0) 
			return( 1-poibinprob(pp,R,V,I)  )
			} )

	logl <- sum(log(p0)) + sum(log(p1))

	if( verbose==TRUE ) {
		cat(paste(" Log: ", as.character(round(logl,3)),  
			"\tBeta: ", paste(as.character(round(beta0,3)), sep=" ", collapse=" | "),"\n"))
		}
	return(-logl)
	}


# This log-likelihood can be used if there is not within-decision variability 
bccdf <- function(x,I,R,V) pbinom(R-1, I, pnorm(x), lower.tail = FALSE) * pnorm(x)^V

blinklikpmp <- function(beta,y,X,I,R,V){
	zeta <- as.vector(X %*% beta)
	y.logic <- as.logical(y)
	n<-length(y)
	lgLik <- rep(NA,n)
	lgLik[y.logic] <- log(bccdf(zeta[y.logic],I,R,V))
	lgLik[!y.logic] <- log(1-bccdf(zeta[!y.logic],I,R,V))
	logl<- sum(lgLik)
	return(-logl)
	}


# Voting rules
# vetomajority_rule <- function(x,R,V,I) {
#     if (V == 0) return( (rowSums(x) >= R) )
#     if (V == 1) return( ( x[,1] == TRUE) & (rowSums(x[,2:I]) >= R) )
#     if (V > 1 & V < I) return( ( rowSums(x[,1:V]) == V ) & (rowSums(x[,(V+1):I]) >= R))
#     if ( V == I) return( ( rowSums(x) == V ) )
#     }
