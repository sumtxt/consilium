#' Computes All Vote Profiles Given a Veto-Majority Voting Rule 
#'
#' \code{VoteProfiles} computes matrices of vote profiles given a voting rule. 
#'
#' @param R (required) the number of members that have to agree 
#' 		to pass a proposal excluding the members that have veto power. 
#' @param V the number of members that have to agree to pass a proposal. 
#' @param I (required) the number of members voting. 
#'
#' @return A \code{list} object with two matrices that list all vote profiles 
#' 	which imply the adoption (\code{yescols}) or rejection (\code{nocols}) of a proposal.  
#'  A column refers to a voting member and a row to a vote profile. 
#' 
#' Notice that this function is time- and memory intensive if \code{I} is large. 
#' 
#' 
#' @export
VoteProfiles <- function(R,V=0,I){
	cols <- mapply(function(x) t(combn(seq(1,I),x, fun=constr,a=I)), seq(1,I), SIMPLIFY=FALSE)
	cols <- do.call(rbind, cols)
	st <- vetomajority_rule(cols, R=R, V=V, I=I)
	yescols <- matrix(cols[st,], ncol=I)
	nocols <- matrix(cols[!st,], ncol=I)
	nocols <- rbind(nocols, rep(FALSE, I))
	out <- list( yescols=yescols, nocols=nocols)
	return(out)
	}



