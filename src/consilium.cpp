#include <RcppArmadillo.h>

using namespace Rcpp ;
// [[Rcpp::depends(RcppArmadillo)]]

#include "rnorm_trunc.h"


NumericVector rbern_v2(NumericVector x, double Q){
	NumericVector draws(x.size());
	for (int i = 0; i < x.size(); i++){
		double u = runif(1,0,1)(0);
		if ( u < (Q * x(i)) ) { 
			draws(i) = 1; 
		} else { 
			draws(i) = 0;
		}
		}
	return draws;
	}


NumericVector rbern_v(NumericVector x){
	NumericVector draws(x.size());
	for (int i = 0; i < x.size(); i++){
		draws(i) = rbinom(1,1,x(i) )(0);
		}
	return draws;
	}


int sum_veto(NumericVector x, int V){
	NumericVector vetos(V);
	for(int i=0; i < V; i++) {
		if (x(i) == 0) vetos(i) += 1;
		}
	return sum(vetos);
	}


int sum_aye(NumericVector x, int V, int I){
	NumericVector ayes(I-V);
	for(int i=0; i < (I-V); i++){
		if (x(V+i) == 1) ayes(i) += 1;
		}
	return sum(ayes);
	}


int sum_weights(NumericVector x, IntegerVector weights, int I){
	NumericVector ayes(I);
	for(int i=0; i < I; i++){
		if ( x(i) == 1 ) ayes(i) += weights(i);
		}
	return(sum(ayes));
	}

// [[Rcpp::export]]
int check_vprofile(NumericVector ydraw, int V, int R, int I, int qvote, int pvote, 
		IntegerVector qweights, IntegerVector pweights){
	int totalvetos = sum_veto(ydraw,V);
	int totalayes = sum_aye(ydraw,V,I);
	int totalqvote = 0;
	int totalpvote = 0;
	if ( qvote != 0 ) totalqvote = sum_weights(ydraw,qweights, I);
	if ( pvote != 0 ) totalpvote = sum_weights(ydraw,pweights, I);
	if ( (totalvetos == 0 & totalayes >= R & totalqvote >= qvote & totalpvote >= pvote) ) {
		return(1);
	} else {
		return(0);
	}
	}


//' Random Draws of Bernoulli Variables Subject to Global Constraints 
//'
//' \code{rbern_constr_v} draws a vector of Bernoulli variables subject to global constraints such as 
//' the number of 1's in each draw.
//'
//'  
//'
//' @param x (required) the vector of probabilities. 
//' @param R (required) how many 1's does the the final draw have to have? 
//' @param V (required) how many of the first \code{V} entries of a draw have the be 1?
//' @param z (required) if \code{z==1} all constrains are obeyed, if \code{z==0} all constraints are not.
//' @param qweights (required) integer numbers. The first weights attached to a draw. 
//' @param pweights (required) integer numbers. The second weights attached to a draw. 
//' @param qvote (required) the smallest sum of \code{qweights} for all entries that are non-zero 
//'		in a draw before the draw is accepted. If \code{qvote=0} the \code{pweights} are ignored. 	
//' @param pvote (required) the smallest sum of \code{pweights} for all entries that are non-zero 
//'		in a draw before the draw is accepted. If \code{qvote=pvote} the \code{pweights} are ignored. 	
//' @param adapt (required) 0 deactivates adaption algorithm 2 in Marbach (2016), 1 enables it. 
//' @param rate (required) schedule for the adaption with algorithm 2 
//' @param step (required) \eqn{\epsilon}-parameter in algorithm 2 
//' @param verbose (required) use \code{0} for no output; use any positive natural number to report progress at each \code{verbose}-th iteration.
//'
//' @details 
//' The function is part of the Gibbs sampler used for \code{BayesPMP()}. The function runs in \code{C++}.
//' The function draws a random vector of Bernoulli variables and checks the global constraints. 
//' If the draw does not obey the constraints  a new draw is obtained until the constraints are met. 
//' 
//' @return \code{vector} of draws. 
//'  
//' @examples 
//'  \dontrun{
//'  # Can be used to simulate vote choices of 3-member committee that has passed a proposal with majority rule
//'  
//'  # Vote choice probabilities 
//'  p <- c(0.2,0.3,0.4)
//'  R <- 2
//'  V <- 0 
//'  z <- 1 
//'  pweights <- qweights <- rep(1,3)
//'  qvote <- pvote <- 0
//'  adapt <- 1
//'  rate <- 200 
//'  step <- 0.05
//'  verbose <- 0
//'  
//'  rbern_constr_v(x=p, R=R, V=V, z=z, qweights=qweights, pweights=pweights, qvote=qvote, 
//'	   pvote=pvote, adapt=adapt, rate=rate, step=step, verbose=verbose)
//'  }
//'  
//'
//' @export
// [[Rcpp::export]]
NumericVector rbern_constr_v(NumericVector x, int R, int V, int z, 
	IntegerVector qweights, IntegerVector pweights, int qvote, int pvote, 
	int adapt, int rate, double step, int verbose){

	// Rprintf("R: %i. qvote: %i. pvote: %i, \n", R, qvote, pvote ); R_FlushConsole(); R_ProcessEvents();
	// Rprintf("qweights[1]: %i. pweights[1]: %i. \n", qweights(0), pweights(0) ); R_FlushConsole(); R_ProcessEvents();

	// RNGScope scope;
	int I = x.size();
	if (I < R) throw exception("R > I");
	if (I < V) throw exception("V > I");
	if (I < V+R) throw exception(" V+R => I ");
	NumericVector ydraw = rbern_v2(x,1);
	int totalvetos = sum_veto(ydraw,V);
	int totalayes = sum_aye(ydraw,V,I);
	int totalqvote = 0;
	int totalpvote = 0;
	if (qvote != 0 ) totalqvote = sum_weights(ydraw,qweights, I); 
	if ( pvote != 0 ) totalpvote = sum_weights(ydraw,pweights, I); 

	double Q = 1;
	double Qmax = 1;
	int count = 1;

	if ( z==1 ){

		if (adapt ==1 ) Qmax = 1/max(x);	
		// sample ydraws for a accepted proposal
		while( !(totalvetos == 0 & totalayes >= R & totalqvote >= qvote & totalpvote >= pvote) ) {
			ydraw = rbern_v2(x,Q);
			totalvetos = sum_veto(ydraw,V);
			totalayes = sum_aye(ydraw,V,I);
			if (qvote != 0 ) totalqvote = sum_weights(ydraw,qweights, I); 
			if ( pvote != 0 ) totalpvote = sum_weights(ydraw,pweights, I); 
			count += 1;
			if ( ((count % rate)==0) & (Q < Qmax) ) Q += step;
			}
	} else { 

		if (adapt ==1 ) Qmax = 1/max(1-x);	
		// sample ydraws for a rejected proposal
		while( totalvetos == 0 & totalayes >= R & totalqvote >= qvote & totalpvote >= pvote ) {
			ydraw = 1-rbern_v2(1-x,Q);
			totalvetos = sum_veto(ydraw,V);
			totalayes = sum_aye(ydraw,V,I);
			if (qvote != 0 ) totalqvote = sum_weights(ydraw,qweights, I); 
			if ( pvote != 0 ) totalpvote = sum_weights(ydraw,pweights, I); 
			count += 1;
			if ( ((count % rate)==0) & (Q < Qmax) ) Q += step;
			}
	}
	if ( (verbose != 0) & (Qmax != Q) ) { Rprintf("Iterations: %i. Q/Qmax: %5.2f/%5.2f \n", count, Q, Qmax ); R_FlushConsole(); R_ProcessEvents(); }			
	return ydraw;
	}


// GIBBS SAMPLING FUNCTIONS //
/////////////////////////////

NumericVector drawAlpha(double omega, IntegerVector nJ, 
                              NumericVector ystarbar, int J){
  NumericVector draws(J);
  double sigma2; double mu;
  for(int j=0; j < J; j++) {
    sigma2 = 1/( omega + nJ(j));
    mu = sigma2 * (ystarbar(j) * nJ(j));
    draws(j) = rnorm(1, mu , sqrt(sigma2))(0);
    }
    return draws;
    }

NumericVector calcYstarBar(arma::colvec e, IntegerVector alpha_idx, 
  IntegerVector nJ, int J, int N){
  NumericVector ystarbar(J,0.0);
  for(int n=0; n < N; n++) {
    ystarbar(alpha_idx(n)) += e(n)/nJ(alpha_idx(n));
  } 
  return(ystarbar);
  }

double drawOmega(double e0,double f0,int J, NumericVector alpha){
  double e1; double f1; double omega;
  e1 = e0 + (J/2);
  f1 = f0 + (sum( pow(alpha,2))/2);
  // Rcpp rgamma has a different param than rgamma in R!
  omega = rgamma(1,e1,1/f1)(0);
  return omega;
  }

arma::colvec makeAlpha(NumericVector alphadraw, 
          IntegerVector alpha_idx, int N){
  arma::colvec alpha(N);
  alpha.fill(arma::datum::nan); 
  for(int n=0; n < N; n++) {
    alpha(n) = alphadraw(alpha_idx(n));
    }
  return alpha;
  }


arma::colvec mvrnormArma(arma::colvec mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   NumericVector YY = rnorm(ncols, 0 , 1);
   arma::colvec Y(YY.begin(), YY.size(), false);
   return mu + arma::chol(sigma).t() * Y ;
   }


//' @export
// [[Rcpp::export]]
List consiliumcA(arma::mat X, NumericVector y, IntegerVector z, 
		 int K, IntegerVector Is, int J, IntegerVector Rs, IntegerVector Vs, 
		 IntegerMatrix qweights, IntegerMatrix pweights,
		 IntegerVector qvote, IntegerVector pvote, IntegerVector pslct, IntegerVector qslct,
		 IntegerVector miss, int NGIBBS, int NTHIN, int NBURN, int VERBOSE, 
		 int adapt, int rate, double step, int VERBOSE2,
		 arma::mat invB0, arma::colvec b0, arma::colvec betastart, 
		 double omegastart, int G, IntegerVector alpha_idx, IntegerVector nG, 
		 double e0, double f0) {

	// CONTROLLS
	RNGScope scope;
	int NSAMP=(NGIBBS/NTHIN);
	int NP = 1;

	// Compute constants 
	arma::mat invB0tXX = inv( invB0 + ( X.t() * X ) );
	
	// Initialize variables
	double lower; double upper;
	
	int N = X.n_rows;
	arma::colvec ystar(N);
	arma::colvec b1;

	arma::mat beta(K, (NSAMP+1));
	beta.fill(arma::datum::nan); 
	beta.col(0) = betastart; 
	
	arma::mat bmu(K, (NSAMP+1));
	bmu.fill(arma::datum::nan); 
	bmu.col(0) = betastart; 
	
	arma::colvec betadraw = betastart;

	double omegadraw = omegastart;  

	arma::colvec omega(NSAMP+1);
	omega.fill(arma::datum::nan); 
	omega(0) = 1/omegastart; 

	NumericMatrix alphaout(G, (NSAMP+1));

	NumericVector ystarbar(G);
 	NumericVector alphadraw(G);
	arma::colvec e(N);

	// Fill with 0 otherwise if G==0 the calculations are wrong!
	arma::colvec alpha(N);
	alpha.fill(0.0);

	// RUN GIBBS
	for(int s=0; s < (NGIBBS + NBURN); s++){

		// Clean storage matrices / variables
	    	ystar.fill(arma::datum::nan);
		int idx = 0; 


		for(int j=0; j < J; j++){ 

			int I = Is(j);

			NumericVector yvec(I);
			NumericVector ystarmu(I);	

			ystarmu.fill(0); 

			// >> Calcualte linear predictor <<
			// Important: ystarmu has to be initialized at 0!
			for (int i=0; i < I; i++) { 
	        	for(int k=0; k < K; k++){
		        	ystarmu(i) += betadraw(k) * X(idx+i,k); 
			       }
			    ystarmu(i) += alpha(idx+i);
			    }

	    		// >> Sample vote profile or obtain from data <<
			if ( miss(j) == 1 ) {
				yvec = rbern_constr_v(pnorm(ystarmu), Rs(j), Vs(j), z(idx), 
					qweights(_,qslct(j)), pweights(_,pslct(j)), qvote(j), pvote(j), 
					adapt, rate, step, VERBOSE2 ) ; 
			} else {
			for (int i=0; i < I; i++) { 
					yvec(i) = y(idx+i);
					}
				}

			// >> Sample latent utility << 	
		    	for(int i=0; i<I; i++){
		    		if ( yvec[i] == 0 ) {  lower = R_NegInf;  } else {  lower = 0; }
			    	if ( yvec[i] == 0 ) {  upper = 0;  } else {  upper = R_PosInf; }
			    	ystar(idx+i) = rnorm_trunc( ystarmu(i), 1, lower, upper) ;
			}

			idx += I;

	  	   	}

	    // >> Sample the coefficients <<
	    b1 = invB0tXX * ( (invB0 * b0) + (X.t() * (ystar-alpha)) );
	  	betadraw = mvrnormArma(b1, invB0tXX);

	    // >> Sample alpha <<< 
	    if (G != 0){
		    e = ystar - (X * betadraw);
		    ystarbar = calcYstarBar(e, alpha_idx, nG, G, N);
		    alphadraw = drawAlpha(omegadraw, nG, ystarbar, G);
		    alpha = makeAlpha(alphadraw,alpha_idx,N);
		    omegadraw = drawOmega(e0,f0,G,alphadraw);
		    }

	    // Flush to console
	    if ( (s >= NBURN) && ((s % NTHIN)==0) ) { 
	    	beta.col(NP) = betadraw; 
	    	bmu.col(NP) = b1; 
	    	if (G != 0){
		    	omega(NP) = omegadraw;
		    	alphaout(_,NP) = alphadraw;
		    	}
	    	NP++; 
	    }
	    if (VERBOSE != 0) if ( ((s % VERBOSE) == 0) ) { Rprintf("%i/%i - ", s, NGIBBS ); R_FlushConsole(); R_ProcessEvents(); }
	    }
	
	if (G == 0){
		return List::create(Named("beta") = beta, 
					  Named("bmu") = bmu );
	} else { 
	return List::create(Named("beta") = beta, 
				  Named("bmu") = bmu, 
				  Named("omega") = omega,
				  Named("alpha") = alphaout );
		}

	}




//' @export
// [[Rcpp::export]]
arma::mat consiliumppc(arma::mat posterior, arma::mat X,  
		 IntegerVector Rs, IntegerVector Vs, IntegerVector Is, 
		 IntegerVector qvote, IntegerVector pvote, IntegerMatrix qweights, 
		 IntegerMatrix pweights, IntegerVector pslct, IntegerVector qslct,
		 int J, int K, int NSAMP, 
		 int SIM, int VERBOSE) {


	// CONTROLLS
	RNGScope scope;

	// Holding-bins
	arma::mat pp(J, NSAMP);
	pp.fill(arma::datum::nan);
	
	NumericVector z(SIM);

	// RUN Simulation
	for(int s=0; s < NSAMP; s++){

		int idx = 0; 
		for(int j=0; j < J; j++){

			// Clean storage matrices / variables
			int I = Is(j);

			NumericVector yvec(I);
			NumericVector ystarmu(I);	
			ystarmu.fill(0); 

			// >> Calcualte linear predictor <<
			// Important: ystarmu has to be initialized at 0!
			for (int i=0; i < I; i++) { 
	        	for(int k=0; k < K; k++){
		        	ystarmu(i) += posterior(s,k) * X(idx+i,k); 
			        }
			    }

			z.fill(arma::datum::nan);
			for(int p=0; p < SIM; p++){    
				    yvec = rbern_v(pnorm(ystarmu)) ; 
				    if ( check_vprofile( yvec, Vs(j),Rs(j), I, qvote(j), pvote(j), 
				    		qweights(_,qslct(j)), pweights(_,pslct(j)) )==1 ) { 
				    	z(p) = 1; 
				    } else { 
				    	z(p) = 0; 
				    }
				    }	
				pp(j,s) = mean(z);
			
			idx += I;
	  	   	}

	  	if (VERBOSE != 0) if ( ((s % VERBOSE) == 0) ) { Rprintf("%i/%i - ", s, NSAMP ); R_FlushConsole(); R_ProcessEvents(); }

	}
	
	return pp;

	}	
