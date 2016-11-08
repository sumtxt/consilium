// Original code by Chris Hans
// Distributed in the R-package 'truncatedNormals', Version 0.5 (27.5.2013) in "trunc_norm.c" 
// Source: http://www.stat.osu.edu/~pfc/software/truncatedNormals/
// Licence:  Attribution-NonCommercial-ShareAlike 4.0, International license.
//
// Changes by Moritz Marbach: 
//    - explicit scope qualifier ("R::") for all distribution functions
//    - replaced norm_rand() with rnorm(0.0,0.1)
//    - droppped rnorm_truncated() - not needed 


// ======================================================================
// norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to be in the interval
// (a,b) via rejection sampling.
// This function should be called by rnorm_truncated (where Get/PutRNGstate
// are invoked)
// ======================================================================
double
norm_rs(double a, double b)
{
   double   x;
   x = R::rnorm(0.0, 1.0);
   while( (x < a) || (x > b) ) x = R::rnorm(0.0, 1.0);
   return x;
}

// ======================================================================
// half_norm_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) (with a > 0) using half normal rejection sampling.
// This function should be called by rnorm_truncated (where Get/PutRNGstate
// are invoked)
// ======================================================================
double
half_norm_rs(double a, double b)
{
   double   x;

   //assert(a >= 0); // check it

   x = fabs(R::rnorm(0.0, 1.0));
   while( (x<a) || (x>b) ) x = fabs(R::rnorm(0.0, 1.0));
   return x;
}

// ======================================================================
// unif_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using uniform rejection sampling. 
// This function should be called by rnorm_truncated (where Get/PutRNGstate
// are invoked)
// ======================================================================
double
unif_rs(double a, double b)
{
   double xstar, logphixstar, x, logu;

   // Find the argmax (b is always >= 0)
   // This works because we want to sample from N(0,1)
   if(a <= 0.0) xstar = 0.0;
   else xstar = a;
   logphixstar = R::dnorm(xstar, 0.0, 1.0, 1);

   x = R::runif(a, b);
   logu = log(R::runif(0.0, 1.0));
   while( logu > (R::dnorm(x, 0.0, 1.0, 1) - logphixstar))
   {
      x = R::runif(a, b);
      logu = log(R::runif(0.0, 1.0));
   }
   return x;
}

// ======================================================================
// exp_rs(a, b)
// generates a sample from a N(0,1) RV restricted to the interval
// (a,b) using exponential rejection sampling.
// This function should be called by rnorm_truncated (where Get/PutRNGstate
// are invoked)
// ======================================================================
double
exp_rs(double a, double b)
{
  double z, u, scale;

  scale = 1.0/a;

   // Generate a proposal on (0, b-a)
   z = R::rexp(scale);
   while(z > (b-a)) z = R::rexp(scale);
   u = R::runif(0.0, 1.0);

   while( log(u) > (-0.5*z*z))
   {
      z = R::rexp(scale);
      while(z > (b-a)) z = R::rexp(scale);
      u = R::runif(0.0,1.0);
   }

   return(z+a);
}



//======================================================================
// rnorm_trunc( mu, sigma, lower, upper)
//
// generates one random normal RVs with mean 'mu' and standard
// deviation 'sigma', truncated to the interval (lower,upper), where
// lower can be -Inf and upper can be Inf.
//======================================================================
double 
rnorm_trunc (double mu, double sigma, double lower, double upper)
{
 int  change;
 double  a, b;
 double  logt1 = log(0.150), logt2 = log(2.18), t3 = 0.725, t4 = 0.45;
 double  z, tmp, lograt;

 change = 0;
 a = (lower - mu)/sigma;
 b = (upper - mu)/sigma;
 
 // First scenario
 if( (a == R_NegInf) || (b == R_PosInf))
   {
     if(a == R_NegInf)
       {
    change = 1;
    a = -b;
    b = R_PosInf;
       }
     
     // The two possibilities for this scenario
     if(a <= 0.45) z = norm_rs(a, b);
     else z = exp_rs(a, b);
     if(change) z = -z;
   }
 // Second scenario
 else if((a * b) <= 0.0)
   {
     // The two possibilities for this scenario
     if((R::dnorm(a, 0.0, 1.0, 1) <= logt1) || (R::dnorm(b, 0.0, 1.0, 1) <= logt1))
       {
    z = norm_rs(a, b);
       }
     else z = unif_rs(a,b);
   }
 // Third scenario
 else
   {
     if(b < 0)
       {
    tmp = b; b = -a; a = -tmp; change = 1;
       }
     
     lograt = R::dnorm(a, 0.0, 1.0, 1) - R::dnorm(b, 0.0, 1.0, 1);
     if(lograt <= logt2) z = unif_rs(a,b);
     else if((lograt > logt1) && (a < t3)) z = half_norm_rs(a,b);
     else z = exp_rs(a,b);
     if(change) z = -z;
   }
 
 return (sigma*z + mu);
}
