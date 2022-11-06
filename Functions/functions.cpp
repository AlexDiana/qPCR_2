#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// // [[Rcpp::depends(RcppParallel)]]
// #include <RcppParallel.h>
  // 
  // using namespace RcppParallel;

#include <cmath>
#include <algorithm>

const double TRUNC = .64;
const double TRUNC_RECIP = 1.0 / .64;

const double log2pi = std::log(2.0 * M_PI);

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392


// [[Rcpp::export]]
double dt2(double x, double mean, double scale, double df){
  double tstat = (x - mean) / scale;
  return(R::dt(tstat, df, 0) / scale);
}

// [[Rcpp::export]]
double rt2(double mean, double scale, double df){
  double tstat = R::rt(df);
  return(tstat * scale + mean);
}

// sample from normal distribution


arma::vec mvrnormArma(arma::vec mu, arma::mat Sigma) {
  int ncols = Sigma.n_cols;
  arma::mat Y = arma::randn(1, ncols);
  return mu + arma::trans(Y * arma::chol(Sigma));
}

// [[Rcpp::export]]
arma::vec mrt2(arma::vec mean, arma::mat Sigma, double df){
  arma::vec zeroVec = arma::zeros(mean.size());
  arma::vec y = mvrnormArma(zeroVec, Sigma);
  double u = R::rchisq(df);
  arma::vec x = sqrt(df / u) * y + mean;
  return x;
}

// [[Rcpp::export]]
double dmt_cpp(arma::vec x, double nu, arma::vec mu, arma::mat Sigma, bool returnLog){
  int p = x.size();
  
  double logratio = R::lgammafn(( nu + p ) / 2) - R::lgammafn( nu / 2 );
  
  arma::vec product = (arma::trans(x - mu) * arma::inv(Sigma) * (x - mu));
  double lognum = (- ( nu + p ) / 2) * log(1 + (1 / nu) * product[0]);
  double logden = (p / 2.0) * log(M_PI * nu) + (0.5) * log(arma::det(Sigma));
  
  double loglikelihood = logratio + lognum - logden;
  
  if(returnLog){
    return loglikelihood;
  } else {
    return exp(loglikelihood);
  }
}

// GAUSSIAN PROCESS FUNCTIONS

double k_cpp(double x1, double x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-(x1-x2)*(x1-x2)/(2*pow(l,2)));
  // return 1;
}

arma::mat K(arma::vec x1, arma::vec x2, double a, double l){
  arma::mat res(x1.size(), x2.size());
  
  for(int i = 0; (unsigned)i < x1.size(); i++){
    for(int j = 0; (unsigned)j < x2.size(); j++){
      res(i,j) = k_cpp(x1[i],x2[j], a, l);
    }  
  }
  
  return res;
}

double k2_cpp(arma::rowvec x1, arma::rowvec x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-( pow(x1[0]-x2[0], 2) + pow(x1[1]-x2[1], 2) ) /(2*pow(l,2)));
}

// [[Rcpp::export]]
arma::mat K2(arma::mat x1, arma::mat x2, double a, double l){
  arma::mat res(x1.n_rows, x2.n_rows);
  
  for(int i = 0; (unsigned)i < x1.n_rows; i++){
    for(int j = 0; (unsigned)j < x2.n_rows; j++){
      res(i,j) = k2_cpp(x1.row(i),x2.row(j), a, l);
    }  
  }
  
  return res;
}

//
  
double aterm(int n, double x, double t) {
    double f = 0;
    if(x <= t) {
      f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
    }
    else {
      f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
    }    
    return (double)exp(f);
  }

double exprnd(double mu) {
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

double truncgamma() {
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

double randinvg(double mu) {
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

double tinvgauss(double z, double t) {
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }    
  return X;
}

double samplepg(double z) {
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
double t = MATH_2_PI;

// Compute p, q and the ratio q / (q + p)
// (derived from scratch; derivation is not in the original paper)
double K = z*z/2.0 + MATH_PI2/8.0;
double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
double logK = (double)std::log(K);
double Kt = K * t;
double w = (double)std::sqrt(MATH_PI_2);

double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
double ratio = 1.0 / (1.0 + p_over_q); 

double u, X;

// Main sampling loop; page 130 of the Windle PhD thesis
while(1) 
{
  // Step 1: Sample X ? g(x|z)
  u = R::runif(0.0,1.0);
  if(u < ratio) {
    // truncated exponential
    X = t + exprnd(1.0)/K;
  }
  else {
    // truncated Inverse Gaussian
    X = tinvgauss(z, t);
  }
  
  // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
  int i = 1;
  double Sn = aterm(0, X, t);
  double U = R::runif(0.0,1.0) * Sn;
  int asgn = -1;
  bool even = false;
  
  while(1) 
  {
    Sn = Sn + asgn * aterm(i, X, t);
    
    // Accept if n is odd
    if(!even && (U <= Sn)) {
      X = X * 0.25;
      return X;
    }
    
    // Return to step 1 if n is even
    if(even && (U > Sn)) {
      break;
    }
    
    even = !even;
    asgn = -asgn;
    i++;
  }
}
return X;
}

// [[Rcpp::export]]
double rpg(int n, double z){
  
  double x = 0;
  for(int i = 0; i < n; i++){
    x += samplepg(z);
  }
  
  return(x);
}

// LAMBERT

const double EPS = 2.2204460492503131e-16;
const double M_1_E = 1.0 / M_E;

double FritschIter(double x, double w_guess){
  double w = w_guess;
  int MaxEval = 5;
  bool CONVERGED = false;
  double k = 2.0 / 3.0;
  int i = 0;
  do {
    double z = std::log(x / w) - w;
    double w1 = w + 1.0;
    double q = 2.0 * w1 * (w1 + k * z);
    double qmz = q - z;
    double e = z / w1 * qmz / (qmz - z);
    CONVERGED = std::abs(e) <= EPS;
    w *= (1.0 + e);
    ++i;
  } while (!CONVERGED && i < MaxEval);
  return(w);
}

// [[Rcpp::export]]
double lambertW0_CS(double x) {
  double result;
  double w;
  if (x == R_PosInf) {
    result = R_PosInf;
  } else if (x < -M_1_E) {
    result = R_NaN;
  } else if (std::abs(x + M_1_E) <= EPS) {
    result = -1.0;
  } else if (x <= M_E - 0.5) {
    if (std::abs(x) <= 1e-16) {
      /* This close to 0 the W_0 branch is best estimated by its Taylor/Pade
      expansion whose first term is the value x and remaining terms are below
      machine double precision. See
      https://math.stackexchange.com/questions/1700919/how-to-derive-the-lambert-w-function-series-expansion
      */
        result = x;
    } else {
      if (std::abs(x) <= 7e-3) {
        /* Use equation (5) in Fritsch */
          w = ((1.33333333333333333 * x + 1.0) * x) /
            ((0.83333333333333333 * x + 2.33333333333333333) * x + 1.0);
      } else {
        /* Use expansion in Corliss 4.22 to create (3, 2) Pade approximant
        Numerator:-10189 / 303840 * p^3 + 40529 / 303840 * p^2 + 489 / 844 * p-1
        Denominator: -14009 / 303840 * p^2 + 355 / 844 * p + 1
        Converted to digits to reduce needed operations
        */
          double p = std::sqrt(2.0 * (M_E * x + 1.0));
        double Numer = ((-0.03353409689310163 * p + 0.1333892838335966) * p +
                          0.5793838862559242) * p - 1.0;
        double Denom = (-0.04610650342285413 * p + 0.4206161137440758) * p + 1.0;
        w = Numer / Denom;
      }
      result = FritschIter(x, w);
    }
  } else {
    /* Use first five terms of Corliss et al. 4.19 */
      w = std::log(x);
      double L_2 = std::log(w);
      double L_3 = L_2 / w;
      double L_3_sq = L_3 * L_3;
      w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
        w + L_3_sq * L_3 / 3.0;
      result = FritschIter(x, w);
  }
  return(result);
}

// [[Rcpp::export]]
double lambertWm1_CS(double x){
  double result;
  double w;
  if (x == 0.0) {
    result = R_NegInf;
  } else if (x < -M_1_E || x > 0.0) {
    result = R_NaN;
  } else if (std::abs(x + M_1_E) <= EPS) {
    result = -1.0;
  } else {
    /* Use first five terms of Corliss et al. 4.19 */
      w = std::log(-x);
      double L_2 = std::log(-w);
      double L_3 = L_2 / w;
      double L_3_sq = L_3 * L_3;
      w += -L_2 + L_3 + 0.5 * L_3_sq - L_3 / w + L_3 / (w * w) - 1.5 * L_3_sq /
        w + L_3_sq * L_3 / 3.0;
      result = FritschIter(x, w);
  }
  return(result);
}

// compute product between a matrix and a  diagonal matrix (summarised in a vector) 

arma::mat diagMatrixProd(arma::mat& X, arma::vec& D){
  
  arma::mat result(X.n_rows, D.size());
  for(int i = 0; i < result.n_rows; i++){
    for(int j = 0; j < result.n_cols; j++){
      result(i, j) = X(i,j) * D(j);
    }  
  }
  
  return(result);
}

// [[Rcpp::export]]
arma::vec sample_beta_cpp(arma::mat& X, arma::mat& B, arma::vec& b, 
                          arma::vec& Omega, arma::vec& k){
  
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat tXOmegaX = tXOmega * X;
  
  arma::mat invXtOmegaXpB = arma::inv(tXOmegaX + arma::inv(B));
  arma::vec XtkpBb = tX * k + arma::inv(B) * b;
  
  arma::vec result = mvrnormArma(invXtOmegaXpB * XtkpBb, invXtOmegaXpB);
  
  // arma::mat L = arma::trans(arma::chol(tXOmegaX + arma::inv(B))); 
  // arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
  // arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
  // 
    // arma::vec result = mvrnormArma(alpha, arma::trans(arma::inv(arma::trimatl(L))));
    
    return(result);
}

arma::vec sample_Omega_cpp(arma::mat& X, arma::vec& beta, arma::vec& n){
  
  int nsize = n.size();
  arma::vec Omega_vec(nsize);
  
  for(int i = 0; i < nsize; i++){
    
    arma::vec b = X.row(i) * beta;
    
    Omega_vec[i] = rpg(n[i], b[0]);
    
  }
  
  return(Omega_vec);
}

// [[Rcpp::export]]
arma::vec sample_betaPG(arma::vec beta, arma::mat X, arma::vec b,
                        arma::mat B, arma::vec n, arma::vec k){
  
  arma::vec Omega = sample_Omega_cpp(X, beta, n);
  
  beta = sample_beta_cpp(X, B, b, Omega, k);
  
  return(beta);
}

// [[Rcpp::export]]
double dmvnorm_cpp(arma::vec data, arma::vec m, arma::mat Sigma, bool returnLog){
  int xdim = data.size();
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  
  double constants = -(xdim/2) * log2pi;
  arma::vec z = rooti * ( data - m) ;  
  
  if(returnLog){
    return (constants - 0.5 * arma::sum(z%z) + rootisum);
  } else {
    return exp(constants - 0.5 * arma::sum(z%z) + rootisum);     
  }
  
}

// [[Rcpp::export]]
arma::mat rmtrnorm(arma::mat mu, arma::mat U, arma::mat V) {
  
  arma::mat Xmat(mu.n_rows, mu.n_cols);
  for(int i = 0; i < Xmat.n_rows; i++){
    for(int j = 0; j < Xmat.n_cols; j++){
      Xmat(i, j) = R::rnorm(0, 1);
    }  
  }
  
  arma::mat A = arma::chol(U);
  arma::mat B = arma::chol(V);
  
  arma::mat W = A * Xmat * B;
  
  arma::mat Y = W + mu;
  
  return(Y);
}

// [[Rcpp::export]]
arma::mat rmtrnorm_chol(arma::mat mu, arma::mat A, arma::mat B) {
  
  arma::mat Xmat(mu.n_rows, mu.n_cols);
  for(int i = 0; i < Xmat.n_rows; i++){
    for(int j = 0; j < Xmat.n_cols; j++){
      Xmat(i, j) = R::rnorm(0, 1);
    }  
  }
  
  arma::mat W = A * Xmat * B;
  
  arma::mat Y = W + mu;
  
  return(Y);
}

/////////////////////////
  ///////// BETA THETA
/////////////////////////
  
double logistic(double x){
    return(1 / (1 + exp(-x)));
  }

double rtexp(double alpha, 
             double mu_bar){
  return(mu_bar + R::rexp(1 / alpha));
}

double tnorm_std(double mu_bar){
  
  double alpha_star = (mu_bar + sqrt(mu_bar * mu_bar + 4)) / 2;
  
  double z = 0;
  double rho_z = 0;
  do {
    z = rtexp(alpha_star, mu_bar);
    
    rho_z = exp(- (z - alpha_star) * (z - alpha_star) / 2);
    
  } while (!(R::runif(0, 1) < rho_z));
  
  return(z);
  
}

double tnorm(double mu, 
             double sigma, 
             double mu_bar){
  return (mu + sigma * tnorm_std((mu_bar - mu) / sigma));
}

arma::vec truncNormal(arma::vec mu, 
                      arma::mat Sigma, 
                      int truncIndex,
                      double truncation,
                      bool updateBetaTheta1){
  
  int p = mu.size();
  
  // sample marginal of truncated variable
  double x1 = 1;
  if(updateBetaTheta1){
    x1 = tnorm(mu[truncIndex],
               sqrt(Sigma(truncIndex, truncIndex)),
               truncation);
  } 
  
  arma::vec mean1 = arma::zeros(p - 1);
  arma::mat Sigma_11 = arma::zeros(p - 1, p - 1);
  arma::vec Sigma_12 = arma::zeros(p - 1);
  for(int k = 0; k < p; k++){
    if(k < truncIndex){
      mean1[k] = mu[k];
      Sigma_12[k] = Sigma(k, truncIndex);
      for(int k2 = 0; k2 < p; k2++){
        if(k2 < truncIndex){
          Sigma_11(k, k2) = Sigma(k, k2);
        }
        if(k2 > truncIndex){
          Sigma_11(k, k2 - 1) = Sigma(k, k2);
        }
      }
    }
    if(k > truncIndex){
      mean1[k - 1] = mu[k];
      Sigma_12[k - 1] = Sigma(k, truncIndex);
      for(int k2 = 0; k2 < p; k2++){
        if(k2 < truncIndex){
          Sigma_11(k - 1, k2) = Sigma(k, k2);
        }
        if(k2 > truncIndex){
          Sigma_11(k - 1, k2 - 1) = Sigma(k, k2);
        }
      }
    }
  }
  
  // x1 = 2;
  
  arma::vec mu_bar = mean1 + 
    arma::trans(arma::trans(Sigma_12) * (1 / Sigma(truncIndex, truncIndex)) * (x1 - mu[truncIndex]));
  // Rcout << mu_bar << std::endl;
  arma::mat Sigma_bar = Sigma_11 - Sigma_12 * (1 / Sigma(truncIndex, truncIndex)) * arma::trans(Sigma_12);
  
  arma::vec x_left = mvrnormArma(mu_bar, Sigma_bar);
  
  arma::vec out = arma::zeros(p);
  
  out[truncIndex] = x1;
  
  for(int k = 0; k < p; k++){
    if(k < truncIndex){
      out[k] = x_left[k];
    }  
    if(k > truncIndex){
      out[k] = x_left[k - 1];
    } 
  }
  
  return out;
  
}


arma::vec sample_beta_cpp_trunc(arma::mat& X, arma::mat& B, 
                                arma::vec& b, arma::vec& Omega, 
                                arma::vec& k, 
                                bool updateBetaTheta1){
  
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat tXOmegaX = tXOmega * X;
  
  arma::mat invXtOmegaXpB = arma::inv(tXOmegaX + arma::inv(B));
  arma::vec XtkpBb = tX * k + arma::inv(B) * b;
  
  arma::vec result = truncNormal(invXtOmegaXpB * XtkpBb, invXtOmegaXpB, 1, 0, updateBetaTheta1);
  // arma::vec result = truncNormal2(invXtOmegaXpB * XtkpBb, invXtOmegaXpB, 0, 0);
  
  return(result);
}

arma::vec sample_betaPG_trunc(arma::vec beta, arma::mat X, arma::vec b,
                              arma::mat B, arma::vec n, arma::vec k,
                              bool updateBetaTheta1){
  
  arma::vec Omega = sample_Omega_cpp(X, beta, n);
  
  beta = sample_beta_cpp_trunc(X, B, b, Omega, k, 1);
  // beta = sample_beta_cpp_trunc(X, B, b, Omega, k, updateBetaTheta1);
  
  return(beta);
}

arma::vec truncNormal2(arma::vec mu, 
                       arma::mat Sigma, 
                       int truncIndex,
                       double truncation){
  
  int p = mu.size();
  
  // sample marginal of truncated variable
  double x1 = tnorm(mu[truncIndex],
                    sqrt(Sigma(truncIndex, truncIndex)),
                    truncation);
  
  arma::vec mean1 = arma::zeros(p - 1);
  arma::mat Sigma_11 = arma::zeros(p - 1, p - 1);
  arma::vec Sigma_12 = arma::zeros(p - 1);
  for(int k = 0; k < p; k++){
    if(k < truncIndex){
      mean1[k] = mu[k];
      Sigma_12[k] = Sigma(k, truncIndex);
      for(int k2 = 0; k2 < p; k2++){
        if(k2 < truncIndex){
          Sigma_11(k, k2) = Sigma(k, k2);
        }
        if(k2 > truncIndex){
          Sigma_11(k, k2 - 1) = Sigma(k, k2);
        }
      }
    }
    if(k > truncIndex){
      mean1[k - 1] = mu[k];
      Sigma_12[k - 1] = Sigma(k, truncIndex);
      for(int k2 = 0; k2 < p; k2++){
        if(k2 < truncIndex){
          Sigma_11(k - 1, k2) = Sigma(k, k2);
        }
        if(k2 > truncIndex){
          Sigma_11(k - 1, k2 - 1) = Sigma(k, k2);
        }
      }
    }
  }
  
  // x1 = 2;
  
  arma::vec mu_bar = mean1 + 
    arma::trans(arma::trans(Sigma_12) * (1 / Sigma(truncIndex, truncIndex)) * (x1 - mu[truncIndex]));
  // Rcout << mu_bar << std::endl;
  arma::mat Sigma_bar = Sigma_11 - Sigma_12 * (1 / Sigma(truncIndex, truncIndex)) * arma::trans(Sigma_12);
  
  arma::vec x_left = mvrnormArma(mu_bar, Sigma_bar);
  
  arma::vec out = arma::zeros(p);
  
  out[truncIndex] = x1;
  
  for(int k = 0; k < p; k++){
    if(k < truncIndex){
      out[k] = x_left[k];
    }  
    if(k > truncIndex){
      out[k] = x_left[k - 1];
    } 
  }
  
  return out;
  
}

arma::vec sample_beta_cpp_trunc2(arma::mat& X, arma::mat& B, 
                                 arma::vec& b, arma::vec& Omega, 
                                 arma::vec& k){
  
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat tXOmegaX = tXOmega * X;
  
  arma::mat invXtOmegaXpB = arma::inv(tXOmegaX + arma::inv(B));
  arma::vec XtkpBb = tX * k + arma::inv(B) * b;
  
  arma::vec result = truncNormal2(invXtOmegaXpB * XtkpBb, invXtOmegaXpB, 1, 0);
  
  return(result);
}

// [[Rcpp::export]]
arma::vec sample_betaPG_trunc2(arma::vec beta, arma::mat X, arma::vec b,
                               arma::mat B, arma::vec n, arma::vec k){
  
  arma::vec Omega = sample_Omega_cpp(X, beta, n);
  
  beta = sample_beta_cpp_trunc2(X, B, b, Omega, k);
  
  return(beta);
}


// arma::vec logisticXb(arma::mat X, arma::vec beta){
//   
//   arma::vec Xbeta = X * beta;
//   
//   arma::vec plogistic = 1 / (1 + exp(-Xbeta));
//   
//   return(plogistic);
// }

// [[Rcpp::export]]
arma::mat updateAlpha_cpp(arma::vec Ct, 
                      arma::vec w_pcr, 
                      arma::vec delta, 
                      arma::vec P, 
                      arma::vec C_star, 
                      arma::vec w_star, 
                      arma::vec delta_star,
                      arma::vec P_star, 
                      arma::vec sigma_y, 
                      double sigma_alpha, 
                      int n_P, 
                      arma::vec priorAlpha){
    
    arma::vec w_all = join_cols(w_pcr, w_star);
    arma::vec delta_all = join_cols(delta, delta_star);
    arma::vec Ct_all = join_cols(Ct, C_star);
    arma::vec P_all = join_cols(P, P_star);
      
    arma::mat mu_alpha = arma::zeros(n_P, 2);
    arma::cube Lambda_alpha = arma::zeros(n_P, 2, 2);
    
    arma::vec Xt_rowi = arma::zeros(2);
    Xt_rowi[0] = 1;
    
    arma::mat priorMatrix = arma::zeros(2, 2);
    priorMatrix(0, 0) = 1 / (sigma_alpha * sigma_alpha);
    priorMatrix(1, 1) = 1 / (sigma_alpha * sigma_alpha);
    
    int n = w_all.size();
    
    for(int i = 0; i < n; i++){
      
      if(delta_all[i] == 1){
        
        int p_i = P_all[i];
        
        Xt_rowi[1] = log(w_all[i]);
        
        arma::mat Xt_i = Xt_rowi * arma::trans(Xt_rowi);
          
        Lambda_alpha.subcube(arma::span(p_i - 1), arma::span(), arma::span()) += Xt_i;
        
        mu_alpha.row(p_i - 1) += arma::conv_to<arma::rowvec>::from(Xt_rowi * Ct_all[i]);
         
      }
      
    }
    
    arma::mat alpha = arma::zeros(n_P, 2);
    
    for(int p = 0; p < n_P; p++){
      
      double sigma_y_p = sigma_y[p];
      
      arma::mat Lambda_alpha_current = Lambda_alpha.subcube(arma::span(p), arma::span(), arma::span());
      
      Lambda_alpha_current = Lambda_alpha_current / (sigma_y_p * sigma_y_p) + priorMatrix;
      
      arma::vec mu_alpha_current = arma::conv_to<arma::vec>::from(mu_alpha.row(p));
      
      mu_alpha_current = mu_alpha_current / (sigma_y_p * sigma_y_p) + priorAlpha / (sigma_alpha * sigma_alpha);
      
      arma::mat Cov_alpha = arma::inv(Lambda_alpha_current);
      
      arma::vec mean_alpha = Cov_alpha * mu_alpha_current;
      
      arma::vec alpha_j = mvrnormArma(mean_alpha, Cov_alpha);
      alpha.row(p) = arma::conv_to<arma::rowvec>::from(alpha_j);
      
    }
    
    return alpha;
  }


//////
  
  // [[Rcpp::export]]
double rinvgamma_cpp(double a, double b){
  return 1 / R::rgamma(a, 1 / b);
}

arma::mat update_betaz_cpp(arma::mat beta_z, arma::mat l_noempty, arma::vec tau,
                           arma::mat X_beta_noempty, double sigma_beta,
                           arma::vec emptySites, arma::mat B_z){
  
  // X_beta_noempty <- X_beta[emptySites == 0,]
  // l_noempty <- logz[emptySites == 0,]
  int S = tau.size();
  
  int ncov_z = X_beta_noempty.n_cols;
  arma::mat tXX = arma::trans(X_beta_noempty) * X_beta_noempty;
  for(int j = 0; j < S; j++){
    
    arma::mat Lambda_beta = (tXX / (tau[j] * tau[j])) + B_z / (sigma_beta * sigma_beta);
    arma::vec mu_beta = arma::zeros(ncov_z);
    
    for(int i = 0; i < X_beta_noempty.n_rows;i++){
      
      arma::rowvec Xbetal = X_beta_noempty.row(i) * l_noempty(i,j) / (tau[j] * tau[j]);
      arma::vec Xbetal_vec = arma::conv_to<arma::vec>::from(Xbetal);
      mu_beta += Xbetal_vec;
      
    }
    
    // Rcout << mu_beta << std::endl;
    
    // mu_beta <- apply(t(sapply(1:nrow(X_beta_noempty), function(i){
      //   X_beta_noempty[i,] * l_noempty[i,j] / tau[j]^2
      // })), 2, sum)
    
    arma::vec beta_zj = mvrnormArma(arma::inv(Lambda_beta) * mu_beta, arma::inv(Lambda_beta));
    beta_z.col(j) = arma::conv_to<arma::colvec>::from(beta_zj);
    // beta_z[,j] <- mvrnorm(1, solve(Lambda_beta) %*% mu_beta, solve(Lambda_beta))
    
  }
  
  return beta_z;
}

// [[Rcpp::export]]
arma::mat update_betaw_cpp(arma::mat beta_w, arma::mat v, arma::mat delta, 
                           arma::mat logz, arma::mat X_w, arma::vec sigma, 
                           double sigma_beta, arma::vec M_site){
  
  int S = beta_w.n_cols;
  int n = logz.n_rows;
  int ncov_w = X_w.n_cols;
  
  arma::mat prior_betaw = (1 / (sigma_beta * sigma_beta)) * arma::eye(ncov_w, ncov_w);
  
  for(int j = 0; j < S; j++){
    
    arma::mat Xw_long = arma::zeros(sum(M_site), ncov_w);
    arma::vec y_long = arma::zeros(sum(M_site));
    
    int l = 0;
    int r_idx = 0;
    for(int i = 0; i < n; i++){
      
      for (int m = 0; m < M_site[i]; m++) {
        
        if(delta(l, j) == 1){
          y_long[r_idx] = v(l, j) - logz(i,j);
          Xw_long.row(r_idx) = X_w.row(l);
          
          r_idx++;  
        }
        
        l++;
      }
      
    }
    
    arma::mat Xw_long2 = arma::zeros(r_idx, ncov_w);
    arma::vec y_long2 = arma::zeros(r_idx);
    arma::vec Xy_long2_colsums = arma::zeros(ncov_w);
    for(int l = 0; l < r_idx; l++){
      Xw_long2.row(l) = Xw_long.row(l);
      y_long2[l] = y_long[l];
      Xy_long2_colsums += arma::trans(Xw_long2.row(l)) * y_long[l];
    }
    
    arma::mat r_longtrlong = arma::trans(Xw_long2) * Xw_long2;
    arma::mat Lambda_beta = (r_longtrlong / pow(sigma[j],2)) + prior_betaw;
    // double mu_beta = sum(y_long % r_long) / pow(sigma[j],2);
    arma::vec mu_beta = Xy_long2_colsums / pow(sigma[j],2);
    
    arma::mat Lambdabeta = arma::inv(Lambda_beta);
    arma::vec mean_beta = Lambdabeta * mu_beta;
    
    arma::vec betaw_j = mvrnormArma(mean_beta, Lambdabeta);
    beta_w.col(j) = arma::conv_to<arma::colvec>::from(betaw_j);
    // alpha[j] = R::rnorm(mu_beta / Lambda_beta, 1 / Lambda_beta);
    
  }
  
  return(beta_w);
}

// DELTA GAMMA

// [[Rcpp::export]]
List updateDeltaGammaCpp(arma::vec Ct, 
                      arma::vec C_star, 
                      arma::vec w_pcr, 
                      arma::vec w_star,
                      double phi_0, 
                      double phi_1, 
                      double p0,
                      arma::vec alpha1,
                      arma::vec alpha2,
                      arma::vec sigma_y,
                      arma::vec P,
                      arma::vec P_star,
                      double lambda,
                      double sigma_lambda){
  
  int n = Ct.size();
  
  arma::vec delta = arma::zeros(n);
  arma::vec gamma = arma::zeros(n);
  
  for (int i = 0; i < n; i++) {
    
    if(Ct[i] > 0){
      
      double w_current = w_pcr[i];
      
      if(w_current > 0){
        
        double Xb_current = phi_0 + phi_1 * log(w_current);
        
        double logprior_delta1 = - log(1 + 1 / exp(Xb_current));
        
        double logprior_delta0 = - Xb_current - log(1 + 1 / exp(Xb_current));
        
        double logprior_gamma_1 = log(p0);
        
        double loglik_delta1 = R::dnorm(Ct[i], alpha1[P[i] - 1] + alpha2[P[i] - 1] * log(w_current), 
                               sigma_y[P[i] - 1], 1);
        
        double loglik_gamma1 = R::dnorm(Ct[i], lambda, sigma_lambda, 1);
        
        double log_delta1 = logprior_delta1 + loglik_delta1;
        
        double log_gamma1 = logprior_delta0 +  logprior_gamma_1 + loglik_gamma1;
        
        double p_delta1 = exp(log_delta1 - log_gamma1) / (1 + exp(log_delta1 - log_gamma1));
        
        delta[i] = R::rbinom(1, p_delta1);
        gamma[i] = 1 - delta[i];
        
      } else {
        
        delta[i] = 0;
        gamma[i] = 1;
        
      }
      
    } else {
      
      delta[i] = 0;
      gamma[i] = 0;
      
    }
    
  }
  
  int n_star = C_star.size();
  
  arma::vec delta_star = arma::zeros(n_star);
  arma::vec gamma_star = arma::zeros(n_star);
  
  for (int i = 0; i < n_star; i++) {
    
    if(C_star[i] > 0){
        
      double w_current = w_star[i];
        
        double Xb_current = phi_0 + phi_1 * log(w_current);
        
        double logprior_delta1 = - log(1 + 1 / exp(Xb_current));
        
        double logprior_delta0 = - Xb_current - log(1 + 1 / exp(Xb_current));
        
        double logprior_gamma_1 = log(p0);
        
        double loglik_delta1 = R::dnorm(C_star[i], alpha1[P_star[i] - 1] + alpha2[P_star[i] - 1] * log(w_current), 
                                     sigma_y[P_star[i] - 1], 1);
        
        double loglik_gamma1 = R::dnorm(C_star[i], lambda, sigma_lambda, 1);
        
        double log_delta1 = logprior_delta1 + loglik_delta1;
        
        double log_gamma1 = logprior_delta0 +  logprior_gamma_1 + loglik_gamma1;
        
        double p_delta1 = exp(log_delta1 - log_gamma1) / (1 + exp(log_delta1 - log_gamma1));
        
        delta_star[i] = R::rbinom(1, p_delta1);
        gamma_star[i] = 1 - delta_star[i];
        
      } else {
        
        delta_star[i] = 0;
        gamma_star[i] = 0;
        
      }
      
    }
  
  return List::create(_["delta"] = delta,
                      _["gamma"] = gamma,
                      _["delta_star"] = delta_star,
                      _["gamma_star"] = gamma_star);
  
}

////////////
// V ///
///////

// [[Rcpp::export]]
double logf_v_cpp(double v_i, 
                  arma::vec Ct_i, 
                  double v_mean, 
                  double sigma, 
                  arma::vec alpha1_i, 
                  arma::vec alpha2_i, 
                  double w_tilde_i, 
                  arma::vec sigma_y, 
                  double p0, 
                  double phi_0, 
                  double phi_1,
                  double lambda, 
                  double sigma_lambda){
  
  double logprior = R::dnorm(v_i, v_mean, sigma, 1);
  
  double loglik = 0;
  if(v_i > 0){
    
    double w_i = exp(v_i) - 1;
    
    for(int i = 0; i < Ct_i.size(); i++){
      
      if(Ct_i[i] > 0){
        
        double prob = logistic(phi_0 + phi_1 * log(w_i));
        
        loglik += log(
          R::dbinom(1, 1, prob, 0) *
                     R::dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_i), sigma_y[i], 0) + 
                             R::dbinom(0, 1, prob, 0) *
                                      p0 * R::dnorm(Ct_i[i], lambda, sigma_lambda, 0));
        
      } else {
        
        double prob = logistic(phi_0 + phi_1 * log(w_i));
        
        loglik += R::dbinom(0, 1, prob, 1) + log(1 - p0);
        
      }
      
    }
    
  } else {
    
    for(int i = 0; i < Ct_i.size(); i++){
      
      if(Ct_i[i] > 0){
        
        double prob = logistic(phi_0 + phi_1 * log(w_tilde_i));
        
        loglik += log(
          R::dbinom(1, 1, prob, 0) * R::dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_tilde_i), sigma_y[i], 0) + 
            R::dbinom(0, 1, prob, 0) * p0 * R::dnorm(Ct_i[i], lambda, sigma_lambda, 0));
        
      } else {
        
        double prob = logistic(phi_0 + phi_1 * log(w_tilde_i));
        
        loglik += R::dbinom(0, 1, prob, 1) + log(1 - p0);
        
      }
      
    }
    
  }
  
  return logprior + loglik;
    
}

// [[Rcpp::export]]
arma::vec updateV_cpp(arma::vec v, 
                      arma::vec Ct, 
                      arma::vec l, 
                      arma::vec M,
                      double beta_w0, 
                      arma::mat X_w, 
                      arma::vec beta_w,
                      arma::vec alpha1, 
                      arma::vec alpha2, 
                      arma::vec P, 
                      arma::vec w_tilde,
                      arma::vec numSampleV,
                      double sigma,
                      arma::vec sigma_y,
                      double p0,
                      double phi_0, 
                      double phi_1, 
                      arma::vec sigmas_v,
                      double lambda, 
                      double sigma_lambda){
  
  arma::vec Xwbetaw = beta_w0 + X_w * beta_w;
  int n = l.size();
    
  int index_m = 0;
  int index_k = 0;
  for(int i = 0; i < n; i++){
    
    double l_i = l[i];
    
    for (int m = 0; m < M[i]; m++) {
      
      int idx_current = m + index_m;
      int numSampleCurrent = numSampleV[idx_current];
      double v_mean = l_i + Xwbetaw[idx_current];
      
      double v_current = v[idx_current];
      
      arma::vec Ct_i = arma::zeros(numSampleCurrent);
      arma::vec alpha1_i = arma::zeros(numSampleCurrent);
      arma::vec alpha2_i = arma::zeros(numSampleCurrent);
      for(int l3 = 0; l3 < numSampleCurrent; l3++){
        Ct_i[l3] = Ct[index_k + l3];
        alpha1_i[l3] = alpha1[P[index_k + l3] - 1];
        alpha2_i[l3] = alpha2[P[index_k + l3] - 1];
      }
      
      // Ct_i <- Ct[idx_sample_K == idx_current]
      // alpha1_i <- alpha1[P[idx_sample_K == idx_current]]
      // alpha2_i <- alpha2[P[idx_sample_K == idx_current]]
      
      double w_tilde_i = w_tilde[idx_current];
      
      double v_star = R::rnorm(v_current, sigmas_v[idx_current]);
        
      double logposterior_current = logf_v_cpp(v_current, Ct_i, v_mean, sigma, 
                                       alpha1_i, alpha2_i, w_tilde_i, 
                                       sigma_y, p0, phi_0, phi_1,
                                       lambda, sigma_lambda);
        
      double logposterior_star = logf_v_cpp(v_star, Ct_i, v_mean, sigma, 
                                    alpha1_i, alpha2_i, w_tilde_i, 
                                    sigma_y, p0, phi_0, phi_1,
                                    lambda, sigma_lambda);

      if(R::runif(0, 1) < exp(logposterior_star - logposterior_current)){
        
        v[idx_current] = v_star;
        
      }
      
      index_k += numSampleCurrent;
      
    }
    
    index_m += M[i];
    
  }
  
  return v;
    
}

////////////
// V TILDE
///////

// [[Rcpp::export]]
double logf_v_tilde_cpp(double v_tilde_i, 
                  arma::vec Ct_i, 
                  double nu, 
                  double sigma_nu, 
                  arma::vec alpha1_i, 
                  arma::vec alpha2_i, 
                  arma::vec sigma_y, 
                  double p0, 
                  double phi_0, 
                  double phi_1,
                  double lambda, 
                  double sigma_lambda){
  
  double logprior = R::dnorm(v_tilde_i, nu, sigma_nu, 1);
  
  double loglik = 0;
  if(v_tilde_i > 0){
    
    double w_tilde_i = exp(v_tilde_i) - 1;
    
    for(int i = 0; i < Ct_i.size(); i++){
      
      if(Ct_i[i] > 0){
        
        double prob = logistic(phi_0 + phi_1 * log(w_tilde_i));
        
        loglik += log(
          R::dbinom(1, 1, prob, 0) *
                     R::dnorm(Ct_i[i], alpha1_i[i] + alpha2_i[i] * log(w_tilde_i), sigma_y[i], 0) + 
                             R::dbinom(0, 1, prob, 0) *
                                      p0 * R::dnorm(Ct_i[i], lambda, sigma_lambda, 0));
        
      } else {
        
        double prob = logistic(phi_0 + phi_1 * log(w_tilde_i));
        
        loglik += R::dbinom(0, 1, prob, 1) + log(1 - p0);
        
      }
      
    }
    
  } else {
    
    for(int i = 0; i < Ct_i.size(); i++){
      
      if(Ct_i[i] > 0){
        
        loglik += log(p0) + R::dnorm(Ct_i[i], lambda, sigma_lambda, 1);
        
      } else {
        
        loglik += log(1 - p0);
        
      }
      
    }
    
  }
  
  return logprior + loglik;
    
}

// [[Rcpp::export]]
arma::vec updateVtilde_cpp(arma::vec v_tilde,
                           arma::vec v,
                           arma::vec Ct,
                           arma::vec l,
                           arma::vec M,
                           double beta_w0,
                           arma::mat X_w,
                           arma::vec beta_w,
                           arma::vec alpha1,
                           arma::vec alpha2,
                           arma::vec P,
                           double nu,
                           double sigma_nu,
                           arma::vec numSampleV,
                           double sigma,
                           arma::vec sigma_y,
                           double p0,
                           double phi_0,
                           double phi_1,
                           arma::vec sigmas_v,
                           double lambda,
                           double sigma_lambda){

  int n = l.size();

  int index_m = 0;
  int index_k = 0;
  for(int i = 0; i < n; i++){

    double l_i = l[i];

    for (int m = 0; m < M[i]; m++) {

      int idx_current = m + index_m;
      int numSampleCurrent = numSampleV[idx_current];
      
      if(v[idx_current] > 0){
        
        double v_tilde_current = v_tilde[idx_current];
        
        arma::vec Ct_i = arma::zeros(numSampleCurrent);
        arma::vec alpha1_i = arma::zeros(numSampleCurrent);
        arma::vec alpha2_i = arma::zeros(numSampleCurrent);
        for(int l3 = 0; l3 < numSampleCurrent; l3++){
          Ct_i[l3] = Ct[index_k + l3];
          alpha1_i[l3] = alpha1[P[index_k + l3] - 1];
          alpha2_i[l3] = alpha2[P[index_k + l3] - 1];
        }
        
        // Ct_i <- Ct[idx_sample_K == idx_current]
        // alpha1_i <- alpha1[P[idx_sample_K == idx_current]]
        // alpha2_i <- alpha2[P[idx_sample_K == idx_current]]
        
        double v_tilde_star = R::rnorm(v_tilde_current, sigmas_v[idx_current]);
      
        double logposterior_current = logf_v_tilde_cpp(v_tilde_current, Ct_i, 
                                                       nu, sigma_nu,
                                                       alpha1_i, alpha2_i, 
                                                       sigma_y, p0, phi_0, 
                                                       phi_1,
                                                       lambda, sigma_lambda);
        
        double logposterior_star = logf_v_tilde_cpp(v_tilde_star, Ct_i, 
                                                    nu, sigma_nu,
                                                    alpha1_i, alpha2_i, 
                                                    sigma_y, p0, phi_0, 
                                                    phi_1,
                                                    lambda, sigma_lambda);
        
        if(R::runif(0, 1) < exp(logposterior_star - logposterior_current)){
          
          v_tilde[idx_current] = v_tilde_star;
          
        }  
        
      }
      

      index_k += numSampleCurrent;

    }

    index_m += M[i];

  }

  return v_tilde;

}

///////////// 
//// SIGMA Y
//////////////

// [[Rcpp::export]]
arma::vec updateSigmaY_cpp(arma::vec Ct, 
                    arma::vec w_pcr, 
                    arma::vec delta, 
                    arma::vec C_star, 
                    arma::vec w_star, 
                    arma::vec delta_star,
                    arma::vec alpha1, 
                    arma::vec alpha2,
                    double sigma_alpha, 
                    arma::vec P,
                    arma::vec P_star,
                    int n_P,
                    double a_sigma_y,
                    double b_sigma_y){
  
  arma::vec w_all = join_cols(w_pcr, w_star);  
  arma::vec delta_all = join_cols(delta, delta_star);  
  arma::vec Ct_all = join_cols(Ct, C_star);  
  arma::vec P_all = join_cols(P, P_star);  
  
  // int n_samples = 0;
  // double sumsq = 0;
  // for(int i = 0; i < w_all.size(); i++){
  //   
  //   if(delta_all[i] == 1){
  //     
  //     double w_i = w_all[i];
  //     double Xtalpha = alpha1[P_all[i] - 1] + alpha2[P_all[i] - 1] * log(w_i);
  //     double Ct_res = pow(Ct_all[i] - Xtalpha, 2);
  //     n_samples += 1;
  //     
  //     sumsq += Ct_res;
  //     
  //   }
  //   
  // }
  
  arma::vec n_samples_p = arma::zeros(n_P);
  arma::vec sumsq_p = arma::zeros(n_P);
  
  for(int i = 0; i < w_all.size(); i++){
    
    if(delta_all[i] == 1){
      
      int p_i = P_all[i];
      
      double w_i = w_all[i];
      double Xtalpha = alpha1[p_i - 1] + alpha2[p_i - 1] * log(w_i);
      double Ct_res = pow(Ct_all[i] - Xtalpha, 2);
      
      n_samples_p[p_i - 1] += 1;
      sumsq_p[p_i - 1] += Ct_res;
      
      // Xt_rowi[1] = log(w_all[i]);
      // 
      // arma::mat Xt_i = Xt_rowi * arma::trans(Xt_rowi);
      // 
      // Lambda_alpha.subcube(arma::span(p_i - 1), arma::span(), arma::span()) += Xt_i;
      // 
      // mu_alpha.row(p_i - 1) += arma::conv_to<arma::rowvec>::from(Xt_rowi * Ct_all[i]);
      
    }
    
  }
  
  arma::vec sigma_y = arma::zeros(n_P);
  
  for(int p = 0; p < n_P; p++){
    
    sigma_y[p] = sqrt(rinvgamma_cpp(a_sigma_y + n_samples_p[p] / 2.0, b_sigma_y + sumsq_p[p] / 2.0));
    
  }
  // double sigma_y = sqrt(rinvgamma_cpp(a_sigma_y + n_samples / 2.0, b_sigma_y + sumsq / 2.0));
  
  return sigma_y;
  
}


///////////// 
//// PHI
//////////////

// [[Rcpp::export]]
double loglik_phi_cpp(arma::vec phi_01,
                      arma::vec X_phi, 
                      arma::vec delta_w){
  
  // Xphi <- X_phi %*% phi_01
  
  double loglik = 0;
  
  for(int i = 0; i < delta_w.size(); i++){
    
    double Xphi_i = phi_01[0] + phi_01[1] * X_phi[i];
    
    if(delta_w[i] == 1){
      
      loglik += (- log(1 + exp(- Xphi_i)));
      
    } else {
      
      loglik += (- Xphi_i - log(1 + exp(- Xphi_i)));
      
    }
    
  }
  
  // sum(sapply(1:nrow(X_phi), function(i){
  //   
  //   if(delta_w[i] == 1){
  //     
  //     return(- log(1 + exp(- Xphi[i])))
  //     
  //   } else {
  //     
  //     return(- Xphi[i] - log(1 + exp(- Xphi[i])))
  //     
  //   }
  //   
  // }))
  
  return loglik;
  
}


// // [[Rcpp::export]]
// double logf_cpp(double l,
//                 double x,
//                 double sigmaj,
//                 arma::vec v_samples,
//                 arma::vec v,
//                 arma::vec y,
//                 double tauj,
//                 double prior_mean){
//   
//   double loglikelihood = - 1 / ( pow(sigmaj,2)) * 
//     sum(- (v_samples - l));
//   
//   double x_all_l = l * x;
//   // double x_all_l = exp(l) * x;
//   
//   // double loglikelihood_v = 0;
//   // for(int i = 0; i < y.size(); i++){
//   //   loglikelihood_v += y[i] * x_all_l -
//   //     x_all_l * exp(v[i] + x_all_l) / (1 + exp(v[i] + x_all_l));
//   // }
//   
//   double loglikelihood_v = 0;
//   for(int i = 0; i < y.size(); i++){
//     if(y[i] == 1){
//       loglikelihood_v += ( x * exp(-v[i] - x_all_l)) / (1 + exp(-v[i] - x_all_l));
//     } else {
//       loglikelihood_v += (- x) / (1 + exp(-v[i] - x_all_l));
//     }
//     // loglikelihood_v += x * (y[i] * exp(l - v[i] - x * exp(l)) +
//     //   (y[i] - 1) * exp(l)) / (1 + exp(-v[i] - x * exp(l)));
//   }
//   
//   double logprior = - (l - prior_mean) / (pow(tauj,2)); 
//   
//   return( loglikelihood + loglikelihood_v + logprior);
// }
// 
// // [[Rcpp::export]]
// double findzero_cpp(double a,
//                     double b,
//                     double tol,
//                     double x_all,
//                     arma::vec y_all,
//                     arma::vec v_all,
//                     arma::vec v_samples,
//                     double tauj,
//                     double sigma_j,
//                     double prior_mean){
//   
//   double c = (a + b) / 2;
//   
//   double fc = logf_cpp(c,
//                        x_all,
//                        sigma_j,
//                        v_samples,
//                        v_all,
//                        y_all,
//                        tauj,
//                        prior_mean);
//   
//   int nsteps = 0;
//   
//   while( (b - a) / 2 > tol  & nsteps < 50){
//     
//     double fa = logf_cpp(a,
//                          x_all,
//                          sigma_j,
//                          v_samples,
//                          v_all,
//                          y_all,
//                          tauj,
//                          prior_mean);
//     
//     if((fc < 0 & fa < 0) | (fc > 0 & fa > 0)){
//       a = c;
//     } else {
//       b = c;
//     }  
//     
//     c = (a + b) / 2;
//     
//     fc = logf_cpp(c,
//                   x_all,
//                   sigma_j,
//                   v_samples,
//                   v_all,
//                   y_all,
//                   tauj,
//                   prior_mean);
//     
//     nsteps++;
//   }
//   
//   return (a + b) / 2;
// }

// [[Rcpp::export]]
arma::mat hes_loglik_phi_cpp(arma::vec phi_01,
                  arma::vec X_phi, 
                  arma::vec delta_w){

  arma::mat hess = arma::zeros(2, 2);
  
  for(int i = 0; i < delta_w.size(); i++){
    
    double expXphi_i = exp(-phi_01[0] - phi_01[1] * X_phi[i]);
    
    hess(0, 0) += (- expXphi_i / pow(1 + expXphi_i, 2));
    hess(1, 0) += X_phi[i] * (- expXphi_i / pow(1 + expXphi_i, 2));
    hess(1, 1) += X_phi[i] * X_phi[i] * (- expXphi_i / pow(1 + expXphi_i, 2));
    
  }
  
  hess(0, 1) = hess(1, 0);

  return(hess);

}
