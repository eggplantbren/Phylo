#include <RcppEigen.h>
#include <Rcpp.h>
//#include <Rmath.h>
using namespace Rcpp;

#define set0M(x) (x < 0.0 ? 0.0 : x);

#define priorMacro(a, b, mu, al, bl, lambda, t, sumRates, phi) ( -(a + 18.0) * log( mu ) + ( al - 1.0 ) * log( lambda ) - (t + b) / mu - lambda / bl + 5.0 * log(phi) - phi * (sumRates + 1.0) );

#define logplusC(x,y) (x>y ? x+log(1.0+exp(y-x)) : y+log(1.0+exp(x-y)) );

// [[Rcpp::export]]  
double logplusvecC( const NumericVector & x){
  
  int    n = x.size();
  
  double r = -DBL_MAX;
  
  for(int i = 0; i < n; i++){
    
      r = logplusC(r, x[i]);
    
  }
  
  return r;
  
}

// Summing 2 columns of 2 matrices 4x4
// [[Rcpp::export]]  
NumericVector sumcolC(const NumericMatrix & P1, const NumericMatrix & P2, 
                      const int & x, const int & y){
  
  NumericVector out(4);
  
  for(int i = 0; i < 4; i++){
    
    out[i] = P1(i,x) + P2(i,y);
    
  }
  
  //return Rcpp::wrap(out);            
  return out;
}

// Summing 2 matrices
// [[Rcpp::export]]
NumericMatrix sumMatC( const NumericMatrix & x, const NumericMatrix & y){
  
  int n = x.nrow(), m = x.ncol();
  
  NumericMatrix out(n,m);
  
  for(int j = 0; j < m; j++){
    
    for(int i = 0; i < n; i++){
      
      out(i, j) = x(i, j) + y(i, j);
      
    }
  }
  
  return out;
  
}

// Summing matrix and vector
// [[Rcpp::export]]
NumericMatrix sumMatVecC(NumericMatrix x, NumericVector y ){
  
  int n = x.nrow(), m = x.ncol();
  
  NumericMatrix out(n,m);
  
  for(int j = 0; j < m; j++){
    for(int i = 0; i < n; i++){
      
      out(i, j) = x(i, j) + y[j];
      
    }
  }
  
  return out;
  
}

// Extracting a k-th row of matrix "x"
// [[Rcpp::export]]
NumericVector ExtRowMatC(const NumericMatrix & x, const int & k){
  
  int m = x.ncol();
  
  if( (k < 0) || (k > x.nrow() - 1) ){
    
    Rcpp::Rcout << "ExtRowMatC: That row does not exist" << std::endl;
    
    return NumericVector::create(NA_REAL);
    
  }
  
  NumericVector out(m);  
  
  for(int j = 0; j < m; j++){
    
    out(j) = x(k, j);
    
  } 
  
  return out;  
  
} 

// Extracting a k-th column of matrix "x"
// [[Rcpp::export]]
NumericVector ExtColMatC(const NumericMatrix & x, const int & k){
  
  int n = x.nrow();  
  
  if( (k < 0) || (k > x.ncol() - 1) ){
    
    Rcpp::Rcout << "ExtColMatC: That column does not exist" << std::endl;
    
    return NumericVector::create(NA_REAL);
    
  }
  
  NumericVector out(n);
  
  for(int j = 0; j < n; j++){
    
    out(j) = x(j, k);
    
  } 
  
  return out;  
  
}

// Summing log row
// [[Rcpp::export]]
NumericVector logplusRowMatC(const NumericMatrix & x){
  
  int n = x.nrow();
  
  NumericVector out(n);
  
  for(int i = 0; i < n; i++){
    
    out(i) = logplusvecC( ExtRowMatC(x, i) );
    
  }
  
  return out;  
  
}

// Summing 2 vectors
// [[Rcpp::export]]
NumericVector sumVecC( NumericVector x, NumericVector y){
  
  int n = x.size(), m = y.size();
  
  if( n != m){
    
    Rcpp::Rcout << "sumVecC: the vectors must have the same dimension" << std::endl;
    
    return NumericVector::create(NA_REAL);
    
  }
  
  NumericVector out(n);
  
  for(int i = 0; i < n; i++){
    
    out(i) = x(i) + y(i);
    
  }
  
  return out;
  
}

// log-likelihood of a site
// [[Rcpp::export]]          
double loglikesiteC(const NumericVector & data, const NumericMatrix & P1, const NumericMatrix & P2,
                    const NumericMatrix & P3, const NumericMatrix & P4, const NumericMatrix & P5,
                    const NumericMatrix & P6, const NumericMatrix & P7, const NumericMatrix & P8,
                    const NumericMatrix & P9, const NumericMatrix & P10, const NumericMatrix & P11,
                    const NumericMatrix & P12, const NumericMatrix & P13, const NumericMatrix & P14,
                    const NumericMatrix & P15, const NumericMatrix & P16, const NumericMatrix & P17,
                    const NumericVector & lbf){ 
  
  NumericVector n11(4), n12(4), n13(4), n14(4), n15(4), n16(4), n17(4), n18(4);
  
  n11    = sumcolC(P1, P2, data(0), data(1));    
  
  n13    = sumVecC( logplusRowMatC( sumMatVecC(P3, n11) ), ExtColMatC(P4, data(2)) ); 
  
  n14    = sumVecC( logplusRowMatC( sumMatVecC(P5, n13) ), ExtColMatC(P6, data(3)) );
  
  n12    = sumcolC(P9, P10, data(5), data(4));
  
  n15    = sumVecC( logplusRowMatC( sumMatVecC(P7, n14) ), 
                    logplusRowMatC( sumMatVecC(P8, n12) ) );   
  
  n16    = sumVecC( logplusRowMatC( sumMatVecC(P11, n15)), ExtColMatC(P12, data(6)) );  
  
  n17    = sumVecC( logplusRowMatC( sumMatVecC(P13, n16)), ExtColMatC(P14, data(7)) );
  
  n18    = sumVecC( logplusRowMatC( sumMatVecC(P15, n17)), ExtColMatC(P16, data(8)) );
  
  n18    = sumVecC(n18, sumVecC(ExtColMatC(P17, data(9)), lbf) );
  
  return logplusvecC(n18);
  
}

// log-likelihood of sites
// [[Rcpp::export]] 
NumericVector loglikesitesC(NumericMatrix data, const NumericMatrix & P1, const NumericMatrix & P2,
                            const NumericMatrix & P3, const NumericMatrix & P4, const NumericMatrix & P5,
                            const NumericMatrix & P6, const NumericMatrix & P7, const NumericMatrix & P8,
                            const NumericMatrix & P9, const NumericMatrix & P10, const NumericMatrix & P11,
                            const NumericMatrix & P12, const NumericMatrix & P13, const NumericMatrix & P14,
                            const NumericMatrix & P15, const NumericMatrix & P16, const NumericMatrix & P17,
                            NumericVector lbf){
  int nsite = data.ncol();                         
  
  NumericVector out(nsite);
  
  for(int j = 0; j < nsite; j++){
    
    out(j) = loglikesiteC(ExtColMatC(data, j), P1, P2, P3, P4, P5, P6, P7, P8,
                             P9, P10, P11, P12, P13, P14, P15, P16, P17, lbf);
  }
  
  return out;
  
}

// sequence from 0 to 1, equally spaced (k+1 elements)
// [[Rcpp::export]]
NumericVector seqC(int k){
  
  NumericVector out(k+1);
  
  out(0)   = 0.0;
  out(k)   = 1.0;
  
  for(int i = 1; i < k; i++ ){
    
    out(i) =  i / double(k); 
    
  }
  
  return out;
  
}

// [[Rcpp::export]]
NumericVector DiscreteGammaC(int k, double lambda){
  
  if (k == 1) {
    
    NumericVector out = NumericVector::create(1.0);
    
    return out;
    
  }
  
  NumericVector quants = qgamma(seqC(k), lambda, 1.0/lambda);
  
  quants = pgamma(quants * lambda, lambda + 1.0);            
  
  return diff(quants) * double(k);
  
}

// logarithm of each value of a matrix
// [[Rcpp::export]]
NumericMatrix LogMatC(const NumericMatrix & x){
  
  int nrow = x.nrow(), ncol = x.ncol();
  
  NumericMatrix out(nrow, ncol);
  
  if( is_true( any(x  < 0.0) ) ){
    
    Rcpp::Rcout << "LogMatC: the values must be positive" << std::endl;
    
    int xsize = out.nrow() * out.ncol();
    for (int i = 0; i < xsize; i++) {
      out[i] = NA_REAL;
    }
    
    return out;
  }
  
  for (int j = 0; j < ncol; j++){
    
    for (int i = 0; i < nrow; i++) {
      
      out(i, j) = log( x(i, j) );
      
    }}
  
  return out;  
  
}

// logarithm of each value of a vector
// [[Rcpp::export]]
NumericVector LogVecC(const NumericVector & x){
  
  int xsize = x.size();
  
  NumericVector out(xsize);
  
  if( is_true( any(x  < 0.0) ) ){
    
    Rcpp::Rcout << "LogVecC: the values must be positive" << std::endl;
    return out;
    
  }
  
  for (int i = 0; i < xsize; i++) {
    
    out[i] = log(x[i]);
    
  }
  
  return out;  
  
}

// [[Rcpp::export]]
NumericMatrix matQC(const NumericVector & x, const NumericVector & bf){
  
  if( x.size() != 6){
    Rcpp::Rcout << "matQC: include 6 rates" << std::endl;
  }
  
  NumericMatrix out(4, 4);
  out(1,0) = x[0] * bf[0]; out(0,1) = x[0] * bf[1]; 
  out(2,0) = x[1] * bf[0]; out(0,2) = x[1] * bf[2];
  out(3,0) = x[2] * bf[0]; out(0,3) = x[2] * bf[3];
  out(2,1) = x[3] * bf[1]; out(1,2) = x[3] * bf[2];
  out(3,1) = x[4] * bf[1]; out(1,3) = x[4] * bf[3];
  out(3,2) = x[5] * bf[2]; out(2,3) = x[5] * bf[3];
  
  for(int i = 0; i < 4; i++){
    out(i,i) = - sum( ExtRowMatC(out, i) );
  }
  
  double scale = 0;
  
  for(int i = 0; i < 4; i++){
    
   scale = out(i, i) * bf[i]+ scale;  // calculating transition rate mean
    
  }
  
  scale = - scale;                    // scaling matrix
  
  for(int j = 0; j < 4; j++){
    
    for(int  i = 0; i < 4; i++){
      
      out(i, j) = out(i, j) / scale;
      
    }
    
  }
  
  return out;  
  
}

// log prior density
// [[Rcpp::export]]
double logpriorC(const NumericVector & t, const double & mu, const double & lambda, 
                 const NumericVector & q, const double & phi, const double & a, 
                 const double & b, const double & al, const double & bl ){
  /*                  
  if( is_true( any(t < 0) ) ) return( R_NegInf );
  if( mu < 0 ) return( R_NegInf );
  if( lambda < 0 ) return( R_NegInf );
  */
  double sumt = sum(t);
  
  double sumRates = sum(q) - 1.0;   // minus the last value
  
  return priorMacro(a, b, mu, al, bl, lambda, sumt, sumRates, phi);  
  
}

////////////////////////
// Exponential Matrix //
////////////////////////

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;                       // 'maps' rather than copies 
using Eigen::MatrixXd;                  // variable size matrix, double precision

// Exponential Matrix in Rcpp format
// [[Rcpp::export]]
NumericMatrix expmC(const NumericMatrix & x){
  const Map<MatrixXd> xx(Rcpp::as<Map<MatrixXd> >(x));
  return Rcpp::wrap( xx.exp() );
}

////////////////////////
////////////////////////
////////////////////////

// Scalar times a matrix
// [[Rcpp::export]]
NumericMatrix scalarXmatrixC(const double & k, const NumericMatrix & x){
  
  int nrow = x.nrow(), ncol = x.ncol();
  
  NumericMatrix out(nrow, ncol);
  
  for(int j = 0; j< ncol; j++){
    
    for(int i = 0; i< nrow; i++){
      
      out(i, j) = x(i, j) * k;
      
    }
  }
  
  return out;
  
}

// [[Rcpp::export]]
NumericMatrix matPC( const double & t, const NumericMatrix & x){
  
  if( (t  < 0.0) ){
    Rcpp::Rcout << "Error: t must be positive " <<std::endl;
  }
  
  int nrow = x.nrow(), ncol = x.ncol();
  
  NumericMatrix out(nrow, ncol);
  
  for(int j = 0; j< ncol; j++){
    
    for(int i = 0; i< nrow; i++){
      
      out(i, j) = x(i, j) * t;
      
    }
  }
  
  out  = expmC(out);
  
  for (int j = 0; j < ncol; j++){
    
    for (int i = 0; i < nrow; i++){
      
      out(i, j) = log( out(i, j) );
      
    }
  }
  
  return out;
  
}

// [[Rcpp::export]]
double loglikeC( const NumericMatrix & data, const NumericVector & t, const NumericVector & q, 
                 const NumericVector & bf, const int & k, const double & lambda){
  
  if( is_true( any(t  < 0) ) ||     // checking inputs
      is_true( any(q  < 0) ) ||
      is_true( any(bf < 0) ) ||
      ( lambda < 0 ) )  return( R_NegInf );
  
  int nsite         = data.ncol();
  
  NumericVector lbf = LogVecC(bf), 
                r   = DiscreteGammaC(k, lambda); 
  
  NumericMatrix Q   = matQC(q, bf), out(nsite, k), 
                      rQ, P1, P2, P3, P4, P5, P6, P7, P8, P9, P10,
                      P11, P12, P13, P14, P15, P16, P17;
  
  for(int j = 0; j < k ; j++){
    
    rQ  = scalarXmatrixC(r[j], Q);  
    P1  = matPC( t[0], rQ);
    P2  = matPC( t[1], rQ);
    P3  = matPC( t[2], rQ);
    P4  = matPC( t[3], rQ);
    P5  = matPC( t[4], rQ);
    P6  = matPC( t[5], rQ);
    P7  = matPC( t[6], rQ);
    P8  = matPC( t[7], rQ);
    P9  = matPC( t[8], rQ);
    P10 = matPC( t[9], rQ);
    P11 = matPC( t[10], rQ);
    P12 = matPC( t[11], rQ);
    P13 = matPC( t[12], rQ);
    P14 = matPC( t[13], rQ);
    P15 = matPC( t[14], rQ);
    P16 = matPC( t[15], rQ);
    P17 = matPC( t[16], rQ);  
    
    for(int i = 0; i < nsite; i++){
      
      out(i, j) = loglikesiteC(ExtColMatC(data, i), P1, P2, P3, P4, P5, P6, 
          P7, P8, P9, P10, P11, P12, P13, P14, P15, 
          P16, P17, lbf); //- log( k )
    }
    
  }
  
  //NumericVector out1 = logplusRowMatC( out ) - log(k);
  
  //double out2 = sum( ExtRowMatC(data, 10) * out1 );
  
  //return out2;
  
  return sum( ExtRowMatC(data, 10) * (logplusRowMatC(out) - log(k)) );
  
}

// [[Rcpp::export]]
NumericVector getSetVecC( const NumericVector & x, const int & a, const int & b){
  NumericVector out(b - a + 1);
  for(int i = a; i < b+1; i++){
    out[i-a] = x[i];
  }
  return out;
}

// theta = ( t, mu, lambda, freq, betas ) 

// a = shape parameter for the mean of the branch lengths
// b = scale parameter for the mean of the branch lengths
// al = shape parameter in the Gamma parameter
// bl = scale parameter in the Gamma parameter
// theta = c(t, mu, lambda, bf, phi, q)
// t      = getSetVecC(theta, 0, 16)    // Branch length
// mu     = theta[17]                   // Mean of branch length 
// lambda = theta[18]                   // shape parameter  
// bf     = getSetVecC(theta, 19, 22)   // frequencies
// phi    = theta[23]
// q      = getSetVecC(theta, 23, m-1)  // values rate matrix 
// length of theta

// power posterior
// [[Rcpp::export]]
NumericVector lpwC(const NumericMatrix & data, NumericVector vec,
                   const double & a, const double & b, const double & al, 
                   const double & bl, const int & k, const double & d){
  
  NumericVector t  = getSetVecC(vec, 0, 16),
                q  = getSetVecC(vec, 24, 29),
                bf = getSetVecC(vec, 19, 22);
  double        mu = vec[17],
            lambda = vec[18], 
               phi = vec[23];
  
  if( ( is_true( any(t  <  0.0 ) ) ) ||
      ( is_true( any(q  <  0.0 ) ) ) ||      
      ( is_true( any(bf <= 0.0 ) ) ) ||
      ( is_true( any(bf >= 1.0 ) ) ) ||
      ( phi < 0.0) || ( mu < 0.0 ) || 
      ( lambda < 0.0 ) ){
    
    NumericVector out = NumericVector::create(R_NegInf, R_NegInf);
    
    return out;
    
  }else{    
    
    double like, prior, pw; 
    
    prior = logpriorC(t, mu, lambda, q, phi, a, b, al, bl);
    
    like  = loglikeC(data, t, q, bf, k, lambda ); 
    
    pw    = d * like + prior;
    
    NumericVector out = NumericVector::create(pw, like);
    
    return out;
  }
}