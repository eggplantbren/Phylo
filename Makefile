default:
	g++ -c -O2 -I/usr/lib/R/site-library/RcppEigen/include -I/usr/lib/R/site-library/Rcpp/include -I/usr/share/R/include Likelihood.cpp Data.cpp
