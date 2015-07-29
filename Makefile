default:
	g++ -O2 -I/usr/lib/R/site-library/RcppEigen/include -I/usr/lib/R/site-library/Rcpp/include -I/usr/share/R/include Likelihood.cpp Data.cpp MyModel.cpp
	g++ -o -L/usr/lib/R/site-library/RcppEigen/libs main main.o -ldnest3 -lrjobject -lboost_thread -lboost_system -lgsl -lgslcblas -lR -lRcppEigen
	rm -f *.o

