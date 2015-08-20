default:
	g++ -c -O2 -I/usr/lib/R/site-library/RcppEigen/include -I/usr/lib/R/site-library/Rcpp/include -I/usr/share/R/include main.cpp Likelihood.cpp Data.cpp MyModel.cpp
	g++ -L/usr/lib/R/site-library/RcppEigen/libs -o main main.o -lrjobject -ldnest3 -lboost_thread -lboost_system -lgsl -lgslcblas -lR
	rm -f *.o

