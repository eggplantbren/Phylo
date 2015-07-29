#include "Data.h"
#include <fstream>
#include <iostream>

using namespace Rcpp;
using namespace std;

Data Data::instance;

Data::Data()
:values(11, 393)
{

}

void Data::load(const char* filename)
{
	fstream fin(filename, ios::in);
	if(!fin)
	{
		cerr<<"# Error. Couldn't open file "<<filename<<"."<<endl;
		return;
	}

	for(int i=0; i<values.nrow(); i++)
		for(int j=0; j<values.ncol(); j++)
			fin>>values(i, j);

	fin.close();
}

