#include "MyModel.h"
#include "RandomNumberGenerator.h"
#include "Utils.h"
#include "Data.h"
#include <cmath>
#include <Rcpp.h>

using namespace Rcpp;
using namespace std;
using namespace DNest3;

MyModel::MyModel()
:branch_lengths(1, 17, true, MyDistribution(0., 50000.))
,transition_rates(1, 5, true, MyDistribution(0., 50000.))
,frequencies(4)
{

}

void MyModel::fromPrior()
{
	// "Good" prior: shape, scale
	lambda = rgamma(1, 10.0, 0.026)[0];
}

double MyModel::perturb()
{
	double logH = 0.;

	lambda = log(lambda);
	logH -= 10.0*lambda - exp(lambda)/0.026;
	lambda += randh();
	logH += 10.0*lambda - exp(lambda)/0.026;
	lambda = exp(lambda);

	return logH;
}

double MyModel::logLikelihood() const
{
	// Get the data
//	const vector<double>& x = Data::get_instance().get_x();

	double logL = 0;

	return logL;
}

void MyModel::print(std::ostream& out) const
{
	out<<lambda<<' ';
}

string MyModel::description() const
{
	return string("");
}

