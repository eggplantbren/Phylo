#ifndef _MyModel_
#define _MyModel_

#include "Model.h"
#include <vector>
#include <RJObject.h>
#include "MyDistribution.h"

class MyModel:public DNest3::Model
{
	private:
		double mu;			// Branch length hyperparameter
		RJObject<MyDistribution> branch_lengths;

		double phi;			// Transition rate hyperparameter
		RJObject<MyDistribution> transition_rates;

		std::vector<double> frequencies;
		double lambda;


	public:
		MyModel();

		// Generate the point from the prior
		void fromPrior();

		// Metropolis-Hastings proposals
		double perturb();

		// Likelihood function
		double logLikelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;
};

#endif

