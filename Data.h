#ifndef _Data_
#define _Data_

#include <Rcpp.h>

class Data
{
	private:
		Rcpp::NumericMatrix values;

	public:
		Data();
		void load(const char* filename);

		// Getters
		const Rcpp::NumericMatrix& get_values() const { return values; }

	// Singleton
	private:
		static Data instance;
	public:
		static Data& get_instance() { return instance; }
};

#endif

