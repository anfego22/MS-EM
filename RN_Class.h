#ifndef RNCLASS_H
#define RNCLASS_H

#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>


using boost::random::cauchy_distribution;
using boost::random::chi_squared_distribution;
using boost::random::uniform_real_distribution;

class randomNb{
public:

	static boost::random::mt19937 rng;
	randomNb();
	static double sampleCauchy(double);
	static double sampleChiSqr(int);
	static double sampleUniform(double);
	static double varCovRg(double);
};



#endif
