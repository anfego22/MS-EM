#ifndef EMCLASSES_H
#define EMCLASSES_H
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>

class randomNb{
	boost::random::mt19937 rng;
	randomNb();
	double sampleCauchy();
	double sampleChiSqr();
	double sampleUniform();
};



#endif
