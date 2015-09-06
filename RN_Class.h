#ifndef EMCLASSES_H
#define EMCLASSES_H
#include <boost/math/distributions/cauchy.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/uniform.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random.hpp>

using boost::random::mt19937;

class randomNb{
public:
	static mt19937 rng;
	randomNb();
	static double sampleCauchy(double);
	static double sampleChiSqr(int);
	static double sampleUniform(double);
	static double varCovRg(double);
};



#endif
