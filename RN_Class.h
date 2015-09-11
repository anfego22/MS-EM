#ifndef RNCLASS_H
#define RNCLASS_H

#include <boost/random/cauchy_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/chi_squared_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <Eigen/Dense>
#include <cmath>

using boost::random::cauchy_distribution;
using boost::random::chi_squared_distribution;
using boost::random::uniform_real_distribution;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class randomNb{
public:

	static boost::random::mt19937 rng;
	randomNb();
	static double sampleCauchy(double);
	static double sampleChiSqr(int);
	static double sampleUniform(double);
	static double varCovRg(double);
};

double multivariateNormal_density(const MatrixXd &,
								  const MatrixXd &,
								  const MatrixXd &);

void multivariateNormal_density(VectorXd &,
								const MatrixXd &,
								const MatrixXd &,
								const MatrixXd &);

#endif
