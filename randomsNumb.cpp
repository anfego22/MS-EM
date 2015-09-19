#include "RN_Class.h"

// Static members should be initialize as:
// and static function?

boost::random::mt19937 randomNb::rng;

randomNb::randomNb(){
	randomNb::rng.seed(time(NULL));
}

double randomNb::sampleCauchy(double c){
	cauchy_distribution <> nd(c, 1);
	return nd(randomNb::rng);
}

double randomNb::sampleChiSqr(int c){
	chi_squared_distribution <> nd(c);
	return nd(randomNb::rng);
}

double randomNb::sampleUniform(double c){
	uniform_real_distribution <> nd(0, 1);
  return nd(randomNb::rng);
}

double randomNb::varCovRg(double c){
	if (c == 0){
		return randomNb::sampleCauchy(c);
	}
	else{
		return randomNb::sampleChiSqr(c);
	}
}

double double_multivariateNormal_density(const MatrixXd &x,
	const MatrixXd &Mu, const MatrixXd &Sigma){
	double expTerm, o;
	using namespace std;
	MatrixXd err;
	err = (x - Mu);
	err*=Sigma.inverse();
	err *=(x - Mu).transpose();
	expTerm = exp(err(0));
	o = 1/sqrt(pow(2*M_PI, Mu.rows())*Sigma.determinant());
	return o*expTerm;
}

VectorXd multivariateNormal_density(const MatrixXd &x,
									const MatrixXd &Mu, const MatrixXd &Sigma){
	if ( x.rows() != Mu.rows())
		std::cerr << "ERROR: X_ROWS != MU_ROWS" << std::endl;
	VectorXd res;
	int T = x.rows();
	res.setZero(T);
	for(int i = 0; i<x.rows(); i++){
		res(i) = double_multivariateNormal_density(x.row(i),
												   Mu.row(i),
												   Sigma);
	}
	return res;
}
