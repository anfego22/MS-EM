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


double multivariateNormal_density(const MatrixXd &x,
	const MatrixXd &Mu, const MatrixXd &Sigma){
	double expTerm, o;
	using namespace std;
	MatrixXd err;
	if ( x.rows() != 1)
		cerr << "X must be a row vector" << endl;
	if ( Mu.rows() != 1)
		cerr << "Mu must be a row vector" << endl;
	err = x - Mu;
	expTerm = exp((err*Sigma.inverse()*err.transpose())(0,0));
	o = 1/sqrt(pow(2*M_PI, Mu.rows())*Sigma.determinant());
	return o*expTerm;
}
