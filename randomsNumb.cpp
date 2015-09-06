#include "RN_Class.h"

using boost::random::cauchy_distribution;
using boost::random::uniform_real_distribution;
using boost::random::chi_squared_distribution;
using boost::random::mt19937;
using boost::math::normal;
using namespace std;


// Static members should be initialize as:
// and static function?
mt19937 randomNb::rng;

randomNb::randomNb(){
	rng.seed(time(NULL));
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
