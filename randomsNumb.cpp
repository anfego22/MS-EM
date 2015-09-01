#include "RN_Class.h"

using boost::random::cauchy_distribution;
using boost::random::uniform_real_distribution;
using boost::random::chi_squared_distribution;
using boost::random::mt19937;
using boost::math::normal;
using namespace std;


randomNb::randomNb(){
	boost::random::mt19937 rng(time(0))
}

double randomNb::sampleCauchy(double c){
  boost::random::cauchy_distribution <> nd(c, 1);
  return nd(rng);
}

double randomNb::sampleChiSqr(int c){
  boost::random::chi_squared_distribution <> nd(c);
  return nd(rng);
}

double randomNb::sampleUniform(double c){
  boost::random::uniform_real_distribution <> nd(0, 1);
  return nd(rng);
}
