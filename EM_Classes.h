#ifndef EMCLASSES_H
#define EMCLASSES_H

#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <stdexcept>
#include <iostream>
#include <thread>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

class Model {
public:
	int N, lagS, lagY;
	bool sigma, beta;
	Model(const int &, const int &, const int &,
		  const bool, const bool);
};

class Data {
public:
	MatrixXd y, x;
	double T, m, mX;
	Data(const MatrixXd & ,
		 const MatrixXd & );
	Data(const MatrixXd & );
};

class EM{
public:
	int N, M, K;
	double convLevel, delta;
	vector<double> likelihood_t;
	vector<MatrixXd> beta, sigma, beta0, sigma0,
		Xit_t, Xit1_t, XiS;
	MatrixXd P, Rho, P0, Rho0, Y, X;
	EM(const Model & ,const Data &,
	   const double &);
	void makeX();
	void setInitialParameters();
	void filterProb();
	void smoothProb();
	void newRho();
	void newP();
	void newBeta();
	void newSigma();
	void newBetaNS();
	void newSigmaNS();
	void iterEM();
};

#endif

