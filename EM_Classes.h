#ifndef EMCLASSES_H
#define EMCLASSES_H
#include <sstream>
#include <string>
#include <vector>
#include "Eigen/Dense"
#include <iostream>


using Eigen::VectorXd;
using Eigen::MatrixXd;
using namespace std;

class Model {
public:
	int N, Nm, lagsY;
	bool sigma, beta, meanCorrected;
	Model(const int &, const int &, const bool &,
		  const bool &, const bool &);
	// N is the original number of regimes
	// Nm is N^{m+1} ficticial regimes
};

class Data {
public:
	MatrixXd Y, X;
	double T, mX, M;
	Data(const MatrixXd & ,
		 const MatrixXd & );
	Data(const MatrixXd & );
};

class Transitions{
	int N, Nm; 
	MatrixXd P, F;
	VectorXd rho;
	Transitions(const Model &);
	void createF();
	void createRho();
	void updateF(const MatrixXd & ,const MatrixXd &,
				 const MatrixXd &);
	void updateRho(const MatrixXd &);
};

class linearParams{
	const Data  &data;
	const Model &model;
	vector<MatrixXd> beta, sigma;
	linearParams(const Model &, const Data &);
	void createB();
	void createS();
	void update();
	void updateBetaS(); // Update if beta has switching
	void updateBetaNS(); // Update if beta has no switching
	void updateSigmaS(); // Update if sigma has switching
	void updateSigmaNS(); // Update if sigma has no switching
};

class Parameters{
	Transitions rhoF;
	linearParams betaSigma;
	// M -> Number of dependent variables (dimensions of Y)
	Parameters(const Model &model, const int & M);
};

class Errors {
	MatrixXd eta;
	Errors(const Data &d, const Model &model, const Parameters &param); 
};

class Xi {
	MatrixXd Xitt, Xit1t, XiS;
	Xi(const Errors &, const Parameters &);
	void filterProb(const Errors &, const Parameters &);
	void smoothProb(const Errors &, const Parameters &);
};

class EM{
public:
	double delta;
	Xi xi;
	Parameters parameters;
	Errors errors;
	vector<double> likelihood_t;
	EM(const Model & ,const Data &,
	   const double & convergence);
	void iterEM();
	void startXtimes(const int & times);
	void startUntilConvergence();
	void defineparIni(const int & gridSize, const int &times);
};

#endif
