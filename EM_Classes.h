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

void embed(MatrixXd &,const MatrixXd &,int);

class Model {
public:
	int N, Nm, lagsY, lagsX;
	bool mean, sigma, betaY, betaX, meanCorrected;
	Model(const int &N_ = 2, const int &lagsY_= 0,
		  const int &lagsX_ = 0,const bool &mean_=true,
		  const bool &sigma_ = true,
		  const bool &betaY_ = true, const bool &betaX_ = true,
		  const bool &meanCorrected = false);
	// N is the original number of regimes
	// Nm is N^{m+1} ficticial regimes
	// N = Nm only if meanCorrected = false
};

class Data {
public:
	MatrixXd Y, X;
	int T, k, M;
	Data(const MatrixXd & ,
		 const MatrixXd & );
	Data(const MatrixXd & );
};

class Transitions{
public:
	const Model &model;
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
public:
	const Data  &data;
	const Model &model;
	vector<MatrixXd> betaY, betaX, sigma;
	MatrixXd mu;
	linearParams(const Model &, const Data &);
	void createB();
	void createS();
	void createMu();
	void update();
	void updateBetaS(); // Update if beta has switching
	void updateBetaNS(); // Update if beta has no switching
	void updateSigmaS(); // Update if sigma has switching
	void updateSigmaNS(); // Update if sigma has no switching
};

class Parameters{
public:
	const Model &model;
	const Data &data;
	Transitions rhoF;
	linearParams lin;
	// M -> Number of dependent variables (dimensions of Y)
	Parameters(const Model &, const Data &);
};

class Errors {
public:
	int T, Nm;
	const Model &model;
	const Data &data;
	const Parameters &param;
	MatrixXd eta, Y, Yt1t, X, means;
	Errors(const Parameters &);
	void designMatrix();
	int permFun(const int &, int ,
				const int &);
	void meansF();
	void addMeans(const int &, MatrixXd &);
	void createMuJ(const MatrixXd &, const MatrixXd &,
				   MatrixXd &);
	void createMuJ(const MatrixXd &, const MatrixXd &,
				   const MatrixXd &, const MatrixXd &,
				   MatrixXd &);
	void createMuJ(const int &, const int &,
				   MatrixXd &);
	void createEta();
	
};

class Xi {
public:
	int Nm;
	const Model &model;
	const Data &data;
	const Parameters &param;
	const Errors &errors;
	double logLike;
	MatrixXd Xitt, Xit1t, XiS;
	VectorXd llikel;
	Xi(const Errors &);
	void filterProb();
	void smoothProb();
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
