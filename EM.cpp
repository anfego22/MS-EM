#include "EM_Classes.h"

Model::Model(const int &N_, const int &lagsS_, const int &lagsY_,
			 const bool sigma_, const bool beta_):
	N(N_), lagsS(lagsS_), lagsY(lagsY_), sigma(sigma_), beta(beta_){}

Model::Model(const int &N_, const int &lagsS_=0, const int &lagsY_=0,
			 const bool sigma_ = TRUE, const bool beta_ = TRUE):
	N(N_), lagsS(lagsS_), lagsY(lagsY_), sigma(sigma_), beta(beta_){}

Data::Data(const MatrixXd & y,const MatrixXd & x):
	y(y), x(x) {
	T = y.rows();
	m = y.cols();
	mX = x.cols();
}

Data::Data(const MatrixXd & y):
	data(y){
	T = y.rows();
	m = y.cols();
	x.setOnes(T, 0);
	mX = x.cols();
}


void EM::makeXY(const Data &data, const Model & model){
	m = model.lagY;
	Y = data.y.block(m-1, 0, data.T-1-m, data.m-1);
	X = data.x.block(m-1, 0, data.T-1-m, data.mX-1);
	T = Y.rows();
	double M = y.cols()*model.lagsY;
	MatrixXd laggsY(y.rows(), M);
	for (int i = 0; i < model.lagsY; i++)
		laggsY.block(0, i*data.m, T-1, (i+1)*data.m) =
			data.y.block(model.lagsY-i, 0, T-1, m-1);
}

void EM::newBetaS(){
	for (int i == 0; i<= N; i++){
		beta[i] = solve(Data.X*Data.X.cwiseProduct(X))*Data.X*Data.Y*XiS;
	}
}
	
void EM::newBeta(){
	if(model.beta)
		newBetaS();
	else
		newBetaNS();
}

void EM::newSigma(){
	if(model.beta)
		newSigmaS();
	else
		newSigmaNS();
}

void EM::iter(){
	filterProb();
	smoothProb();
	newRho();
	newP();
	newBeta();
	newSigma();
}

EM::EM(const Model &model,const Data &d, const double &c):
	convLevel(0.0), delta(c), Y(){
	setInitialParameters();
	while(delta > convLevel){
		iter();
	}
}


