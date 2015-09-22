#include "EM_Classes.h"
#include "RN_Class.h"

using Eigen::MatrixXd;


/* This Functions takes an integer number with range
   0:(N^{m+1}-1) and other integer j with range
   0:m. It returns the state of mu_{s*_{t}} that correspond
   to mu_{s_{t-j}}. */

int Errors:: permFun(const int &N, int i,
					 const int &j){
	i = i/pow(N, j);
	i %= N;
	return i;
}

// This function makes a (N^{m+1}Mx1) elements 
void Errors::meansF(){
	if(model.meanCorrected){
		means.setZero(model.Nm, data.M);
		int sel0, sel;
		MatrixXd Phi, muJ;
		Phi.setZero(data.M, data.M*(model.lagsY +1));
		muJ.setZero(data.M*(model.lagsY +1), 1);
		for (int i = 0; i < model.Nm; i++){
			sel0 = permFun(model.N, i, 0);
			for (int j = 0; j < model.lagsY; j++){
				sel = permFun(model.N, i, j);
				muJ.middleRows(j*data.M, data.M) =
					param.lin.mu.row(sel).transpose();
			}
			Phi << MatrixXd::Identity(data.M, data.M),
				param.lin.betaY[sel0];
			means.row(i) = (Phi*muJ).transpose();
		}
	} else 
		means = param.lin.mu;
	
}

void Errors::designMatrix(){
	embed(Y, param.data.Y, param.model.lagsY);  
	Yt1t = Y.rightCols(Y.cols()-data.M);
	Y = Y.leftCols(data.M);
	T = Y.rows();
	embed(X, param.data.X, param.model.lagsX);
	if( T < X.rows()){
		X = X.topRows(T);
	}
	if( T > X.rows()){
		cerr << "ERROR: X_ROWS < Y_ROWS" << endl; 
	}
}

void Errors::createEta(){
	MatrixXd MuJ;
	for (int j = 0, i = 0; j < model.Nm; j++, i++){
		if(i >= model.N)
			i = 0;
		if (model.betaX){
			MuJ = Yt1t*param.lin.betaY[i].transpose()
				+ X*param.lin.betaX[i].transpose();
		} else {
			MuJ = Yt1t*param.lin.betaY[i].transpose()
				+ X*param.lin.betaX[0].transpose();
		}
		if (model.mean){
			for(int h = 0; h < MuJ.rows(); h++){
				MuJ.row(h) += means.row(j);
			}
		} else{
			for(int h = 0; h < MuJ.rows(); h++){
				MuJ.row(h) += means;
			}
		}
		if (model.sigma){
			eta.row(j) = multivariateNormal_density(Y, MuJ,
													param.lin.sigma[i]);
		} else{
			eta.row(j) = multivariateNormal_density(Y, MuJ,
													param.lin.sigma[0]);
		}
		
	}
}



Errors::Errors(const Parameters &param):param(param),
										model(param.model),
										data(param.data){
	
	
	designMatrix();
	meansF();
	eta.setZero(model.Nm, T);		   
	// 
	createEta();
}

