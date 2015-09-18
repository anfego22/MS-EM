#include "EM_Classes.h"

Model::Model(const int &N_, const int &lagsY_,
			 const int &lagsX_,
			 const bool &sigma_, const bool &betaY_,
			 const bool &betaX_, const bool &meanCorrected):
	N(N_), lagsY(lagsY_), lagsX(lagsX_),
	sigma(sigma_), betaY(betaY_),
	betaX(betaX_), meanCorrected(meanCorrected){
	Nm = (meanCorrected == true) ? std::pow(N, lagsY + 1):N;
	}

Data::Data(const MatrixXd & y,const MatrixXd & x):
	Y(y){
	if (y.rows() < y.cols())
		Y.transposeInPlace();
	T = Y.rows();
	M = Y.cols(); // if M = 1, Univariate model
	// T can't be different from x.rows()
	if (x.rows() > x.cols()){
		k = x.cols();		   
		X.setZero(T, k);
		X = x;
	} else {
		k = x.rows();
		X.setZero(T, k);
		X = x.transpose();
	}
}

Data::Data(const MatrixXd & y):
	Y(y){
	if(y.cols() > y.rows())
		Y.transposeInPlace();
	T = Y.rows();
	M = Y.cols();
	k = 0;
}

void Data::embed(MatrixXd &Result, const MatrixXd &Y,
				  int m){
	const int M = Y.cols();
	MatrixXd ytm(Y.rows()-m, Y.cols());
	Result.resize(Y.rows()-m, Y.cols()*(m+1));
	for (int i = 0; i<=m; i++){
		ytm = Y.block(i, 0, Y.rows() - m, Y.cols());
		Result.middleCols((m-i)*M, M) = ytm;
	}
}
