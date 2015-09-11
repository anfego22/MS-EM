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
		k = x.cols()+1;		   
		X.setZero(T, k);
		X.col(0) = MatrixXd::Ones(T, 1);
		X.rightCols(k-1) = x;
	} else {
		k = x.rows()+1;
		X.setZero(T, k);
		X.col(0) = MatrixXd::Ones(T, 1);
		X.rightCols(k-1) = x.transpose();
	}
}

Data::Data(const MatrixXd & y):
	Y(y){
	if(y.cols() > y.rows())
		Y.transposeInPlace();
	T = Y.rows();
	M = Y.cols();
	k = 1;
	X.setOnes(T, k);
}

void Data::embed(MatrixXd *Result, const MatrixXd &Y,
				  int m){
	MatrixXd ytm(Y.rows()-m, Y.cols());
	Result->resize(Y.rows()-m, Y.cols()*(m+1));
	for (int i = 0; i<=m; i++){
		ytm = Y.block(i, 0, Y.rows() - m, Y.cols());
		Result->middleCols((m-i)*Y.cols(), Y.cols()) = ytm;
	}
}
