#include "EM_Classes.h"

Model::Model(const int &N_, const int &lagsY_,
			 const bool &sigma_, const bool &beta_,
			 const bool &meanCorrected):
	N(N_), lagsY(lagsY_), sigma(sigma_), beta(beta_),
	meanCorrected(meanCorrected){
	Nm = std::pow(N, lagsY + 1);
	}
	



Data::Data(const MatrixXd & y,const MatrixXd & x):
	Y(y){
	if (y.rows() < y.cols())
		Y.transposeInPlace();
	T = Y.rows();
	M = Y.cols(); // if M = 1, Univariate model
	// T can't be different from x.rows()
	if (x.rows() > x.cols()){
		mX = x.cols();
		X.setZero(T, mX+1);
		X.col(0) = MatrixXd::Ones(T, 1);
		X.rightCols(mX) = x;
	} else {
		mX = x.rows();
		X.setZero(T, mX+1);
		X.col(0) = MatrixXd::Ones(T, 1);
		X.rightCols(mX) = x.transpose();
	}
}

Data::Data(const MatrixXd & y):
	Y(y){
	if(y.cols() > y.rows())
		Y.transposeInPlace();
	T = Y.rows();
	M = Y.cols();
	X.setOnes(T, 1);
	mX = 1;
}
