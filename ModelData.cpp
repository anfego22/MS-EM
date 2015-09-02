#include "EM_Classes.h"

Model::Model(const int &N_, const int &lagsS_, const int &lagsY_,
			 const bool &sigma_, const bool &beta_,
			 const bool &meanCorrected):
	N(N_), lagsY(lagsY_), sigma(sigma_),
	beta(beta_), meanCorrected(meanCorrected){
	if (meanCorrected == TRUE)
		Nm = std::pow(N, lagsY + 1);
}

Model::Model(const int &N_ = 2, const int &lagsY_= 0,
			 const bool &sigma_ = TRUE, const bool &beta_ = TRUE,
			 const bool &meanCorrected = FALSE):
	N(N_), lagsY(lagsY_), sigma(sigma_), beta(beta_),
	meanCorrected(meanCorrected){}

Data::Data(const MatrixXd & y,const MatrixXd & x):
	Y(y){
	MatrixXd c;
	if (y.rows() > y.cols())
		Y.transposeInPlace();
	if (x.rows() > x.cols())
		x.transposeInPlace();
	T = y.rows();
	M = y.cols(); // if M = 1, Univariate model
	// T can't be different from x.rows()
	Mx = x.cols();
	X << c.setOnes(T, 0) << x;
}

Data::Data(const MatrixXd & y):
	Y(y){
	if(y.cols > y.rows)
		y.transposeInPlace()
	T = y.rows();
	M = y.cols();
	x.setOnes(T, 0);
	mX = x.cols();
}
