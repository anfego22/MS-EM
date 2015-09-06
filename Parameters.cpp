#include "EM_Classes.h"
#include "RN_Class.h"

void Transitions::createF(){
	F.setZero(Nm, Nm);
	P.setZero(N, N);
	randomNb rnb;
	P = P.unaryExpr(std::ptr_fun(rnb.sampleUniform));
	VectorXd normTerm(N);
	normTerm = P.colwise().sum();
	for (int i = 0; i <= N - 1; i++){
		P.col(i) = P.col(i)/normTerm(i);
	}
	if (Nm != N){
		MatrixXd Fi(Nm, N);
		Fi.setZero(Nm, N);
		for(unsigned int i = 0, j = 0;
			i < N; i++, j+=N){
			Fi.col(i).segment(j, N) = P.col(i);
		}
		for (unsigned int i = 0; i<Nm; i+=N){
			F.block(0,i,Nm,N) = Fi;
		}
	} else
		F = P;
}

// rho = (A'A)^{-1}A'e
// Since A'A is positive semidefinite matrix we use ldlt method of matrix

void Transitions::createRho(){
	rho.setZero(Nm);
	MatrixXd A(Nm+1, Nm+1);
	MatrixXd I;
	I.setIdentity(Nm, Nm);
	A << (I - F), MatrixXd::setOnes(0,d);
	rho = (A.transpose()*A).ldlt.solve(A.transpos()*I.col(Nm));
}

// updateF needs the smooth probabilities
void Transitions::updateF(const MatrixXd & Xtt, const MatrixXd & Xt1t,
						  const MatrixXd & XiS, const int &T,
						  const int & N){
	for(unsigned int j = 0; j < Nm; ++j){
		for(unsigned int i = 0;  i < Nm; ++i){
			xS = XiS.row(j).tail(T - 1);
			xtt_lag1 = Xitt.row(i).head(T - 1);
			xt1t = Xit1t.row(j).head(T).tail(T - 1);
			xS_lag1 = XiS.row(i).head(T - 1);
			F(j, i) = F(j, i)*((xS.cwiseProduct(xtt_lag1).cwiseQuotient(xt1t)).sum())/(xS_lag1.sum());
		}
	}
}

void Transitions::updateRho(const MatrixXd & XiS){
	rho = Xis.row(0);
}

Transitions::Transitions(const Model &model): N(model.N), Nm(model.Nm){
	createF();
	createRho();
}

void linearParams::createS(){
	if (model.sigma == TRUE)
		int N = model.N;
	else
		int N = 1;
	
	int M = data.Y.ncols();
	MatrixXd Sig;
	randomNb rnb;
	sigma.reserve(N);
	for (int i = 0; i <N; i++){
		Sig.setIdentity(M,M)*=M;
		Sig = Sig.unaryExpr(std::ptr_fun(rnb.varCovRg));
		Sig.triangularView<Eigen::Lower>() = Sig.triangularView<Eigen::Upper>().transpose();
		sigma.push_back(Sig);
		cout << sigma[i] << endl;
	}
}


void linearParams::createB(){
	int N;
	if (model.Beta == TRUE)
		N = model.N;
	else
		N = 1;
	// The concatenated phi coefficients in a single matrix
	// plus the intercept term
	int mM = M*Nm;
	beta.reserve(N);
	MatrixXd Be;
	for (int i = 0; i<N; i++){
		Be.setRandom(M,mM);
		Be.unaryExpr(std::ptr_fun(rnb.sampleCauchy));
		beta.push_back(Be);
	}
}

linearParams::linearParams(const Model &model,
						   const Data &data): data(data), model(model){
	createS();
	createB();
}
