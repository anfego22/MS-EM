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
	if (model.meanCorrected){
		int i = 0,jj = 0,j = 0;
		while(jj < Nm){
		F.col(jj).segment(j, N) = P.col(i);
		if (i >= N-1)
			i = 0;
		else
			i++;
		if (j < Nm-N)
			j +=N;
		else
			j = 0;
		jj++;
	}
	} else
		F = P;
}

// rho = (A'A)^{-1}A'e
// dim -> 
// Since A'A is positive semidefinite matrix we use ldlt method of matrix

void Transitions::createRho(){
	rho.setZero(Nm);
	MatrixXd A(Nm+1, Nm);
	MatrixXd I;
	I.setIdentity(Nm+1, Nm+1);
	A.setZero(Nm+1, Nm);
	A.topRows(Nm) = MatrixXd::Identity(Nm, Nm) - F;
	A.row(Nm) = VectorXd::Ones(Nm);
	rho = (A.transpose()*A).ldlt().solve(A.transpose()*I.col(Nm));
}

// updateF needs the smooth probabilities
void Transitions::updateF(const MatrixXd & Xitt, const MatrixXd & Xit1t, const MatrixXd & XiS){
	int T = XiS.cols();
	MatrixXd xS, xtt_lag1, xS_lag1, xt1t;
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
	rho = XiS.row(0);
}



Transitions::Transitions(const Model &model): model(model){
	N = model.N;
	Nm = model.Nm;
	createF();
	createRho();
}

void linearParams::createS(){
	int N;

	if (model.sigma)
		N = model.N;
	else
		N = 1;

	int M = data.Y.cols();
	MatrixXd Sig;
	randomNb rnb;
	sigma.reserve(N);

	for (int i = 0; i <N; i++){
		Sig.setIdentity(M,M)*=M;
		Sig = Sig.unaryExpr(std::ptr_fun(rnb.varCovRg));
		Sig.triangularView<Eigen::Lower>() = Sig.triangularView<Eigen::Upper>().transpose();
		sigma.push_back(Sig);
	}
	
}

void linearParams::createB(){
	
	int Ny, Nx;
	Ny = (model.betaY == true)? model.N:1;
	Nx = (model.betaX == true)? model.N:1;
	// The concatenated phi coefficients in a single matrix
	// plus the intercept term
	int betaLags = data.M*(model.lagsY);
	betaY.reserve(Ny);
	betaX.reserve(Nx);
	MatrixXd beY, beX;		
	randomNb rnb;
	
	for (int i = 0; i<Ny; i++){
		beY.setRandom(data.M,betaLags);
		beY.unaryExpr(std::ptr_fun(rnb.sampleCauchy));
		betaY.push_back(beY);
	}
	if(data.k){
		for (int i = 0; i<Nx; i++){
			beX.setRandom(data.M, data.k);
			beX.unaryExpr(std::ptr_fun(rnb.sampleCauchy));
			betaX.push_back(beX);
		}
		
	}

}

void linearParams::createMu(){
	MatrixXd mu;
	mu.setRandom(model.N, data.M);
}
	
linearParams::linearParams(const Model &model,
						   const Data &data): data(data), model(model){
	createS();
	createB();
}

Parameters::Parameters(const Model &model,
					   const Data &data): rhoF(model), lin(model, data), model(model), data(data){}
