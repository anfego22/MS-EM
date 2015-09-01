#include "EM_Classes.h"

Transitions::Transitions(const Model &model): N(model.N), Nm(model.Nm){
	F.setZero(Nm, Nm);
}

// updateF needs the smooth probabilities
void updateF(const MatrixXd & Xtt, const MatrixXd & Xt1t,
			 const MatrixXd & XiS, const int &T, const int & N){
	for(unsigned int j = 0; j <= parSelect.N - 1; ++j){
		for(unsigned int i = 0;  i <= parSelect.N - 1; ++i){
			xS = XiS.row(j).tail(T - 1);
			xtt_lag1 = Xitt.row(i).head(T - 1);
			xt1t = Xit1t.row(j).head(T).tail(T - 1);
			xS_lag1 = XiS.row(i).head(T - 1);
			F(j, i) = F(j, i)*((xS.cwiseProduct(xtt_lag1).cwiseQuotient(xt1t)).sum())/(xS_lag1.sum());
		}
	}
}

void updateRho(const MatrixXd & Xtt, const MatrixXd & Xt1t,
			   const MatrixXd & XiS){
	
}
