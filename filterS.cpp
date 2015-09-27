#include "EM_Classes.h"

void Xi::filterProb(){
	Xitt.setZero(Nm, errors.T);
	Xit1t.setZero(Nm, errors.T+1);
	Xit1t.col(0) = param.rhoF.rho;
	llikel.setZero(errors.T);
	
	for (int t = 0; t< errors.T; t++){
		llikel[t] = Xit1t.col(t).dot(errors.eta.col(t));
		Xitt.col(t) = Xit1t.col(t).cwiseProduct(errors.eta.col(t));
		Xitt.col(t) /= llikel(t);
		Xit1t.col(t+1) = param.rhoF.F*Xitt.col(t);
		llikel(t) = log(llikel(t));
	}
	try {
    logLike = llikel.sum();
    if(std::isinf(logLike) | std::isnan(logLike))
      throw runtime_error("likelihood is infinity or is nan");
  } catch (runtime_error err){
		err.what();
  }
}

void Xi::smoothProb(){
	XiS.setZero(Nm, errors.T);
	XiS.col(errors.T - 1) = Xitt.col(errors.T - 1);
	for (int t = errors.T - 2; t != -1; t--){
		XiS.col(t) = Xitt.col(t).cwiseProduct(
			param.rhoF.F.transpose()*
			(XiS.col(t+1).cwiseQuotient(Xit1t.col(t+1)))
			);
	}
}


Xi::Xi(const Errors &err):errors(err), param(err.param),
						  data(err.data), model(err.model){
	Nm = errors.eta.rows();
	filterProb();
	smoothProb();
}

