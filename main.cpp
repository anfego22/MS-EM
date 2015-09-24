#include "EM_Classes.h"
#include "RN_Class.h"

using namespace std;
using Eigen::MatrixXd;


int main(int argc, char *argv[]){
	vector<int> parInt = {2, 0, 0};
	//                      mean, sigma,betaY,betaX,meanCorrec
	vector<bool> parBool = {false, false, true, true, true};
	
	if ( argc > 1){
		for (int i = 0; i < argc-1; i++){
			if (i <= 2)
				parInt[i] = atoi(argv[i+1]);
			else{
				if (atoi(argv[i+1]) == 1)
					parBool[i-3] = true;
			else 
				parBool[i-3] = false;
			}
		}
	}
	
    //                             
	cout << "The call model is call with" << endl;
	cout << parInt[0] << parInt[1] << parInt[2]<< endl;
	cout << "and the switching" << endl;
	cout << parBool[0] << parBool[1] << parBool[2]  << parBool[3] << parBool[4] << endl;
	
	Model myModel(parInt[0], parInt[1], parInt[2],
				  parBool[0], parBool[1], parBool[2],
				  parBool[3], parBool[4]);
	
	cout << "This is N^{m+1}" << endl;
	cout << myModel.Nm << endl;
	Eigen::Matrix<double, 2, 10> Y;
	Eigen::Matrix<double, 2, 10> X;
	Y.setRandom();
	X.setRandom();
	Data myData(Y, X);
	cout << "This is Y" << endl;
	cout << Y << endl;
	cout << "This is X" << endl;
	cout << X << endl;
	cout << "This is Y from the Data class" << endl;
	cout << myData.Y << endl;
	cout << "And X from the Data class" << endl;
	cout << myData.X << endl;

	Parameters myParam(myModel, myData);
	cout << "Take a look at rho" << endl;
	cout << myParam.rhoF.rho << endl;
	cout << myParam.rhoF.rho.sum() << endl;
	cout << "This is F" << endl;
	cout << myParam.rhoF.F << endl << endl;
	cout << "Take a look at mu" << endl;
	cout << myParam.lin.mu << endl;
	cout << "Take a look at sigma" << endl;
	for (int i = 0; i< myParam.lin.sigma.size(); i++){
		cout << myParam.lin.sigma[i] << endl;
		cout << endl;
	}
	cout << "Take a look at Beta" << endl;
	for (int i = 0; i< myParam.lin.betaY.size(); i++){
		cout << myParam.lin.betaY[i] << endl;
		cout << endl;
	}
	cout << "Take a look at BetaX" << endl;
	for (int i = 0; i< myParam.lin.betaX.size(); i++){
		cout << myParam.lin.betaX[i] << endl;
		cout << endl;
	}
	cout << "Take a look at the design matrix" << endl;
	cout << "Y_t\t" << "Y_{t-1}\t\t" << endl;
	
	MatrixXd resY, resX, MuJ;
	embed(resY, myData.Y, parInt[1]);
	cout << resY << endl;
	if (myData.k){
		embed(resX, myData.X, parInt[2]);
		cout << "What happend if we embed a matrix with 0" << endl;
		cout << resX << endl;
	} else
		cout << "There is no X" << endl;

	Errors err(myParam);
	cout << "This is eta" << endl;
	cout << err.eta << endl;

	Xi xi(err);
	cout << "This is XiS" << endl;
	cout << xi.XiS << endl;
	
}

