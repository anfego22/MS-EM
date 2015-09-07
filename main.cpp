#include "EM_Classes.h"
#include "RN_Class.h"

using namespace std;
using Eigen::MatrixXd;


int main(){
	Model myModel(2, 0, false, false, false);
	cout << "This is N^{m+1}" << endl;
	cout << myModel.Nm << endl;
	
	Eigen::Matrix<double, 2, 10> Y, X;
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
	
	cout << "Take a look at sigma" << endl;
	for (int i = 0; i< myParam.betaSigma.sigma.size(); i++){
		cout << myParam.betaSigma.sigma[i] << endl;
		cout << endl;
	}
	
	cout << "Take a look at Beta" << endl;
	for (int i = 0; i< myParam.betaSigma.beta.size(); i++){
		cout << myParam.betaSigma.beta[i] << endl;
		cout << endl;
	}
	
}

