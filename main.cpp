#include "EM_Classes.h"
#include "RN_Class.h"

using namespace std;
using Eigen::MatrixXd;


int main(){
	Model myModel;
	cout << "This is N^{m+1}" << endl;
	cout << myModel.Nm << endl;

	Eigen::Matrix<double, 2, 10> Y;
	Y.setRandom();

	Data myData(Y);
	cout << "This is Y" << endl;
	cout << Y << endl;
	cout << "This is Y from the Data class" << endl;
	cout << myData.Y << endl;

	cout << "And X from the Data class" << endl;
	cout << myData.X << endl;

	Parameters myParam(myModel, myData);
	cout << "Take a look at sigma" << endl;
	for (int i = 0; i< myModel.N; i++){
		cout << myParam.betaSigma.sigma[i] << endl;
		cout << endl;
	}
	
	cout << "Take a look at Beta" << endl;
	for (int i = 0; i< myModel.N; i++){
		cout << myParam.betaSigma.beta[i] << endl;
		cout << endl;
	}

	
	
}
