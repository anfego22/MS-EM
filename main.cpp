#include "EM_Classes.h"
#include "RN_Class.h"

using namespace std;
using Eigen::MatrixXd;


int main(){
	int lagsY = 3;
	int lagsX = 0;
	Model myModel(2, lagsY, lagsX, true, true, true, false);
	cout << "This is N^{m+1}" << endl;
	cout << myModel.Nm << endl;
	Eigen::Matrix<double, 3, 10> Y;
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
	embed(resY, myData.Y, lagsY);
	cout << resY << endl;
	embed(resX, myData.X, lagsX);
	cout << "What happend if we embed a matrix with 0" << endl;
	cout << resX << endl;

	MatrixXd &SigmaJ = myParam.lin.sigma[0];
	MatrixXd &beta1 = myParam.lin.betaY[0];
	MatrixXd etaRow;
	MuJ = resY.rightCols(resY.cols() - myData.Y.cols())*beta1.transpose(); 
	cout << "MuJ*Yt1t" <<endl;
	cout << MuJ << endl;
	etaRow.setZero(1, resY.rows());
	multivariateNormal_density(etaRow,
							   resY.leftCols(myData.Y.cols()),
							   MuJ, SigmaJ);
	cout << "This is the first row of eta" << endl;
	cout << etaRow << endl;
}

