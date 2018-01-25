/*
 * CrankNicholson.cpp
 *
 *  Created on: 31 oct. 2017
 *      Author: yoann
 */

#include "CrankNicolson.h"
#include "vector.h"

using namespace std;

CrankNicolson::CrankNicolson(double dx, double dt, double L, double T, double D, double Tsur, double Tin, int npes, int myRank ):
								Solver (dx, dt, L, T, D, Tsur, Tin,npes,myRank){}


Matrix CrankNicolson::computeSolution(){

	Matrix m = getComputedSolution();//Matrix containing the result
	int nRows = m.getNrows();
	int nCols = m.getNcols()-2;//We remove the boundary conditions, that are not going to change through time
	double const C = D*dt/(2*dx*dx);


	Vector bottomDiagonal(nCols-1);//The upper diagonal of the triangle matrix A of the system A.F = D
	Vector upDiagonal(nCols-1);//The upper diagonal of the triangle matrix A of the system A.F = D
	Vector diagonal(nCols);//The main diagonal of the triangle matrix A of the system A.F = D
	Vector resultVector(nCols);//The D vector of the system A.F = D;
	Vector f(nCols);//The f vector of the previus iteration


	//initialisation
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols+2 ; j++) {
			if(j == 0 || j == nCols+1){
				m[i][j] = Tsur;
			}else{
				m[i][j] = Tin;
			}
		}
	}

	for(int i = 0; i< nCols ; i++){
		diagonal[i] = 1+2*C;
	}

	for(int i = 0; i< nCols-1 ; i++){
		upDiagonal[i] = -C;
	}

	for(int i = 0; i< nCols-1 ; i++){
		bottomDiagonal[i] = -C;
	}

	//f is the fist time step at initialization
	for(int j = 0 ; j < nCols; j++){
		f[j] = Tin;
	}



	//Calcul
	for (int timeStep = 1; timeStep < nRows; timeStep++) {

		//Filling resultVector with init value
		resultVector[0] =  C*f[1]+(1-2*C)*f[0]+C*Tsur;
		resultVector[nCols - 1] =  C*Tsur+(1-2*C)*f[nCols - 1]+C*f[nCols - 2];
		resultVector[0] += C*Tsur;
		resultVector[nCols-1] += C*Tsur;
		for(int i = 1 ; i < nCols-1 ; i++){
			resultVector[i] = C*f[i+1]+(1-2*C)*f[i]+C*f[i-1];
		}

		resolveOneStep(bottomDiagonal,diagonal,upDiagonal,resultVector,f);

		for(int i = 0; i < nCols ; i++){
			m[timeStep][i+1] = f[i];
		}
	}

	(*this).computedSolution = m;
	return m;
}

Vector CrankNicolson::resolveOneStep(Vector bottomDiagonal,
		Vector diagonal,
		Vector upDiagonal,
		Vector resultDiagonal,
		Vector &f){

	//Initialization
	int nCols = this->getComputedSolution().getNcols() - 2;
	Vector tempUp(nCols-1);
	Vector tempRes(nCols);

	tempUp[0] = upDiagonal[0]/diagonal[0];
	tempRes[0] = resultDiagonal[0]/diagonal[0];

	for(int i= 1 ; i < nCols-1 ; i++){
		tempUp[i] = upDiagonal[i]/(diagonal[i]-bottomDiagonal[i]*tempUp[i-1]);
	}

	tempRes[0] = resultDiagonal[0]/diagonal[0];
	for(int i= 1 ; i < nCols-1 ; i++){
		tempRes[i] = (resultDiagonal[i]-bottomDiagonal[i]*tempRes[i-1])/(diagonal[i]-bottomDiagonal[i]*tempUp[i-1]);
	}
	tempRes[nCols-1] = (resultDiagonal[nCols-1]-bottomDiagonal[nCols-2]*tempRes[nCols-2])/(diagonal[nCols-1]-bottomDiagonal[nCols-2]*tempUp[nCols-2]);


	//reverse resolving the system

	f[nCols-1] = tempRes[nCols-1];
	for(int i = nCols-2;i >= 0 ; i--){
		f[i] = tempRes[i]-tempUp[i]*f[i+1];
	}

	//	cout << "bottomDiagonal: " << bottomDiagonal;
	//	cout << "diagonal: " << diagonal;
	//	cout << "upDiagonal: " << upDiagonal;
	//	cout << "resultDiagonal: " << resultDiagonal;
	//	cout << "f: " << f;

	return f;

}

CrankNicolson::~CrankNicolson() {
	// TODO Auto-generated destructor stub
}

