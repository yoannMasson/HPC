/*
 * DufortFrankel.cpp
 *
 *  Created on: 27 oct. 2017
 *      Author: yoann
 */

#include "DufortFrankel.h"

#include <cmath>



DufortFrankel::DufortFrankel(double dx, double dt, double L, double T, double D, double Tsur, double Tin,int npes, int myRank ):
Solver (dx, dt, L, T, D, Tsur, Tin,npes,myRank){}

Matrix DufortFrankel::computeSolution(){

	Matrix m = getComputedSolution();
	int nRows = m.getNrows();
	int nCols = m.getNcols();

	//initialisation
	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < nCols ; j++) {
			if(j == 0 || j == nCols-1){
				m[i][j] = Tsur;
			}else{
				m[i][j] = Tin;
			}
		}
	}


	double r = D*dt/(dx*dx);
	//Calcul
	for (int i = 0; i < nRows-1; i++) {
		if(i == 0){//Special case where m[i-1] makes not sense
			for (int j = 1; j < nCols-1; j++) {
				m[i+1][j] = (Tin +2*r*(m[i][j+1]-Tin+m[i][j-1]))/(1+2*r);
			}
		}else{
			for (int j = 1; j < nCols-1; j++) {
				m[i+1][j] = (m[i-1][j] +2*r*(m[i][j+1]-m[i-1][j]+m[i][j-1]))/(1+2*r);
			}
		}

	}
	(*this).computedSolution = m;
	return m;
}

DufortFrankel::~DufortFrankel() {
	// TODO Auto-generated destructor stub
}

