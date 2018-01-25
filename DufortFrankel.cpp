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

	//INDEX ATTRIBUTION
	int startIndex,lastIndex,offset ;

	if(npes < m.getNcols()){//Dividing work between proc
		offset = ((double)1/npes)*m.getNcols();
		startIndex = offset*myRank ;
		if(myRank == npes -1 ){
			lastIndex = m.getNcols()-1;
		}else{
			lastIndex = startIndex + offset-1;
		}
	}else{//one processor for one space step
		startIndex = myRank;
		lastIndex = myRank;
	}

	std::cout << "From processor " << myRank << " I start at " << startIndex << "and finish at " << lastIndex << "\n";


	double r = D*dt/(dx*dx);
	//CALCUL
	for (int i = 0; i < nRows-1; i++) {
		if(i == 0){//Special case where m[i-1][j] makes not sense
			for (int j = startIndex; j <= lastIndex ; j++) {
				if(j == 0){
					if(startIndex == lastIndex){
						//MPI_ISEND();
					}
				}else{
					if(j == startIndex){
						//MPI_RECIEVE()
					}
					m[i+1][j] = (Tin +2*r*(m[i][j+1]-Tin+m[i][j-1]))/(1+2*r);
					if(j == lastIndex){
						//MPI_ISEND();
					}
				}
			}
		}else{

			for (int j = startIndex; j <= lastIndex ; j++) {
				if(j == 0){
					if(startIndex == lastIndex){
						//MPI_ISEND();
					}
				}else{
					if(j == startIndex){
						//MPI_RECIEVE()
					}
					m[i+1][j] = (m[i-1][j] +2*r*(m[i][j+1]-m[i-1][j]+m[i][j-1]))/(1+2*r);
					if(j == lastIndex){
						//MPI_ISEND();
					}
				}



			}

		}
		(*this).computedSolution = m;
		return m;
	}

	DufortFrankel::~DufortFrankel() {
		// TODO Auto-generated destructor stub
	}

