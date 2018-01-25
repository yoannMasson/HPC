#include "Analytic.h"
#include <cmath>


//Constructor of analytic object, calls the upper class constructor: Solver
Analytic::Analytic(double dx, double dt, double L, double T, double D, double Tsur, double Tin, int npes, int myRank ):
Solver (dx, dt, L, T, D, Tsur, Tin,npes,myRank){}

//Computes and returns the solution, using the analytic equation
Matrix Analytic::computeSolution(){

	const double PI = 3.14159265;
	Matrix m = getComputedSolution();
	int nRows = m.getNrows();
	int nCols = m.getNcols();
	double sum (0.0),t(0.0),x(0.0);
	for (int i = 0; i < nCols ; i++) {
		if(i == 0 ||i == nCols-1){
			m[0][i] = Tsur;
		}else{
			m[0][i] = Tin;
		}
	}


	for (int i = 1; i < nRows; i++) {
		x = 0;
		t += dt;
		m[i][0] = Tsur;
		m[i][nCols-1] = Tsur;
		for (int j = 1; j < m.getNcols()-1; j++) {

			sum = 0.0;
			x += dx;
			for(int z = 1; z<=50 ; z++){
				sum += exp(-D*(z*PI)*(z*PI)*t)*((1-(pow(-1,z)))/(z*PI))*sin(z*PI*x);
			}
			m[i][j] = Tsur+2*(Tin-Tsur)*sum;

		}

	}

	(*this).computedSolution = m;
	return m;

}



