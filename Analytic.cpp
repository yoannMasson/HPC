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
	for (int i = 1; i < nRows; i++) {
		x = 0;
		t += dt;
		m[i][0] = Tsur;
		m[i][nCols-1] = Tsur;
		for (int j = startIndex; j < lastIndex; j++) {
			if(j > 0 && j < m.getNcols()-1){
				sum = 0.0;
				x += dx;
				for(int z = 1; z<=50 ; z++){
					sum += exp(-D*(z*PI)*(z*PI)*t)*((1-(pow(-1,z)))/(z*PI))*sin(z*PI*x);
				}
				m[i][j] = Tsur+2*(Tin-Tsur)*sum;
			}
		}

	}

	(*this).computedSolution = m;
	return m;

}



