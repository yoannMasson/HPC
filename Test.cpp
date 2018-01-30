#include <iostream>
#include "Analytic.h"
#include "Laasonen.h"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "CrankNicolson.h"
#include <mpi.h>
#include "FTCS.h"


using namespace std;

const double THICNESS = 1.0;
const double TIN = 100.0;
const double TSUR = 300.0;
const double T = 0.5;
const double D = 0.1;
const double DX = 0.05;
const double DT = 0.01;

int main(int argc, char *argv[]) {

	char choice('3');
	double t0,t1;

	MPI_Init(&argc, &argv);

	MPI_Status status;
	int myRank,npes,send[10],recieve[10];
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	ofstream f;

	Analytic analytic(DX,DT,THICNESS,T,D,TSUR,TIN,npes,myRank);//Use to compute errors
	FTCS dufortFrankel(DX,DT,THICNESS,T,D,TSUR,TIN,npes,myRank);
	Laasonen laasonen(DX,DT,THICNESS,T,D,TSUR,TIN,npes,myRank);
	CrankNicolson crankNicholson(DX,DT,THICNESS,T,D,TSUR,TIN,npes,myRank);


	cout << fixed;
	cout << setprecision(2);

	Solver *solution;

	//analytic.computeSolution();

	switch (choice) {

	case '1':
		solution = &dufortFrankel;
		f.open ("DufortFrankel");
		break;

	case '2':
		solution = &laasonen;
		f.open ("Laasonen");
		break;

	case '3':
		solution = &crankNicholson;
		f.open ("CrankNicholson");
		break;

	case '4':
		solution = &analytic;
		f.open ("Analytic");
		break;
	}

	t0 = MPI_Wtime();
	solution->computeSolution();
	t1 = MPI_Wtime();



	if(myRank == 0){
		cout << solution->getComputedSolution();
		cout << setprecision(8);
		cout << t1-t0;

		/*Matrix error = solution->getComputedSolution()-analytic.getComputedSolution();
		cout << "it took "<< t1-t0 << " ms to compute" << endl;
		cout << "ERRORS: "<< endl;
		cout << "one norm: "<< error.one_norm()<< endl;
		cout << "second norm: "<< error.two_norm()<< endl;
		cout << "uniform norm: "<< error.uniform_norm()<< endl;*/
		f << fixed;
		f << setprecision(2);
		f << *solution;
		f.close();
	}


	MPI_Finalize();
	return 0;
}

