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


	int myRank,npes;
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	ofstream f;

	FTCS ftcs(DX,DT,THICNESS,T,D,TSUR,TIN,npes,myRank);
	Laasonen laasonen(DX,DT,THICNESS,T,D,TSUR,TIN,npes,myRank);
	CrankNicolson crankNicholson(DX,DT,THICNESS,T,D,TSUR,TIN,npes,myRank);

	//Compute FTCS

	t0 = MPI_Wtime();
	ftcs.computeSolution();
	t1 = MPI_Wtime();
	//
	//	if(myRank == 0){
	//		cout << setprecision(4);
	//		cout << "FTCS Result took " << (t1-t0)*1000 << "ms to compute with "<<npes<<" processors"<<endl;
	//		cout << setprecision(2);
	//		cout << fixed;
	//		cout << ftcs;
	//	}
	//
	if ((npes & (npes - 1)) == 0 && npes > 1 ){
		//		//Compute Crank-Nicholson
		//		t0 = MPI_Wtime();
		//		crankNicholson.computeSolution();
		//		t1 = MPI_Wtime();
		//
		//		if(myRank == 0){
		//			cout << "CrankNicholson Result took " << (t1-t0)*1000 << "ms to compute with "<<npes<<" processors"<<endl;
		//			cout << setprecision(2);
		//			cout << fixed;
		//			cout << crankNicholson << endl;;
		//		}

		//Compute Laasonen
		t0 = MPI_Wtime();
		laasonen.computeSolution();
		t1 = MPI_Wtime();

		if(myRank == 0){
			cout << "Laasonen Result took " << (t1-t0)*1000 << "ms to compute with "<<npes<<" processors"<<endl;
			cout << setprecision(2);
			cout << fixed;
			cout << laasonen << endl;
		}

	}else{
		if(myRank==0)
			cout << " SORRY THE IMPLEMENTED ALGORITHM FOR LAASONEN AND CRANKNICHOLSON ONLY WORK WITH NUMBER OF PROCESSOR AS A POWER OF 2"<<endl;
	}





	MPI_Finalize();
	return 0;
}

