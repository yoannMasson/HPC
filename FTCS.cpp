/*
 * DufortFrankel.cpp
 *
 *  Created on: 27 oct. 2017
 *      Author: yoann
 */

#include "FTCS.h"

#include <cmath>
#include "mpi.h"

using namespace std;
FTCS::FTCS(double dx, double dt, double L, double T, double D, double Tsur, double Tin,int npes, int myRank ):
														Solver (dx, dt, L, T, D, Tsur, Tin,npes,myRank){}

Matrix FTCS::computeSolution(){

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
		offset = ((double)1/npes)*(m.getNcols()-2);
		startIndex = offset*myRank+1 ;
		if(myRank == npes -1 ){
			lastIndex = m.getNcols()-2;
		}else{
			lastIndex = startIndex + offset-1;
		}
	}else{//one processor for one space step
		startIndex = myRank+1;
		lastIndex = myRank+1;
	}

	std::cout << "From processor " << myRank << " I start at " << startIndex << "and finish at " << lastIndex << "\n";


	double r = D*dt/(dx*dx);

	//CALCUL
	MPI_Request request;
	MPI_Status status;
	int tag;
	double previous,following;

	for (int i = 0; i < m.getNrows()-1; i++) {
		if(i == 0){//Special case where don't need result from previous line
			for(int j = startIndex; j <= lastIndex;j++){
				m[i+1][j] = r*(m[i][j+1]-2*m[i][j]+m[i][j-1])+m[i][j];
				if(myRank != 0){//Send result back to proc 0
					tag = (i+1)*m.getNcols()+j;
					MPI_Isend(&m[i+1][j],1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&request);
					if(j==startIndex){
						MPI_Isend(&m[i+1][j],1,MPI_DOUBLE,myRank-1,tag,MPI_COMM_WORLD,&request);
					}
				}
				if(j == lastIndex && myRank != npes -1){
					tag = (i+1)*m.getNcols()+j;
					MPI_Isend(&m[i+1][j],1,MPI_DOUBLE,myRank+1,tag,MPI_COMM_WORLD,&request);
				}
			}


		}else{//When the current j is at the "processor boundaries" we need to send and recieve

			for (int j = startIndex; j <= lastIndex ; j++) {

				previous = m[i][j-1];
				following = m[i][j+1];

				if(j == startIndex && myRank != 0){
					tag = i*m.getNcols()+j-1;
					MPI_Recv(&previous,1,MPI_DOUBLE,myRank-1,tag,MPI_COMM_WORLD,&status);

				}
				if(j == lastIndex && myRank != npes-1){
					tag = i*m.getNcols()+j+1;
					MPI_Recv(&following,1,MPI_DOUBLE,myRank+1,tag,MPI_COMM_WORLD,&status);

				}
				m[i+1][j] = r*(following-2*m[i][j]+previous)+m[i][j];
				tag = (i+1)*m.getNcols()+j;
				if(myRank != 0){//Send result back to proc 0
					MPI_Isend(&m[i+1][j],1,MPI_DOUBLE,0,tag,MPI_COMM_WORLD,&request);
				}
				if(j == startIndex && myRank != 0){
					MPI_Isend(&m[i+1][j],1,MPI_DOUBLE,myRank-1,tag,MPI_COMM_WORLD,&request);
				}
				if(j == lastIndex && myRank != npes -1){
					MPI_Isend(&m[i+1][j],1,MPI_DOUBLE,myRank+1,tag,MPI_COMM_WORLD,&request);
				}
				//	cout << m[i][j+1]<<endl;
			}



		}

	}

	if(myRank == 0){//Waiting for results
		int waitStartIndex,waitLastIndex;
		double result;
		for(int rank = 1; rank < npes ; rank++){
			waitStartIndex = rank*offset+1;
			if( rank == npes -1){
				waitLastIndex = m.getNcols()-2;
			}else{
				waitLastIndex = waitStartIndex + offset - 1;
			}
			//cout << "Waiting for proc number " << i << " to give me result from " << waitStartIndex << " to " << waitLastIndex << "\n";
			for(int i = 1 ; i < m.getNrows(); i ++ ){
				for( int j = waitStartIndex; j <= waitLastIndex ; j ++){
					tag = j+i*m.getNcols();
					cout << "Waiting for i:"<<i<<" j: "<<j<<" from rank: "<<rank << endl;
					MPI_Recv(&result,1,MPI_DOUBLE,rank,tag,MPI_COMM_WORLD,&status);
					m[i][j] = result;
				}

			}

		}
	}

	(*this).computedSolution = m;
	return m;
}

FTCS::~FTCS() {
	// TODO Auto-generated destructor stub
}

