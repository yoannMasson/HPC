/*
 * CrankNicholson.cpp
 *
 *  Created on: 31 oct. 2017
 *      Author: yoann
 */

#include "Laasonen.h"
#include "vector.h"
#include "mpi.h"
#include "math.h"

using namespace std;

Laasonen::Laasonen(double dx, double dt, double L, double T, double D, double Tsur, double Tin, int npes, int myRank ):
																																																		Solver (dx, dt, L, T, D, Tsur, Tin,npes,myRank){}


Matrix Laasonen::computeSolution(){

	Matrix m = getComputedSolution();//Matrix containing the result
	int nRows = m.getNrows();
	int nCols = m.getNcols();//We remove the boundary conditions, that are not going to change through time
	double const C = D*dt/(2*dx*dx);


	double *bottomDiagonal = new double[nCols];//The upper diagonal of the triangle matrix A of the system A.F = D
	double *upDiagonal = new double[nCols];//The upper diagonal of the triangle matrix A of the system A.F = D
	double *diagonal = new double[nCols];//The main diagonal of the triangle matrix A of the system A.F = D
	double *resultVector = new double[nCols];//The D vector of the system A.F = D;
	double *f = new double [nCols];//The f vector of the previus iteration


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

	for(int i = 0; i< nCols ; i++){
		diagonal[i] = 1+2*C;
	}

	for(int i = 0; i< nCols ; i++){
		upDiagonal[i] = -C;
	}

	for(int i = 0; i< nCols ; i++){
		bottomDiagonal[i] = -C;
	}

	for(int i = 0 ; i < nCols ; i++){
		resultVector[i] = m[0][i];
	}

	//f is the fist time step at initialization
	for(int i = 0 ; i < nCols; i++){
		f[i] = m[0][i];
	}

	//Calcul
	for (int timeStep = 1; timeStep < nRows; timeStep++) {

		for(int i = 0 ; i < nCols ; i++){
			resultVector[i] = f[i];
		}


		resultVector[0] += C*Tsur;
		resultVector[nCols -1 ] += C*Tsur;
		ThomasAlgorithm_P(nCols-1,bottomDiagonal,diagonal,upDiagonal,f,resultVector);
		f[0] = Tsur;
		f[nCols-1] = Tsur;
		if(myRank == 0){
			//			for(int i = 0 ; i < nCols ; i++)
			//				cout << bottomDiagonal[i]<< " ";
			//			cout << endl;
			//			for(int i = 0 ; i < nCols ; i++)
			//				cout  << diagonal[i]<< " ";
			//			cout << endl;
			//			for(int i = 0 ; i < nCols ; i++)
			//				cout << upDiagonal[i]<< " ";
			//			cout << endl;
			for(int i = 0; i < nCols /2 ;i++){
				f[nCols-i-1] = f[i];
			}
			//			for(int i = 0 ; i < nCols ; i++)
			//				cout  << resultVector[i]<< " ";
			//			cout << endl;
			//			for(int i = 0 ; i < nCols ; i++)
			//				cout << f[i] << " ";
			//			cout << endl;
		}

		for(int i = 0; i < nCols ; i++){
			m[timeStep][i] = f[i];
		}



	}

	(*this).computedSolution = m;
	return m;
}

//Thanks to George Em Karniadakis and Robert M. Kirby for the algorithm, http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.706.6668&rep=rep1&type=pdf
void Laasonen::ThomasAlgorithm_P( int N, double *b,double *a, double *c, double *x, double *q){
	int i;
	int rows_local,local_offset;
	double S[2][2],T[2][2],s1tmp,s2tmp;
	double *l,*d,*y;
	MPI_Status status;
	l = new double[N];
	d = new double[N];
	y = new double[N];
	for(i=0;i<N;i++)
		l[i] = d[i] = y[i] = 0.0;
	S[0][0] = S[1][1] = 1.0;
	S[1][0] = S[0][1] = 0.0;
	rows_local = (int) floor(N/npes);
	local_offset = myRank*rows_local;
	// Form local products of R_k matrices
	if(myRank==0){
		s1tmp = a[local_offset]*S[0][0];
		S[1][0] = S[0][0];
		S[1][1] = S[0][1];
		S[0][1] = a[local_offset]*S[0][1];
		S[0][0] = s1tmp;
		for(i=1;i<rows_local;i++){
			s1tmp = a[i+local_offset]*S[0][0] -
					b[i+local_offset-1]*c[i+local_offset-1]*S[1][0];
			s2tmp = a[i+local_offset]*S[0][1] -
					b[i+local_offset-1]*c[i+local_offset-1]*S[1][1];
			S[1][0] = S[0][0];
			S[1][1] = S[0][1];
			S[0][0] = s1tmp;
			S[0][1] = s2tmp;
		}
	}
	else{
		for(i=0;i<rows_local;i++){
			s1tmp = a[i+local_offset]*S[0][0] -
					b[i+local_offset-1]*c[i+local_offset-1]*S[1][0];
			s2tmp = a[i+local_offset]*S[0][1] -
					b[i+local_offset-1]*c[i+local_offset-1]*S[1][1];
			S[1][0] = S[0][0];
			S[1][1] = S[0][1];
			S[0][0] = s1tmp;
			S[0][1] = s2tmp;
		}
	}
	// Full-recursive doubling algorithm for distribution
	double np(npes);
	for(i=0; i<=log2(np);i++){
		if(myRank+pow(2,i) < npes)
			MPI_Send(S,4,MPI_DOUBLE,int(myRank+pow(2,i)),0,MPI_COMM_WORLD);
		if(myRank-pow(2,i)>=0){
			MPI_Recv(T,4,MPI_DOUBLE,int(myRank-pow(2,i)),0,MPI_COMM_WORLD,&status);
			s1tmp = S[0][0]*T[0][0] + S[0][1]*T[1][0];
			S[0][1] = S[0][0]*T[0][1] + S[0][1]*T[1][1];
			S[0][0] = s1tmp;
			s1tmp = S[1][0]*T[0][0] + S[1][1]*T[1][0];
			S[1][1] = S[1][0]*T[0][1] + S[1][1]*T[1][1];
			S[1][0] = s1tmp;
		}
	}
	//Calculate last d_k first so that it can be distributed,
	//and then do the distribution.
	d[local_offset+rows_local-1] = (S[0][0] + S[0][1])/
			(S[1][0] + S[1][1]);
	if(myRank == 0){
		MPI_Send(&d[local_offset+rows_local-1],1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
	}
	else{
		MPI_Recv(&d[local_offset-1],1,MPI_DOUBLE,myRank-1,0,MPI_COMM_WORLD,&status);
		if(myRank != npes-1)
			MPI_Send(&d[local_offset+rows_local-1],1,MPI_DOUBLE,myRank+1,0,MPI_COMM_WORLD);
	}
	// Compute in parallel the local values of d_k and l_k
	if(myRank == 0){
		l[0] = 0;
		d[0] = a[0];
		for(i=1;i<rows_local-1;i++){
			l[local_offset+i] = b[local_offset+i-1]/
					d[local_offset+i-1];
			d[local_offset+i] = a[local_offset+i] -
					l[local_offset+i]*c[local_offset+i-1];
		}
		l[local_offset+rows_local-1] = b[local_offset+rows_local-2]/
				d[local_offset+rows_local-2];
	}
	else{
		for(i=0;i<rows_local-1;i++){
			l[local_offset+i] = b[local_offset+i-1]/
					d[local_offset+i-1];
			d[local_offset+i] = a[local_offset+i] -
					l[local_offset+i]*c[local_offset+i-1];
		}
		l[local_offset+rows_local-1] = b[local_offset+rows_local-2]/
				d[local_offset+rows_local-2];
	}
	/***************************************************************/
	if(myRank>0)
		d[local_offset-1] = 0;
	// Distribute d_k and l_k to all processes
	double * tmp = new double[N];
	for(i=0;i<N;i++)
		tmp[i] = d[i];
	MPI_Allreduce(tmp,d,N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	for(i=0;i<N;i++)
		tmp[i] = l[i];
	MPI_Allreduce(tmp,l,N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	delete[] tmp;
	if(myRank ==0){
		/* Forward Substitution [L][y] = [q] */
		y[0] = q[0];
		for(i=1;i<N;i++)
			y[i] = q[i] - l[i]*y[i-1];
		/* Backward Substitution [U][x] = [y] */
		x[N-1] = y[N-1]/d[N-1];
		for(i=N-2;i>=0;i--)
			x[i] = (y[i] - c[i]*x[i+1])/d[i];
	}
	delete[] l;
	delete[] y;
	delete[] d;

	return;
}

Laasonen::~Laasonen() {
	// TODO Auto-generated destructor stub
}

