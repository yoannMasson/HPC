#include <stdio.h>

#include <mpi.h>
int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	int size = 20000000;
	int localSum = 0;
	int globalSum = 0;
	int npes,myRank,a[20000000],b[20000000],upperBound, underBound;
	double t0,t1;
	MPI_Request request;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	
	
	for(int i =  0; i < size ; i++){//Populating vectors
		a[i] = 1;
	}
	
	for(int i =  0; i < size ; i++){//Populating vectors
		b[i] = 2;
	}
	
	
	underBound = myRank*size/npes;
	upperBound = underBound + size/npes;
	
	
	t0 = MPI_Wtime();
	
	
	for(int i = underBound ; i < upperBound ; i ++ ){
		localSum += a[i]*b[i];
	}
	
	MPI_Reduce(&localSum,&globalSum,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	if(myRank == 0){
		printf("Hey, this is process %d, the sum is %d, it took %f secondes to compute\n",myRank,globalSum,t1-t0);
	}
	MPI_Finalize();
	 
}
