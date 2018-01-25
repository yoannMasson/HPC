#include <stdio.h>

#include <mpi.h>
int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);

	int npes,myRank;
	double localValue, globalValue, startingIntervalle, endingIntervalle;
	double t0,t1;
	MPI_Request request;
	MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	
	
	t0 = MPI_Wtime();
	
	startingIntervalle = (double)myRank / (double)npes ;	
	localValue =  4.0/(1+startingIntervalle*startingIntervalle); 
	printf("Hey, this is process %d,localValue is: %lf\n",myRank,localValue);
	MPI_Reduce(&localValue,&globalValue,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
	t1 = MPI_Wtime();

	if(myRank == 0){
		globalValue = globalValue/npes;
		printf("Hey, this is process %d, the sum is %lf, it took %lf secondes to compute\n",myRank,globalValue,t1-t0);
	}
	
	MPI_Finalize();
	 
}

