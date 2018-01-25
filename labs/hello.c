#include <mpi.h>

int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	
	int npes, myRank;
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	
	printf("Hey, I am rank number: %d/%d\n",(myRank+1),npes);
	
	MPI_Finalize();
	 
}
