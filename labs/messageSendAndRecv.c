#include <mpi.h>
#include <time.h>
int main(int argc, char *argv[]){

	MPI_Status status;
	
	clock_t t;
	int myRank,size;
	MPI_Init(&argc, &argv);
	size = 10000000;
	int a [size],b[size];
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	
	for(int i =  0; i < size ; i++){
		a[i] = i;
	
	}
	
	
   	t = clock();
    	
    	
	if(myRank == 0){
		
		MPI_Sendrecv(a,size,MPI_INT,1,0,b,size,MPI_INT,1,0,MPI_COMM_WORLD,&status); 
		 
		
	}else if (myRank == 1){
		
		MPI_Sendrecv(a,size,MPI_INT,0,0,b,size,MPI_INT,0,0,MPI_COMM_WORLD,&status); 
		 
	}
	printf("lol\n");
	
	t = clock() - t;
    	double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
 
    	printf("Send And Recv took %f seconds to execute \n", time_taken);
		
	MPI_Finalize();
	return 0;
	 
}
