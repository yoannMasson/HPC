


#include <mpi.h>
int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	int n = 100;
	int myRank;
	int a[100];
	double b,t0,t1;
	
	MPI_Request request[4000];
	MPI_Status status;
	
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	
	for(int i =  0; i < 10 ; i++){
		a[i] = i+0.5;
	}
	
	t0 = MPI_Wtime();
	for(int i = 0; i < 1000; i++ ){
		if(myRank == 0){
			MPI_Send(&a,n,MPI_DOUBLE,1,i,MPI_COMM_WORLD); 
			MPI_Recv(&b,n,MPI_DOUBLE,1,i,MPI_COMM_WORLD,&status); 
		}else if (myRank == 1){
			MPI_Recv(&b,n,MPI_DOUBLE,0,i,MPI_COMM_WORLD,&status);  
			MPI_Send(&a,n,MPI_DOUBLE,0,i,MPI_COMM_WORLD); 
			
		}
	}
	
	
	t1 = MPI_Wtime();
	
	if(myRank == 0){
		printf ("Hey this is procees %d, sending 1000 messages of size %d took %lf secondes",myRank, n, t1-t0);		
	}	
	MPI_Finalize();
	 
}
