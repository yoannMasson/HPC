#include <mpi.h>
int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	
	int myRank,a[10],b[10];
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	
	if(myRank == 0){
		for(int i =  0; i < 10 ; i++){
			a[i] = i;
			b[i] = i*10;
		}
	}else{
		for(int i =  0; i < 10 ; i++){
			a[i] = i+5;
			b[i] = i+5;
		}	
	}
	MPI_Status status;
	
	
	if(myRank == 0){
		MPI_Send(a,10,MPI_INT,1,1,MPI_COMM_WORLD); 
		MPI_Send(b,10,MPI_INT,1,2,MPI_COMM_WORLD); 
	}else if (myRank == 1){
		MPI_Recv(a,10,MPI_INT,0,1,MPI_COMM_WORLD,&status); 
		MPI_Recv(b,10,MPI_INT,0,2,MPI_COMM_WORLD,&status); 
	}
	
	if(myRank == 1 ){
		for(int i =  0; i < 10 ; i++){
			printf ("a[%d] = %d \n",i,a[i]);
		}
		for(int i = 0 ; i < 10 ; i++){
			printf ("b[%d] = %d \n",i,b[i]);
		}
	}	
	MPI_Finalize();
	 
}
