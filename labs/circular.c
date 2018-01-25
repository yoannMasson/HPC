

#include <mpi.h>
int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	
	MPI_Status status;
	int myRank,npes,send[10],recieve[10];
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	
	for(int i =  0; i < 10 ; i++){
		send[i] = myRank;
	}
	
	
	if(myRank%2 ==  0){
		printf("Hey I am proc %d Im sending\n",myRank);
		MPI_Send(send,10,MPI_INT,(myRank+1)%npes,0,MPI_COMM_WORLD); 
		printf("Hey I am proc %d Im recieving\n",myRank);
		MPI_Recv(recieve,10,MPI_INT,(myRank-1+npes)%npes,0,MPI_COMM_WORLD,&status); 
	}else{
		printf("Hey I am proc %d Im recieving\n",myRank);
		MPI_Recv(recieve,10,MPI_INT,(myRank-1+npes)%npes,0,MPI_COMM_WORLD,&status); 
		printf("Hey I am proc %d Im sending\n",myRank);
		MPI_Send(send,10,MPI_INT,(myRank+1)%npes,0,MPI_COMM_WORLD); 
	}
	
	
	for(int i =  0; i < 10 ; i++){
		printf ("Hey I am proc %d, I recieve from %d -> b[%d] = %d \n",myRank,(myRank-1+npes)%npes,i,recieve[i]);
	}
		
	MPI_Finalize(); 
}
