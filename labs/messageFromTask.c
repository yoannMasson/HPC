#include <stdio.h>

#include <mpi.h>
int main(int argc, char *argv[]){

	MPI_Init(&argc, &argv);
	
	int npes,myRank,a[10];
	double t0,t1;
	MPI_Request request;
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &npes);
	
	if(myRank == 0){
		for(int i =  0; i < 10 ; i++){//Populating the array
			a[i] = i;
		}
	}
	MPI_Status status;
	
	t0 = MPI_Wtime();
	if(myRank == 0){
		for(int i = 1;i < npes ; i++){
			MPI_Isend(a,10,MPI_INT,i,1,MPI_COMM_WORLD,&request); 
		}		
	}else{
		MPI_Recv(a,10,MPI_INT,0,1,MPI_COMM_WORLD,&status); 
		printf("Hey, this is process %d, here is the first element of the array: a[0] = %d\n",myRank,a[0]);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	t1 = MPI_Wtime();
	if(myRank == 0){
		printf("Hey, this is process %d, it took %f secondes to send everything\n",myRank,t1-t0);
	}
	int b;
	scanf("%d",&b);
	printf("\n %d \n",b);
	scanf("%d",&b);
	printf("\n %d \n",b);
	scanf("%d",&b);
	printf("\n %d \n",b);
	MPI_Finalize();
	 
}
