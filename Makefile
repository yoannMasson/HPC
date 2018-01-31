all: 
	mpicxx ./*.cpp 

run: all
	mpirun -n 2 ./a.out
