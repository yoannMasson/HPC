all: 
	mpicxx ./*.cpp 

run: all
	mpirun -n 4 ./a.out
