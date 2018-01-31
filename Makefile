all: 
	mpicxx -o test ./*.cpp 

run: all
	mpirun -n 2 ./test
