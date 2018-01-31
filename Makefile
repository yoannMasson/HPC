all: 
	mpicxx -o test ./*.cpp 

run: all

	mpirun -n 4 ./test
