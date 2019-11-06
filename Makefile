# compile the program

compile: mpi-APSP.c
	mpicc mpi-APSP.c -o fox -lm
