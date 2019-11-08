# compile the program

compile: mpi-APSP.c fox
	rm -f -- fox
	mpicc mpi-APSP.c -o fox -lm
