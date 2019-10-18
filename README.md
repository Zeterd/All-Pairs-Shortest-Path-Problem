# All-Pairs-Shortest-Path-Problem

## Description

### About

The problem to be solved consists in determining all shortest paths between pairs of nodes in a given graph.

For more information about the implementation problem, see the file "Description" in this repository

### Language

The programming language used is C with the [MPI library][1].

[1]: http://condor.cc.ku.edu/~grobe/docs/intro-MPI-C.shtml

## Made by

A program made by [José Pedro][2] for a college assignment of computer science from Faculty of Science of University of Porto.

<!-- [2]: Biographie link about the author --->

## How to run

To run the program, you must have the MPI software installed in your device. If you don't have You can go to [MPI Home Page][3] and see the process to installation.

After install the software, you can compile first the program by:

```shell
mpicc program.c
```

And to run the program after compile:

```shell
mpirun -np 2 a.out
```

The -np flag means how many processors(cores) you want to use to run the program(you only can run the number of max the cores of your CPU have or it will fail the execution).

[3]: https://mpitutorial.com/tutorials/
