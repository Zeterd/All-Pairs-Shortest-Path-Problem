#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "aux.h"

#define ROOT 0
#define INFINITO 99999

int** createMatrix(int mat_d){
    int **matrix;

    matrix = malloc(sizeof(int*) * mat_d+1);

    for(i = 1; i <= mat_d; i++) {
      matrix[i] = (int*)malloc(sizeof(int*) * mat_d+1);
    }

    for(i=1; i<=mat_d; i++){
        for(j=1; j<=mat_d; j++){
          int num;
          scanf("%d", &num);

          if(num == 0 && i != j){
            matrix[i][j] = INFINITO; //Meter no futuro +INFINITO
          }
          else{
            matrix[i][j] = num;
          }
        }
    }

    return matrix;
}

//Verify if is possible to construct matix
int flag_func(int p, int n){
    if(n%sqrt(p) != 0 || floor(sqrt(p)) != sqrt(p)){
        printf("Algorithm not apply, Aborting!!!\n");
        return 0;
    }
    else
      return 1;
}

int create_grid(GRID_TYPE *grid){
    MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));

    //???

}

MATRIX* malloc_matrix(int dim){
    MATRIX* tmp = (MATRIX*malloc(sizeof(MATRIX)));
    tmp->dim = dim;
    //Falta alocar a matrix em memoria

    return tmp;
}

int main(int argc, char *argv[]) {
    int n_procs, rank, q, mat_d, **matrix, f=0;

    MPI_Init(&argc, &argv);                                                                                                               │
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == ROOT){
      scanf("%d", &mat_d);

      if(flag_func(n_procs, mat_d) == 1)
        f=1;
    }


    //set the flag to the others processors
    MPI_Bcast(&f, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

    if(f != 1){
        MPI_Bcast(&mat_d, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
        create_grid(&grid);

        if(my_rank == ROOT){
            matrix = createMatrix(mat_d);
        }
    }

    MPI_Bcast(&(matrix[0][0]), mat_d, MPI_INT, ROOT, MPI_COMM_WORLD);

    half_dim = mat_d/grid->q;
    matrix1 = malloc_matrix(half_dim);





    //Operaçao para todos os processos fazer(multiplicaçao)

    //======

    //Operaçao de agregaçao dos dados de cada processo para a ROOT --> MPI_GATHER








}
