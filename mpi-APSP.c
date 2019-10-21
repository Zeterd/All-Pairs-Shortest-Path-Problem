#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

#define ROOT 0

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
            matrix[i][j] = MAX_VALUE; //Meter no futuro +INFINITO
          }
          else{
            matrix[i][j] = num;
          }
        }
    }

    return matrix;
}

int main(int argc, char *argv[]) {
    int numprocs, rank, q, mat_d, **matrix;

    MPI_Init(&argc, &argv);                                                                                                               │
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == ROOT){
      scanf("%d", &mat_d);
      createMatrix(mat_d);
    }

    q = sqrt(numprocs);

    if(mat_d%q == 0){
      printf("Error in modulo\n");
      return;
    }

    for(int i=1; i<=mat_d; i++){
        
    }




}
