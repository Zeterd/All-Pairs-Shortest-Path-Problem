#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include "aux.h"

#define ROOT 0
#define INFINITO 99999

void malloc_matrix(int ***matrix, int dim){
    //Allocate dim^2 contiguous items
    int *p = (int*)malloc(dim*dim*sizeof(int*));

    //Allocate the row pointers into the memory
    *matrix = (int**)malloc(dim*sizeof(int*));

    //Set uo pointers
    int i;
    for(i=0; i<dim; i++){
        *matrix[i] = &(p[i*dim]);
    }
}


int** createMatrix(int mat_d){
    int **matrix;
    int i,j;

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
    if((n%(int)sqrt(p)) != 0 || floor(sqrt(p)) != sqrt(p)){
        printf("Algorithm not apply, Aborting!!!\n");
        return 0;
    }
    else
      return 1;
}

//Creation of the grid
int create_grid(GRID_TYPE *grid){
    int dim[2], vect_flag[2], coor[2], coor2[2];

    MPI_Comm_size(MPI_COMM_WORLD, &(grid->procs));

    grid->q = (int)sqrt(grid->procs);
    dim[0] = grid->q;
    dim[1] = grid->q;

    vect_flag[0] = 1;
    vect_flag[1] = 1;

    //Makes a new communicator to which topology information has been attached
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, vect_flag, 1, &(grid->comm));
    MPI_Comm_rank(grid->comm, &(grid->my_rank));

    //Determines process coords in cartesian topology given rank in group
    MPI_Cart_coords(grid->comm, grid->my_rank, 2, coor);

    grid->my_row = coor[0];
    grid->my_col = coor[1];

    //ROW_COMM_WORLD
    coor2[0] = 0;
    coor2[1] = 1;

    //Partitions a communicator into subgroups which form lower-dimensional cartesian subgrids
    MPI_Cart_sub(grid->comm, coor2, &(grid->row_comm));

    //COLUMN_COMM_WORLD
    coor2[0] = 1;
    coor2[1] = 0;
    MPI_Cart_sub(grid->comm, coor2, &(grid->col_comm));
}

//Allocation of the matrix type MATRIX structer
MATRIX* malloc_MATRIX(int dim){
    MATRIX* tmp = (MATRIX*)malloc(sizeof(MATRIX));
    tmp->dim = dim;
    malloc_matrix(&(tmp->entries), dim);
    return tmp;
}

//Implementation of the FOX algorithm
void fox(GRID_TYPE* grid, MATRIX* m1, MATRIX* m2, MATRIX* m3, int dim){
    int org, dest, i, root;


    org = (grid->my_row+1) % grid->q;
    dest = (grid->my_row + grid->q-1) % grid->q;

    m_tmp = malloc_MATRIX(dim);
    for(i=0; i<grid->q-1; i++){
        root = (grid->my_row+i)%grid->q;
        if(root == grid->my_col){
            MPI_Bcast(&(m1->entries[0][0]), dim*dim, MPI_INT, root, grid->row_comm);
            operation_multiply(m1,m2,m3);
        }
        else{
            MPI_Bcast(&(m_tmp->entries[0][0]), dim*dim, MPI_INT, root, grid->row_comm);
            operation_multiply(m_tmp,m2,m3);
        }

        MPI_Send(&(m2->entries[0][0]), dim*dim, MPI_INT, dest, 0, grid->col_comm);
        MPI_Recv(&(m2->entries[0][0]), dim*dim, MPI_INT, org, 0, grid->col_comm, MPI_STATUS_IGNORE);
  
    }
}


//This is the main :/
int main(int argc, char *argv[]) {
    GRID_TYPE grid;
    MATRIX *matrix_1, *matrix_2;
    int half_dim;
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

        if(rank == ROOT){
            matrix = createMatrix(mat_d);
        }
    }

    MPI_Bcast(&(matrix[0][0]), mat_d, MPI_INT, ROOT, MPI_COMM_WORLD);

    half_dim = mat_d/grid.q;
    matrix_1 = malloc_MATRIX(half_dim);
    matrix_2 = malloc_MATRIX(half_dim);
    matrix_3 = malloc_MATRIX(half_dim);

    get_submatrix(matrix, &grid, matrix_1)

    matrix_2->entries = copy_matrix(matrix_1, half_dim);
    fill_matrix(matrix_3, INFINITO);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    //the "start" of execution of the algorithm
    int i;
    for(i=1; i<mat_d-1; i*=2){
        fox(&grid, matrix_1, matrix_2, matrix_3, half_dim);



    }

}
