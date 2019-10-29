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

//The operation to multiply the values
void operation_multiply(MATRIX* m1, MATRIX* m2, MATRIX* m3){
    int i,j,k;

    for(i=0; i<m1->dim; i++){
        for(j=0; j<m1->dim; j++){
            for(k=0; k<m2->dim; k++){
              if((m1->entries[i][k]+m2->entries[k][j])<m3->entries[i][j])
                m3->entries[i][j] = m1->entries[i][k]+m2->entries[k][j];
            }
        }
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


    if((n%(int)sqrt(p)) != 0 || floor(sqrt(p)) != (int)sqrt(p)){
        printf("Algorithm not apply, Aborting!!!\n");
        return 0;
    }
    else
      return 1;
}

//Creation of the grid
void create_grid(GRID_TYPE *grid){
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
    MATRIX* m_tmp;


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
//print the matrix
void print_matrix(int **m, int dim){
    int i,j;
    for(i=0;i<dim;i++){
        printf("%d", m[i][0]);
        for(j=0;j<dim; j++){
            printf(" %d", m[i][j]);
        }
        printf("\n");
    }
}
//Basacly the name explains it self
int** copy_matrix(int **m, int dim){
    int **m_tmp, i, j;
    malloc_matrix(&m_tmp, dim);

    for(i=0; i<dim; i++){
        for(j=0; j<dim; j++){
            m_tmp[i][j] = m[i][j];
        }
    }
    return m_tmp;
}
//ask
void submatrix(int **m, GRID_TYPE* grid, MATRIX* matrix_aux){
    int i,j;
    for(i=0; i<matrix_aux->dim; i++){
        for(j=0; j<matrix_aux->dim; j++){
            int index_i = grid->my_row*matrix_aux->dim+i;
            int index_j = grid->my_col*matrix_aux->dim+j;
            matrix_aux->entries[i][j] = m[index_i][index_j];
        }
    }
}
//fill the matrix with a value, used to fill the result matrix with infinit
void fill_matrix(MATRIX* m, int v){
    int i,j;

    for(i=0; i<m->dim; i++){
      for(j=0; j<m->dim; j++){
          m->entries[i][j] = v;
      }
    }
}

//This is the main :/
int main(int argc, char *argv[]) {
    GRID_TYPE grid;
    MATRIX *matrix_1, *matrix_2, *matrix_3;
    int half_dim;
    double finish, start;
    int n_procs, rank, mat_d, **matrix, **matrix_res, f=0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank == ROOT){
      scanf("%d", &mat_d);

      if(flag_func(n_procs, mat_d) == 1)
        f=1;
    }

    if(f == 1){
      MPI_Finalize();
      return 0;
    }

    //set the flag to the others processors
    MPI_Bcast(&f, 1, MPI_INT, ROOT, MPI_COMM_WORLD);


    MPI_Bcast(&mat_d, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
    create_grid(&grid);

    if(rank == ROOT){
        matrix = createMatrix(mat_d);
    }


    MPI_Bcast(&(matrix[0][0]), mat_d, MPI_INT, ROOT, MPI_COMM_WORLD);

    half_dim = mat_d/grid.q;
    matrix_1 = malloc_MATRIX(half_dim);
    matrix_2 = malloc_MATRIX(half_dim);
    matrix_3 = malloc_MATRIX(half_dim);

    submatrix(matrix, &grid, matrix_1);

    matrix_2->entries = copy_matrix(matrix_1->entries, half_dim);
    fill_matrix(matrix_3, INFINITO);

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    //the "start" of execution of the algorithm
    int i;
    for(i=1; i<mat_d-1; i*=2){
        fox(&grid, matrix_1, matrix_2, matrix_3, half_dim);
        matrix_1->entries = copy_matrix(matrix_3->entries, half_dim);
        matrix_2->entries = copy_matrix(matrix_3->entries, half_dim);
    }

    malloc_matrix(&matrix_res, half_dim);
    int j;

    if(rank == ROOT){
        MPI_Send(&(matrix_3->entries[0][0]), half_dim*half_dim, MPI_INT, ROOT, 0, MPI_COMM_WORLD);

    }

    else {
        for(i=0; i<half_dim; i++){
              for(j=0; j<half_dim; j++){
                  if(matrix_3->entries[i][j] == INFINITO)
                    matrix[i][j] = 0;
                  else
                    matrix[i][j] = matrix_3->entries[i][j];
              }
        }
        for(i=1;i<n_procs;i++){
            MPI_Recv(&(matrix_res[0][0]), half_dim*half_dim, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int r = i/grid.q;
            int c = i%grid.q;


            int k;
            for(k=0; k<half_dim; k++){
              for(j=0; j<half_dim; j++){
                if(matrix_res[k][j] == INFINITO)
                    matrix[r*half_dim+k][c*half_dim+j] = 0;
                else
                  matrix[r*half_dim+k][c*half_dim+j] = matrix_res[k][j];
              }
            }

        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();

    if(rank == ROOT){
        printf("Execution time: %lf\n", finish-start);
        print_matrix(matrix, mat_d);
    }

    MPI_Finalize();
    return 0;
}
