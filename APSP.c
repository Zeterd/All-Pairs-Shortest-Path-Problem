/****
   All Pairs Shortest Paths
   Given a directed, connected weighted graph G(V,E), for each edge ⟨u,v⟩∈E, a weight w(u,v) is associated with the edge.
The all pairs of shortest paths problem (APSP) is to find a shortest path from u to v for every pair of vertices u and v in V.

    Made by Jose Pedro
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_VALUE 9999999


int **matrix_sum(int **matrix, int **matrix2, int dim);
int minValue(int x, int y);
int maxValue(int x, int y);
void printMatriz(int mat_d, int **new_mat);

int main() {
    int i=0, j=0;
    int mat_d;

    scanf("%d", &mat_d);

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

    int m = 1;
    int n_rows = mat_d;
    int **new_mat;

    new_mat = matrix;
    //printMatriz(mat_d, new_mat);

    while(m<n_rows-1){
        new_mat = matrix_sum(new_mat, new_mat, mat_d);

        m = 2*m;
    }

    printMatriz(mat_d, new_mat);
    return 0;
}

int **matrix_sum(int **matrix, int **matrix2, int dim){
    int i=1, j=1, k=1;
    int **matrix3;

    matrix3 = malloc(sizeof(int*) * dim+1);


    for(i = 1; i <= dim; i++) {
        matrix3[i] = (int*)malloc(sizeof(int*) * dim+1);
    }

    //printf("ESTOU AQUUII\n");

    for(i = 1; i <= dim; i++){
        for(j = 1; j <= dim; j++){
            matrix3[i][j] = MAX_VALUE;
            //printf("%d\n", matrix3[i][j]);
            for(k = 1; k <= dim; k++){
                //printf("cij:%d  aik:%d   bkj:%d\n", matrix3[i][j], matrix[i][k], matrix2[k][j]);
                matrix3[i][j] = minValue(matrix3[i][j], (matrix[i][k]+matrix2[k][j]));
                //printf("minvalue:%d\n\n\n", matrix3[i][j]);
            }
        }
    }
    //printMatriz(dim, matrix3);
    return matrix3;

}


int minValue(int x, int y){
    return (x<y) ? x : y;
}
int maxValue(int x, int y){
    return (x>y) ? x : y;
}

void printMatriz(int mat_d, int **new_mat){
  int i, j;
  for(i=1;i<=mat_d;i++){
    for(j=1;j<=mat_d;j++){
      if(new_mat[i][j] == MAX_VALUE){
        new_mat[i][j] = 0;
        printf("%d", new_mat[i][j]);
      }
      else
        printf("%d", new_mat[i][j]);
      printf(" ");
    }
    printf("\n");
  }
  printf("\n");

}
