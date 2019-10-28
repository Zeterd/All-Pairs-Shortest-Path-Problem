//Struct of a grid
typedef struct {
    int       p;        //processors of the grid
    MPI_Comm  comm;     
    MPI_Comm  row_comm;
    MPI_Comm  col_comm;
    int       q;
    int       my_row;
    int       my_col;
    int       my_rank;
} GRID_INFO_TYPE;

//struct of a Matriz
typedef struct {
    int  n_bar;
    int  **entries;
} MATRIX;
