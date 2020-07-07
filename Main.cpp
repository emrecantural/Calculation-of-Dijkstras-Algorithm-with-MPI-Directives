#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define INFINITY 1000000

int matrix_size( MPI_Comm comm);
MPI_Datatype Build_blk_col_type(int n, int loc_n);
void Read_matrix(int loc_mat[], int n, int loc_n, MPI_Datatype blk_col_mpi_t, int rank, MPI_Comm comm);
void Dijkstra_Init(int loc_mat[], int loc_pred[], int loc_dist[], int loc_known[],int rank, int loc_n);
void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n, MPI_Comm comm);
int find_min_distance(int loc_dist[], int loc_known[], int loc_n);
void print_matrix(int mat[]);
void print_distances(int distance[], int n);
int * matrix();

int main(int argc, char** argv) {
    int* loc_mat, * loc_dist, * loc_pred, * distance = NULL, * global_pred = NULL;
    int rank, p, loc_n, n;
     MPI_Comm comm;
    MPI_Datatype blk_col_mpi_t;
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);
    n = matrix_size(comm);
    loc_n = n / p;
    loc_mat =(int*) malloc(n * loc_n * sizeof(int));
    loc_dist = (int*)malloc(loc_n * sizeof(int));
    loc_pred = (int*)malloc(loc_n * sizeof(int));
    blk_col_mpi_t = Build_blk_col_type(n, loc_n);

    if (rank == 0) {
        distance = (int*)malloc(n * sizeof(int));
        global_pred = (int*)malloc(n * sizeof(int));
    }
    Read_matrix(loc_mat, n, loc_n, blk_col_mpi_t, rank, comm);
    Dijkstra(loc_mat, loc_dist, loc_pred, loc_n, n, comm);
    MPI_Gather(loc_dist, loc_n, MPI_INT, distance, loc_n, MPI_INT, 0, comm);
    MPI_Gather(loc_pred, loc_n, MPI_INT, global_pred, loc_n, MPI_INT, 0, comm);
    if (rank == 0) {

        print_distances(distance, n);
        free(distance);
        free(global_pred);
    }
    free(loc_mat);
    free(loc_pred);
    free(loc_dist);
    MPI_Type_free(&blk_col_mpi_t);
    MPI_Finalize();
    getchar();
    return 0;
}
int matrix_size( MPI_Comm comm) {
    int n;
    n = 8;
    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    return n;
}
MPI_Datatype Build_blk_col_type(int n, int loc_n) {
    MPI_Aint lb, extent;
    MPI_Datatype block_mpi_t;
    MPI_Datatype first_bc_mpi_t;
    MPI_Datatype blk_col_mpi_t;
    MPI_Type_contiguous(loc_n, MPI_INT, &block_mpi_t);
    MPI_Type_get_extent(block_mpi_t, &lb, &extent);
    MPI_Type_vector(n, loc_n, n, MPI_INT, &first_bc_mpi_t);
    MPI_Type_create_resized(first_bc_mpi_t, lb, extent, &blk_col_mpi_t);
    MPI_Type_commit(&blk_col_mpi_t);
    MPI_Type_free(&block_mpi_t);
    MPI_Type_free(&first_bc_mpi_t);
    return blk_col_mpi_t;
}
void Read_matrix(int loc_mat[], int n, int loc_n,
    MPI_Datatype blk_col_mpi_t, int rank, MPI_Comm comm) {
    int* mat= NULL;
    if (rank == 0) {
        mat = (int*)matrix();
    }
    printf("\n");
    printf(" This project has been prepared by Emre Can Tural:\n");
    printf("\n");
    print_matrix(mat);
    MPI_Scatter(mat, 1, blk_col_mpi_t, loc_mat, n * loc_n, MPI_INT, 0, comm);
}
int * matrix() {
    static int  mat[64];
    mat[0] = mat[9] = mat[10] = mat[17] = mat[18] = mat[27] = mat[36] = 0;
    mat[45] = mat[47] = mat[54] = mat[55] = mat[61] = mat[62] = mat[63] = 0;
    mat[3] = mat[4] = mat[5] = mat[7] = mat[12] = mat[13] = mat[15] = mat[22] = INFINITY;
    mat[24] = mat[29] = mat[31] = mat[32] = mat[33] = mat[40] = mat[41] = INFINITY;
    mat[43] = mat[46] = mat[50] = mat[53] = mat[56] = mat[57] = mat[59] = INFINITY;
    mat[20] = mat[28] = mat[34] = mat[35] = 10;
    mat[2] = mat[16] = 15;
    mat[19] = mat[26] = 20;
    mat[14] = mat[49] = 25;
    mat[37] = mat[39] = mat[44] = mat[60] = 30;
    mat[6] = mat[48] = 35;
    mat[1] = mat[8] = 40;
    mat[30] = mat[51] = 45;
    mat[21] = mat[23] = mat[38] = mat[42] = mat[52] = mat[58] = 50;
    mat[11] = mat[25] = 100;
    return mat;
}
void Dijkstra_Init(int loc_mat[], int loc_pred[], int loc_dist[], int loc_known[],
    int rank, int loc_n) {
    int loc_v;

    if (rank == 0)
        loc_known[0] = 1;
    else
        loc_known[0] = 0;

    for (loc_v = 1; loc_v < loc_n; loc_v++)
        loc_known[loc_v] = 0;

    for (loc_v = 0; loc_v < loc_n; loc_v++) {
        loc_dist[loc_v] = loc_mat[0 * loc_n + loc_v];
        loc_pred[loc_v] = 0;
    }
}

void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n,
    MPI_Comm comm) {

    int i, loc_v, loc_u, glbl_u, new_dist, rank, dist_glbl_u;
    int* loc_known;
    int my_min[2];
    int glbl_min[2];

    MPI_Comm_rank(comm, &rank);
    loc_known = (int*)malloc(loc_n * sizeof(int));

    Dijkstra_Init(loc_mat, loc_pred, loc_dist, loc_known, rank, loc_n);
    for (i = 0; i < n - 1; i++) {
        loc_u = find_min_distance(loc_dist, loc_known, loc_n);

        if (loc_u != -1) {
            my_min[0] = loc_dist[loc_u];
            my_min[1] = loc_u + rank * loc_n;
        }
        else {
            my_min[0] = INFINITY;
            my_min[1] = -1;
        }
        MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC, comm);
        dist_glbl_u = glbl_min[0];
        glbl_u = glbl_min[1];
        if (glbl_u == -1)
            break;
        if ((glbl_u / loc_n) == rank) {
            loc_u = glbl_u % loc_n;
            loc_known[loc_u] = 1;
        }
        for (loc_v = 0; loc_v < loc_n; loc_v++) {
            if (!loc_known[loc_v]) {
                new_dist = dist_glbl_u + loc_mat[glbl_u * loc_n + loc_v];
                if (new_dist < loc_dist[loc_v]) {
                    loc_dist[loc_v] = new_dist;
                    loc_pred[loc_v] = glbl_u;
                }
            }
        }
    }
    free(loc_known);
}

int find_min_distance(int loc_dist[], int loc_known[], int loc_n) {
    int loc_u = -1, loc_v;
    int shortest_dist = INFINITY;

    for (loc_v = 0; loc_v < loc_n; loc_v++) {
        if (!loc_known[loc_v]) {
            if (loc_dist[loc_v] < shortest_dist) {
                shortest_dist = loc_dist[loc_v];
                loc_u = loc_v;
            }
        }
    }
    return loc_u;
}


void print_matrix(int mat[]) {
    int i, j;
    printf("\n"); 
    printf("  Distance matrix:\n");
    printf("\n");   
    printf( "    N      1        2       3      4       5       6        7      8");
    printf("\n");  
    printf("\n");
    for (i = 0; i < 8; i++) {
        printf("%5d  ",i+1);
        for (j = 0; j < 8; j++)
        {
            if (mat[i * 8 + j] == INFINITY)
                printf("%5s\t", "Inf");
            else
                printf("%5d\t", mat[i * 8 + j]);
        }
        printf("\n");
        printf("\n");
    }

    printf("\n");
}


void print_distances(int distance[], int n) {
    int v;

    printf("\n");
    printf("Minimum distances from node 1:\n");
    printf("\n");

    for (v = 0; v < n; v++) {
        if (distance[v] == INFINITY) {
            printf("%3d       %5s\n", v+1, "inf");
        }
        else
            printf("%3d       %4d\n", v+1, distance[v]);
    }
    printf("\n");
}
