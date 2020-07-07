### Description of the Codes
The program started with the addition of libraries as follows. Besides the standard libraries, the required "mpi.h" library for MPI has been added. 1000000 value has been determined for infinity. If there is no edge between any two vertices the weight is the constant

```sh
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#define INFINITY 1000000
 ```
The following section has definitions of functions. These functions will be defined where they are created.
```sh
int matrix_size(int rank, MPI_Comm comm);
MPI_Datatype Build_blk_col_type(int n, int loc_n);
void Read_matrix(int loc_mat[], int n, int loc_n, MPI_Datatype blk_col_mpi_t, int rank, MPI_Comm comm);
void Dijkstra_Init(int loc_mat[], int loc_pred[], int loc_dist[], int loc_known[],int rank, int loc_n);
void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n, MPI_Comm comm);
int find_min_distance(int loc_dist[], int loc_known[], int loc_n);
void print_matrix(int mat[]);
void print_distances(int distance[], int n);
int * matrix();
```

This section contains explanations about the main method of the program. In order to use MPI, it is necessary to make definitions at the beginning of the program.
```sh
int main(int argc, char** argv) {
int* loc_mat, * loc_dist, * loc_pred, * distance = NULL, * global_pred = NULL;
int rank, p, loc_n, n;
MPI_Comm comm; // This is a MPI method. This is used to determine the number of processes in a computation.
MPI_Datatype blk_col_mpi_t;
MPI_Init(NULL, NULL); 
```
We use this code to initiate an MPI computation. Parameters are required only in the C language binding where they are the main program arguments. We put NULL for the parameter requirement.
```sh
comm = MPI_COMM_WORLD; // A set of attributes that describe the execution environment are attached to the communicator MPI_COMM_WORLD .
MPI_Comm_rank(comm, &rank);// In MPI, this code used to determine the process identifier of the current process.Comn means communicator , rank is number of processes in the group.
MPI_Comm_size(comm, &p); // This line number of processes in computation
```

In the other part of the code, there is the continuation of the main method. Variables begin with "loc_" to show that they are locally defined within the main.
```sh
n = matrix_size(rank, comm);// This part is the matrix size we use. A constant value for our project. But it is defined by the method for the case of wanting to be determined from the outside. Explanation will be made in the method section.
loc_n = n / p;
loc_mat =(int*) malloc(n * loc_n * sizeof(int)); //It means matrix.
loc_dist = (int*)malloc(loc_n * sizeof(int)); //It means distance
loc_pred = (int*)malloc(loc_n * sizeof(int)); //It means predecessor.
blk_col_mpi_t = Build_blk_col_type(n, loc_n); // This is a method. The method has a description where it is located.
if (rank == 0) // If the process rank is 0 it is necessary to limit the size by int.
    {
     distance = (int*)malloc(n * sizeof(int));
     global_pred = (int*)malloc(n * sizeof(int));
    }
Read_matrix(loc_mat, n, loc_n, blk_col_mpi_t, rank, comm);
Dijkstra(loc_mat, loc_dist, loc_pred, loc_n, n, comm);}

```
They are methods. The method has a description where it is located.To accumulate an output array, Gather the results from Dijkstra. 32 and 33 are used for this. 
```sh
	MPI_Gather(loc_dist, loc_n, MPI_INT, distance, loc_n, MPI_INT, 0, comm);
	MPI_Gather(loc_pred, loc_n, MPI_INT, global_pred, loc_n, MPI_INT, 0, comm);
```
In the section below, the results are printed on the screen. Then, the variables determined at the beginning are made free.
```sh
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
return 0;}
```

The explanations of the Main method are like this. the next section has a description of the methods used. The first method used is matrix_size. Here is the matrix size that we defined at the beginning of the Main method. The size of our matrix is 8x8.
```sh
int matrix_size( MPI_Comm comm) {
     int n;
     n = 8;
     MPI_Bcast(&n, 1, MPI_INT, 0, comm);
     return n; }
```

MPI_BCAST is used to broadcast the problem size parameter (size) from process 0 to all number of processes. Size is included in this section. Because of this the definition is made here. Next is the Build_blk_col_type method, defined by MPI_Datatype.Purpose of this method is building an MPI_Datatype that represents a block column of a matrix. In method, n is number of rows in the matrix and the block column. Also,   loc_n = n/p is number cols in the block column.

```sh
MPI_Datatype Build_blk_col_type(int n, int loc_n) {
MPI_Aint lb, extent;
.MPI_Datatype block_mpi_t;
MPI_Datatype first_bc_mpi_t;
MPI_Datatype blk_col_mpi_t; //This code is MPI_Datatype that represents a block column.
MPI_Type_contiguous(loc_n, MPI_INT, &block_mpi_t);
MPI_Type_get_extent(block_mpi_t, &lb, &extent);
MPI_Type_vector(n, loc_n, n, MPI_INT, &first_bc_mpi_t);// Parameters needed by the command like (numblocks, elts_per_block, stride, oldtype, *newtype) 
MPI_Type_create_resized(first_bc_mpi_t, lb, extent, &blk_col_mpi_t); // This call is needed to get the right extent of the new data type.
MPI_Type_commit(&blk_col_mpi_t);
MPI_Type_free(&block_mpi_t);
MPI_Type_free(&first_bc_mpi_t);
return blk_col_mpi_t; }
```
Next is the Read_matrix method. In this section, the matrix with the values to be processed is printed on the screen. MPI_SCATTER distribute an input arra from process 0 to other processes.

```sh
void Read_matrix(int loc_mat[], int n, int loc_n,
MPI_Datatype blk_col_mpi_t, int rank, MPI_Comm comm) {
int* mat= NULL;
if (rank == 0) {
            mat = (int*)matrix();// This is a method. The method has a description where it is located.
        }
printf("\n");
printf(" This project has been prepared by Emre Can Tural (151180061):\n");
printf("\n");
print_matrix(mat);
MPI_Scatter(mat, 1, blk_col_mpi_t, loc_mat, n * loc_n, MPI_INT, 0, comm);}
```
Next, there are explanations of the matrix method we use. We turn the two-dimensional matrix into one dimension so that we can make the calculation. We do this in the order in the table below. Each element in our main matrix is equal to the sequence number printed on the new matrix.
- Table 1.1: Transformation table of the matrix.

| N| 1 | 2| 3|  4| 5| 6| 7| 8|
| ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | ------ | 
| 1| 	0| 	40| 	15| 	Inf| 	Inf| 	Inf| 	35| 	Inf| 
| 2| 	40	| 0| 	0| 	100| 	Inf| 	Inf| 	25| 	Inf| 
| 3| 	15	| 0	| 0| 	20| 	10| 	50| 	Inf| 	50| 
| 4| 	Inf	| 100| 	20	| 0| 	10	| Inf| 	45| 	Inf| 
| 5| 	Inf	| Inf| 	10| 	10	| 0	| 30	| 50	| 30| 
| 6| 	Inf	| Inf| 	50| 	Inf	| 30| 	0	| Inf| 	0| 
| 7| 	35	| 25| 	Inf| 	45	| 50| 	Inf	| 0	| 0| 
| 8| 	Inf	| Inf| 	50| 	Inf	| 30| 	0	| 0| 	0| 


Next is the Dijkstra_Init method. Purpose of this method is initializing all the matrices for  Dijkstra's shortest path. In the method, loc_mat is local matrix containing edge costs between verticesi, loc_dist is loc_dist[v] = shortest distance from the source to each vertex v,loc_pred is  loc_pred[v] = predecessor of v on a shortest path from source to v, loc_known isloc_known[v] = 1 if vertex has been visited, 0 else.
```sh
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
```
Now it is our Dijkstra method, which is our main purpose. This method compute all the shortest paths from the source vertex 0 to all vertices. In the method, loc_mat is local matrix containing edge costs between vertices, loc_n is    local number of vertices and comm is the communicator, n is total global number of vertices, loc_dist is  loc_dist[v] = shortest distance from the source to each vertex.
```sh
void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n,
MPI_Comm comm) {
int i, loc_v, loc_u, glbl_u, new_dist, rank, dist_glbl_u;
int* loc_known;
int my_min[2];
int glbl_min[2];
MPI_Comm_rank(comm, &rank);
loc_known = (int*)malloc(loc_n * sizeof(int));
Dijkstra_Init(loc_mat, loc_pred, loc_dist, loc_known, rank, loc_n);
 for (i = 0; i < n - 1; i++) { //Loop run n - 1 times since we already know the shortest path to global vertex 0 from global vertex 0.
        loc_u = find_min_distance(loc_dist, loc_known, loc_n);
        if (loc_u != -1) {
            my_min[0] = loc_dist[loc_u];
            my_min[1] = loc_u + rank * loc_n;
        }
        else {
            my_min[0] = INFINITY;
            my_min[1] = -1;
        }  
```
MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC, comm);// MPI_Allreduce determine the maximum of a set of localer values computed at the different processes and to distribute this maximum value to each process.  It get the minimum distance found by the processes and store that distance and the global vertex in glbl_min.
    
```sh
        dist_glbl_u = glbl_min[0];
        glbl_u = glbl_min[1];
        if (glbl_u == -1)// This code is to assure that loc_known is not accessed with -1. In the other section, code check if global u belongs to process, and if so  it update loc_known.
            break;
        if ((glbl_u / loc_n) == rank) {
            loc_u = glbl_u % loc_n;
            loc_known[loc_u] = 1;
        }
```
For each local vertex (global vertex = loc_v + my_rank * loc_n), this code update the distances from source vertex (0) to loc_v.       
```sh
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
```
After all these steps, we will find the minimum local distance from the source to the assigned vertices of the process that calls the method. We use the find_min_distance method for this task. In the method, In args is   loc_dist:  array with distances from source 0   , loc_known is array with values 1 if the vertex has been visited, loc_n is  local number of vertices, return value  loc_u is the vertex with the smallest value in loc_dist (-1 if all vertices are already known). Also, loc_u = -1 is not supposed to be used when this function returns.

```sh
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
```
The next method is print matrix. The output in Figure 1.2 belongs to this method. The previously defined series in the method was arranged according to the desired form. If the value is INFINITY, it is named as "Inf" as desired.
```sh
oid print_matrix(int mat[]) {
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
```
The last method is the print_distances method. If the distance was infinite, "Inf" would be written. But there is no infinite value for our matrix.

```sh
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
```
