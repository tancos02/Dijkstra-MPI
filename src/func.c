#include <stdio.h>
#include <stdlib.h>
#include </usr/include/mpi/mpi.h>
#include <time.h>
#include "func.h"
#define INFINITY 1000000

/*
    Fungsi untuk membaca n (jumlah node yang akan digunakan)
 */
int Read_n(int my_rank, MPI_Comm comm) {
    int n;
    if (my_rank == 0){
        scanf("%d", &n);
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    return n;
}

/*
    Membuat MPI_Datatype yang merepresentasikan matriks
 */
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

/*
    Membuat sebuah matriks nxn yang mewakili graf.
    Nilai M[i][j] = jarak dari node i ke j.
    Jarak diinisialisasi dengan nilai random yang memiliki seed = 13517111.
 */
void Make_matrix(int loc_mat[], int n, int loc_n,
                 MPI_Datatype blk_col_mpi_t, int my_rank, MPI_Comm comm) {
    int *mat = NULL, i, j;
    int upper = 1000;
    int lower = 0;
    srand(13517111);

    if (my_rank == 0) {
        mat = malloc(n * n * sizeof(int));
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if(j>i) {
                    int x = (rand() % (upper - lower + 1)) + lower; 
                    mat[i * n + j] = x;
                    mat[i + j * n] = x;
                }
                else if (i==j){
                    mat[i * n + j] = 0;
                }
            }
        }
    }

    MPI_Scatter(mat, 1, blk_col_mpi_t, loc_mat, n * loc_n, MPI_INT, 0, comm);

    if (my_rank == 0) free(mat);
}

/*
    Menginisialisasi setiap matriks yang ada pada nodes agar algoritmanya bisa berjalan.
 */
void Dijkstra_Init(int loc_mat[], int loc_pred[], int loc_dist[], int loc_known[],
                   int my_rank, int loc_n) {
    int loc_v;

    if (my_rank == 0)
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

/*
    Algoritma dijkstra untuk menghitung jarak terpendek dari simpul 0 ke setiap simpul lainnya.
 */
void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n,
              MPI_Comm comm) {

    int i, loc_v, loc_u, glbl_u, new_dist, my_rank, dist_glbl_u;
    int *loc_known;
    int my_min[2];
    int glbl_min[2];

    MPI_Comm_rank(comm, &my_rank);
    loc_known = malloc(loc_n * sizeof(int));

    Dijkstra_Init(loc_mat, loc_pred, loc_dist, loc_known, my_rank, loc_n);

    /* 
        Membuat loop sebanyak n-1 kali karena jarak terpendek node 0 ke 0 adalah 0
    */
    for (i = 0; i < n - 1; i++) {
        loc_u = Find_min_dist(loc_dist, loc_known, loc_n);

        if (loc_u != -1) {
            my_min[0] = loc_dist[loc_u];
            my_min[1] = loc_u + my_rank * loc_n;
        }
        else {
            my_min[0] = INFINITY;
            my_min[1] = -1;
        }

        /* 
            Menyimpan jarak terpendek global dan nodenya dalam glbl_min
        */
        MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC, comm);

        dist_glbl_u = glbl_min[0];
        glbl_u = glbl_min[1];

        if (glbl_u == -1)
            break;

        if ((glbl_u / loc_n) == my_rank) {
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






/*
    Mencari jarak minimum lokal dari titik 0 ke sebuah simpul
 */
int Find_min_dist(int loc_dist[], int loc_known[], int loc_n) {
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

/*
    Menuliskan setiap jarak dari titik 0 ke setiap simpul ke out.txt
 */
void Print_dists(int global_dist[], int n) {
    int v;
    FILE *fptr;
    fptr = fopen("out.txt","w");
    if(fptr == NULL){
       printf("Error!");   
       exit(1);             
    }

    fprintf(fptr,"  v    dist 0->v\n");
    fprintf(fptr,"----   ---------\n");

    for (v = 1; v < n; v++) {
        if (global_dist[v] == INFINITY) {
            fprintf(fptr,"%3d       %5s\n", v, "inf");
        }
        else
            fprintf(fptr,"%3d       %4d\n", v, global_dist[v]);
    }
    fprintf(fptr,"\n");
    fclose(fptr);
}