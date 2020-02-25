#include <stdio.h>
#include <stdlib.h>
#include </usr/include/mpi/mpi.h>
#include <time.h>
#include "func.h"
#define INFINITY 1000000

int main(int argc, char **argv) {
    // Inisialisasi setiap variabel yang digunakan
    int *loc_mat, *loc_dist, *loc_pred, *global_dist = NULL, *global_pred = NULL;
    int my_rank, p, loc_n, n;
    double *time_mat;
    MPI_Comm comm;
    MPI_Datatype blk_col_mpi_t;
    clock_t start, end;
    double cpu_time_used;
    // Inisiasi openmpi dan size dari variabel yang digunakan
    MPI_Init(NULL, NULL);
    comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &p);
    n = Read_n(my_rank, comm);
    loc_n = n / p;
    loc_mat = malloc(n * loc_n * sizeof(int));
    loc_dist = malloc(loc_n * sizeof(int));
    loc_pred = malloc(loc_n * sizeof(int));
    blk_col_mpi_t = Build_blk_col_type(n, loc_n);

    if (my_rank == 0) {
        global_dist = malloc(n * sizeof(int));
        global_pred = malloc(n * sizeof(int));
    }

    // Membuat graf dalam bentuk matriks
    Make_matrix(loc_mat, n, loc_n, blk_col_mpi_t, my_rank, comm);
    // Memulai menghitung waktu
    start = clock();
    // Memulai dijkstra untuk setiap node
    Dijkstra(loc_mat, loc_dist, loc_pred, loc_n, n, comm);

    /* Mengumpulkan hasil dari setiap node */
    MPI_Gather(loc_dist, loc_n, MPI_INT, global_dist, loc_n, MPI_INT, 0, comm);
    MPI_Gather(loc_pred, loc_n, MPI_INT, global_pred, loc_n, MPI_INT, 0, comm);

    /* Menuliskan hasil ke dalam file out.txt */
    if (my_rank == 0) {
        Print_dists(global_dist, n);
        free(global_dist);
        free(global_pred);
    }
    free(loc_mat);
    free(loc_pred);
    free(loc_dist);
    MPI_Type_free(&blk_col_mpi_t);
    MPI_Finalize();
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("Process in node %d use time : %f\n",my_rank,cpu_time_used);
    return 0;
}