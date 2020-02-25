#include <stdio.h>
#include <stdlib.h>
#include </usr/include/mpi/mpi.h>
#include <time.h>

#ifndef FUNC_H
#define FUNC_H

int Read_n(int my_rank, MPI_Comm comm);
MPI_Datatype Build_blk_col_type(int n, int loc_n);
void Make_matrix(int loc_mat[], int n, int loc_n, MPI_Datatype blk_col_mpi_t,
                 int my_rank, MPI_Comm comm);
void Dijkstra_Init(int loc_mat[], int loc_pred[], int loc_dist[], int loc_known[],
                   int my_rank, int loc_n);
void Dijkstra(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int n,
              MPI_Comm comm);
int Find_min_dist(int loc_dist[], int loc_known[], int loc_n);
void Print_dists(int global_dist[], int n);

#endif