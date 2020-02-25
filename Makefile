all: mpi_dijkstra

mpi_dijkstra: mpi_dijkstra.o func.o
	mpicc mpi_dijkstra.o func.o -o mpi_dijkstra

mpi.o: mpi_dijkstra.c
	mpicc -c mpi_dijkstra.c -o mpi_dijkstra.o

func.o: func.c
	mpicc -c func.c -o func.o

clean:
	rm -f mpi_dijkstra.o func.o mpi_dijkstra core *~