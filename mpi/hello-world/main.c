#include "mpi.h"
#include "stdio.h"

int rank;
int nproc;

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Process: %d out of %d:\nHello World!\n", rank, nproc);
    MPI_Finalize();
}
