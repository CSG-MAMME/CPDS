#include "mpi.h"
#include "stdio.h"

int rank;
int nproc;

int main(int argc, char* argv[])
{
    int i, myid, numprocs, num_steps;
    double x, pi, step, sum = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 0) scanf("%d\n", 

    if (rank == 0)
    {
        isbuf = 1;
        printf("Before sending... \n");
        MPI_Send(&isbuf, 1, MPI_INTEGER, 1, 1, MPI_COMM_WORLD);
        printf("Sent %d!\n", isbuf);
    }
    else if (rank == 1)
    {
        printf("Before receiving...\n");
        MPI_Recv(&irbuf, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, &status);
        printf("Received %d!\n", irbuf);
    }

    MPI_Finalize();
}
