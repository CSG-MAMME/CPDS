#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int rank;
int nproc;

int main(int argc, char* argv[])
{
    int nproc, rank;
    int gsize, localbuff[100];
    int root = 0, *rootbuff;
    unsigned int i,j;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == root)
        rootbuff = (int *) malloc(sizeof(int) * 100 * nproc);

    MPI_Scatter(rootbuff, 100, MPI_INT, localbuff, 100, MPI_INT, root,
        MPI_COMM_WORLD);

    for (i = 0; i < 100; i++)
        localbuff[i] = pow(10, rank) + i;

    MPI_Gather(localbuff, 100, MPI_INT, rootbuff, 100, MPI_INT, root,
        MPI_COMM_WORLD);

    if (rank == root)
    {
        for (j = 0; j < nproc; j++)
        {
            for (i = 0; i < 100; i++)
                printf("M[%i][%i] = %i ", j, i, rootbuff[j*100 + i]);
            printf("\n");
        }
    }

    MPI_Finalize();
}
