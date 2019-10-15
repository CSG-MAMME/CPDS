EXECS=mpi_scatter_gather
MPICC?=mpicc

all: ${EXECS}

mpi_scatter_gather: main.c
	${MPICC} -o mpi_scatter_gather main.c -lm

clean:
	rm -f ${EXECS}
