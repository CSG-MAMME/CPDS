EXECS=mpi_blocking_comms
MPICC?=mpicc

all: ${EXECS}

mpi_blocking_comms: main.c
	${MPICC} -o mpi_blocking_comms main.c

clean:
	rm -f ${EXECS}
