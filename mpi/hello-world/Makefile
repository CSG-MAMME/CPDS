EXECS=mpi_hello_world
MPICC?=mpicc

all: ${EXECS}

mpi_hello_world: main.c
	${MPICC} -o mpi_hello_world main.c

clean:
	rm -f ${EXECS}
