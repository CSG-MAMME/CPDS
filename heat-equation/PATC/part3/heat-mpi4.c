/*
 * Iterative solver for heat distribution
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "heat.h"

void usage( char *s )
{
    fprintf(stderr, 
	    "Usage: %s <input file> [result file]\n\n", s);
}

int main( int argc, char *argv[] )
{
    unsigned iter;
    FILE *infile, *resfile;
    char *resfilename;
    int myid, numprocs;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (myid == 0) {
      printf("I am the master (%d) and going to distribute work to %d additional workers ...\n", myid, numprocs-1);

    // algorithmic parameters
    algoparam_t param;
    int np;

    double runtime, flop;
    double residual=0.0;

    // check arguments
    if( argc < 2 )
    {
	usage( argv[0] );
	return 1;
    }

    // check input file
    if( !(infile=fopen(argv[1], "r"))  ) 
    {
	fprintf(stderr, 
		"\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);
      
	usage(argv[0]);
	return 1;
    }

    // check result file
    resfilename= (argc>=3) ? argv[2]:"heat.ppm";

    if( !(resfile=fopen(resfilename, "w")) )
    {
	fprintf(stderr, 
		"\nError: Cannot open \"%s\" for writing.\n\n", 
		resfilename);
	usage(argv[0]);
	return 1;
    }

    // check input
    if( !read_input(infile, &param) )
    {
	fprintf(stderr, "\nError: Error parsing input file.\n\n");
	usage(argv[0]);
	return 1;
    }
    print_params(&param);

    // set the visualization resolution
    
    param.u     = 0;
    param.uhelp = 0;
    param.uvis  = 0;
    param.visres = param.resolution;
   
    if( !initialize(&param) )
	{
	    fprintf(stderr, "Error in Solver initialization.\n\n");
	    usage(argv[0]);
            return 1;
	}

    // full size (param.resolution are only the inner points)
    np = param.resolution + 2;
    
    // starting time
    runtime = wtime();

    // send to workers the necessary data to perform computation
    int first_row = 1;
    int last_row;
    for (int i=0; i<numprocs; i++) {
        last_row = first_row + param.resolution/numprocs - 1;
        last_row = (i < (param.resolution%numprocs)) ? last_row+1 : last_row;
        if (i>0) {
                MPI_Send(&param.maxiter, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&param.resolution, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&param.algorithm, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&first_row, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&last_row, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&param.u[(first_row-1)*np], (last_row-first_row+3)*np, 
			 MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(&param.uhelp[(first_row-1)*np], (last_row-first_row+3)*np, 
			 MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
        }
        first_row = last_row + 1;
    } 

    first_row = 1;
    last_row = first_row + param.resolution/numprocs - 1;
    last_row = ((param.resolution%numprocs)>0) ? last_row+1 : last_row;
    int rows = last_row-first_row+1;

#if 1
    iter = 0;
    while(1) {
	switch( param.algorithm ) {
	    case 0: // JACOBI
	            residual = relax_jacobi(param.u, param.uhelp, rows+2, np);
		    MPI_Sendrecv(&param.uhelp[last_row*np], np, MPI_DOUBLE, 1, 0,
		     		 &param.uhelp[(last_row+1)*np], np, MPI_DOUBLE, 1, 0,
		     		 MPI_COMM_WORLD, &status);
		    // Copy uhelp into u
		    for (int i=first_row-1; i<last_row+2; i++)
    		        for (int j=0; j<np; j++)
	    		    param.u[ i*np+j ] = param.uhelp[ i*np+j ];
		    break;
	    case 1: // RED-BLACK
		    residual = relax_redblack(param.u, np, np);
		    break;
	    case 2: // GAUSS
		    residual = relax_gauss(param.u, np, np);
		    break;
	    }

        iter++;

        // solution good enough ?
        double gresidual=0.0;
        MPI_Reduce(&residual, &gresidual, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&gresidual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (gresidual < 0.00005) break;

        // max. iteration reached ? (no limit with maxiter=0)
        if (param.maxiter>0 && iter>=param.maxiter) break;
    }
#endif

    // receive from workers the portions of matrix u computed
    first_row = 1;
    for (int i=0; i<numprocs; i++) {
        last_row = first_row + param.resolution/numprocs - 1;
        last_row = (i < (param.resolution%numprocs)) ? last_row+1 : last_row;
	if (i>0) MPI_Recv(&param.u[first_row*np], (last_row-first_row+1)*(np), 
		          MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
        first_row = last_row + 1;
    } 

    // Flop count after iter iterations
    flop = iter * 11.0 * param.resolution * param.resolution;
    // stopping time
    runtime = wtime() - runtime;

    fprintf(stdout, "Time: %04.3f ", runtime);
    fprintf(stdout, "(%3.3f GFlop => %6.2f MFlop/s)\n", 
	    flop/1000000000.0,
	    flop/runtime/1000000);
    fprintf(stdout, "Convergence to residual=%f: %d iterations\n", residual, iter);

    // for plot...
    coarsen( param.u, np, np,
	     param.uvis, param.visres+2, param.visres+2 );
  
    write_image( resfile, param.uvis,  
		 param.visres+2, 
		 param.visres+2 );

    finalize( &param );

    fprintf(stdout, "Process %d finished computing with residual value = %f\n", myid, residual);

    MPI_Finalize();

    return 0;

} else {

    printf("I am worker %d and ready to receive work to do ...\n", myid);

    // receive information from master to perform computation locally

    int columns, rows, np;
    int first_row, last_row;
    int iter, maxiter;
    int algorithm;
    double residual;

    MPI_Recv(&maxiter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&columns, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&algorithm, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&first_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&last_row, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);

    printf("Worker %d will receive rows from %d to %d\n", myid, first_row, last_row);

    rows = last_row-first_row+1;
    np = columns + 2;

    // allocate memory for worker
    double * u = calloc( sizeof(double), (rows+2)*np );
    double * uhelp = calloc( sizeof(double), (rows+2)*np );
    if( (!u) || (!uhelp) )
    {
        fprintf(stderr, "Error: Cannot allocate memory\n");
        return 0;
    }
    
    // fill initial values for matrix with values received from master
    MPI_Recv(u, (rows+2)*np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(uhelp, (rows+2)*np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

#if 1
    iter = 0;
    while(1) {
	switch( algorithm ) {
	    case 0: // JACOBI
	            residual = relax_jacobi(u, uhelp, rows+2, np);
		    if (myid < (numprocs-1)) {
			MPI_Sendrecv(&uhelp[rows*np], np, MPI_DOUBLE, myid+1, 0,
				     &uhelp[(rows+1)*np], np, MPI_DOUBLE, myid+1, 0,
			   	     MPI_COMM_WORLD, &status);
			}
	 	    MPI_Sendrecv(&uhelp[np], np, MPI_DOUBLE, myid-1, 0,
		     		 &uhelp[0], np, MPI_DOUBLE, myid-1, 0,
		     		 MPI_COMM_WORLD, &status);
		    // Copy uhelp into u
		    for (int i=0; i<rows+2; i++)
    		        for (int j=0; j<np; j++)
	    		    u[ i*np+j ] = uhelp[ i*np+j ];
		    break;
	    case 1: // RED-BLACK
		    residual = relax_redblack(u, np, np);
		    break;
	    case 2: // GAUSS
		    residual = relax_gauss(u, np, np);
		    break;
	    }

        iter++;

        // solution good enough ?
        MPI_Reduce(&residual, 0, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Bcast(&residual, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (residual < 0.00005) break;

        // max. iteration reached ? (no limit with maxiter=0)
        if (maxiter>0 && iter>=maxiter) break;
    }
#endif

    // send values computed for matrix to the master
    MPI_Send(&u[np], rows*np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

    if( u ) free(u); if( uhelp ) free(uhelp);

    fprintf(stdout, "Process %d finished computing %d iterations with residual value = %f\n", myid, iter, residual);

    MPI_Finalize();
    exit(0);
  }
}
