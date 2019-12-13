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
    int myid, numprocs, myrows;
    double global_residual;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    fprintf(stderr, "Hello error\n");

    if (myid == 0) 
    {
        printf("I am the master (%d) and going to distribute work to %d additional workers ...\n", myid, numprocs-1);

        // algorithmic parameters
        algoparam_t param;
        int np;
        double runtime, flop;
        double residual = 0.0;

        // check arguments
        if(argc < 2)
        {
            usage(argv[0]);
            return 1;
        }

        // check input file
        if(!(infile=fopen(argv[1], "r"))) 
        {
            fprintf(stderr, 
                "\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);
            usage(argv[0]);
            return 1;
        }

        // check result file
        resfilename = (argc>=3) ? argv[2] : "heat.ppm";

        if(!(resfile=fopen(resfilename, "w")))
        {
            fprintf(stderr, 
                "\nError: Cannot open \"%s\" for writing.\n\n", 
                resfilename);
            usage(argv[0]);
            return 1;
        }

        // check input
        if(!read_input(infile, &param))
        {
            fprintf(stderr, "\nError: Error parsing input file.\n\n");
            usage(argv[0]);
            return 1;
        }
        printf("MPI Num Procs: %d\n", numprocs);
        print_params(&param);

        // set the visualization resolution
        
        param.u     = 0;
        param.uhelp = 0;
        param.uvis  = 0;
        param.visres = param.resolution;
       
        if(!initialize(&param))
        {
            fprintf(stderr, "Error in Solver initialization.\n\n");
            usage(argv[0]);
                return 1;
        }

        fprintf(stderr, "Master Initialized\n");
        // full size (param.resolution are only the inner points)
        np = param.resolution + 2;
        // CSG: Define these two helper parameters
        myrows = param.resolution / numprocs;
        fprintf(stderr, "np: %i myrows: %i\n", np, myrows);
        
        // starting time
        runtime = wtime();

        // send to workers the necessary data to perform computation
        for (int i = 0; i < numprocs; i++)
        {
            if (i > 0) 
            {
                MPI_Send(&param.maxiter, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&param.resolution, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&param.algorithm, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                //MPI_Send(&param.u[0], (np)*(np), MPI_DOUBLE, i, 0,
                //         MPI_COMM_WORLD);
                //MPI_Send(&param.uhelp[0], (np)*(np), MPI_DOUBLE, i, 0, 
                //         MPI_COMM_WORLD);
                // CSG: Send the right data to everyone
                fprintf(stderr, "Master Distributing: %d\n", i);
                fprintf(stderr, "Accessing Value: %d\n", myrows * np * i);
                MPI_Send(&param.u[myrows * np * i], (myrows + 2) * (np),
                         MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
                MPI_Send(&param.uhelp[myrows * np * i], (myrows + 2) * (np),
                         MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
            }
        }

        fprintf(stderr, "Master Distributed\n");

        iter = 0;
        while(1) {
        switch( param.algorithm ) {
            case 0: // JACOBI
                // CSG: We need to change the sizex to match new resolution
                residual = relax_jacobi(param.u, param.uhelp, myrows + 2, np);
                // Copy uhelp into u
                for (int i = 0; i < (myrows + 2); i++)
                    for (int j = 0; j < np; j++)
                        param.u[i * np + j] = param.uhelp[i * np + j];
                // CSG: We need to exchange our lower boudnry (myid == 0)
                // with the next process. Down => Send then Recv
                // CSG: this check is to prevent deadlocks when running w/
                // only one process.
                if (myid != numprocs - 1)
                {
                    MPI_Send(&param.u[myrows * np], np, MPI_DOUBLE, myid + 1,
                             iter, MPI_COMM_WORLD);
                    MPI_Recv(&param.u[(myrows + 1) * np], np, MPI_DOUBLE,
                             myid + 1, iter, MPI_COMM_WORLD, &status);
                }
                break;
            case 1: // RED-BLACK
                residual = relax_redblack(param.u, np, np);
                break;
            case 2: // GAUSS
                // CSG: In this case, no need to recv before computin,
                // first row is dependency free.
                residual = relax_gauss(param.u, myrows + 2, np);
                if (myid != numprocs -1)
                {
                    MPI_Send(&param.u[myrows * np], np, MPI_DOUBLE, myid + 1,
                             iter, MPI_COMM_WORLD);
                    MPI_Recv(&param.u[(myrows + 1) * np], np, MPI_DOUBLE,
                             myid + 1, iter, MPI_COMM_WORLD, &status);
                }
                break;
            }

            iter++;

            // CSG: We now have to compute a global residual. We reduce in the
            // master and broadcast the value.
            // MPI Reduce + MPI Broadcast <=> MPI_Allreduce
            // https://mpitutorial.com/tutorials/mpi-reduce-and-allreduce/
            switch  (param.algorithm)
            {
                case 0: // Jacobi
                case 2: // Gauss
                    MPI_Allreduce(&residual, &global_residual, 1, MPI_DOUBLE,
                                  MPI_SUM, MPI_COMM_WORLD);
                    break;
                //case 2: // CSG: Broken version
                //
                //    MPI_Reduce(&residual, &global_residual, 1, MPI_DOUBLE,
                //                  MPI_SUM, 0, MPI_COMM_WORLD);
                //    MPI_Bcast(&global_residual, 1, MPI_DOUBLE, 0,
                //              MPI_COMM_WORLD);
                //    break;
                default:
                    global_residual = residual;
                    break;
            }

            // solution good enough ?
            //if (residual < 0.00005) break;
            // CSG: Use the global value
            if (global_residual < 0.00005) break;

            // max. iteration reached ? (no limit with maxiter=0)
            if (param.maxiter>0 && iter>=param.maxiter) break;
        }

        // CSG: Receive results
        for (int i = 1; i < numprocs; i++)
        {
            MPI_Recv(&param.u[(myrows * i + 1) * np], myrows * np, MPI_DOUBLE,
                     i, param.maxiter + 1, MPI_COMM_WORLD, &status);
        }

        // Flop count after iter iterations
        flop = iter * 11.0 * param.resolution * param.resolution;
        // stopping time
        runtime = wtime() - runtime;

        fprintf(stdout, "Time: %04.3f ", runtime);
        fprintf(stdout, "(%3.3f GFlop => %6.2f MFlop/s)\n", 
                flop/1000000000.0,
                flop/runtime/1000000);
        fprintf(stdout, "Convergence to residual=%f: %d iterations\n",
                global_residual, iter);

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
        int iter, maxiter;
        int algorithm;
        double residual;

        // CSG: Missing i
        MPI_Recv(&maxiter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&columns, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(&algorithm, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        printf("Worker %d alive!\n", myid);

        // CSG: we modify the rows to be resolution/numprocs
        // rows = columns;
        // Note that param.resolution (second MPI_Send) matches the "columns"
        // variable, and numprocs is available to all processes.
        myrows = columns / numprocs;
        np = columns + 2;
        //myoffset = myid * myrows * np;

        // allocate memory for worker
        // CSG: Modifying rows, we also alter the memory we allocate.
        double *u = calloc(sizeof(double), (myrows + 2) * np);
        double *uhelp = calloc(sizeof(double), (myrows + 2) * np);
        if( (!u) || (!uhelp) )
        {
            fprintf(stderr, "Error: Cannot allocate memory\n");
            return 0;
        }
        
        // fill initial values for matrix with values received from master
        MPI_Recv(&u[0], (myrows + 2) * np, MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, &status);
        MPI_Recv(&uhelp[0], (myrows + 2) * np, MPI_DOUBLE, 0, 0,
                 MPI_COMM_WORLD, &status);
        iter = 0;
        while(1) {
        switch( algorithm ) {
            case 0: // JACOBI
                // residual = relax_jacobi(u, uhelp, np, np);
                // CSG: We need to change the sizex to match new resolution
                residual = relax_jacobi(u, uhelp, myrows + 2, np);
                // Copy uhelp into u
                // CSG: This was a tricky one to catch
                for (int i = 0; i < (myrows + 2); i++)
                    for (int j = 0; j < np; j++)
                        u[i * np + j] = uhelp[i * np + j];
                // CSG: First we go down. Master process sends then receives.
                // So we receive then send.
                MPI_Recv(&u[0], np, MPI_DOUBLE, myid - 1, iter,
                         MPI_COMM_WORLD, &status);
                if (myid != numprocs -1)
                {
                    MPI_Send(&u[myrows * np], np, MPI_DOUBLE, myid + 1, iter,
                             MPI_COMM_WORLD);
                } 
                // CSG: Second we go up. We start from the bottom so we Send
                // we send our second row (element np)
                MPI_Send(&u[np], np, MPI_DOUBLE, myid - 1, iter,
                         MPI_COMM_WORLD);
                if (myid != numprocs - 1)
                {
                    MPI_Recv(&u[(myrows + 1) * np], np, MPI_DOUBLE, myid + 1,
                             iter, MPI_COMM_WORLD, &status);
                }
                break;
            case 1: // RED-BLACK
                residual = relax_redblack(u, np, np);
                break;
            case 2: // GAUSS
                // CSG: The main difference here is that given that we, 
                // EXPLICITLY modify u at each iteration. We must receive the
                // freshest values of u.
                // CSG: This means, the upper boundry computed at this round,
                // and the lower in the previous (next block will be computed
                // strictly after).
                MPI_Recv(&u[0], np, MPI_DOUBLE, myid - 1, iter,
                         MPI_COMM_WORLD, &status);
                // CSG: receiving from a previous iteration, enables us to
                // cascade parallelism.
                if (iter > 0 && myid != numprocs - 1)
                {
                    MPI_Recv(&u[(myrows + 1) * np], np, MPI_DOUBLE, myid + 1,
                             iter - 1, MPI_COMM_WORLD, &status);
                }
                // CSG: Now we can update the values of our u chunk.
                residual = relax_gauss(u, myrows + 2, np);
                // CSG: Now we need to send the new results down.
                // CSG: It is important to send down first to prevent deadlocks.
                if (myid != numprocs - 1)
                {
                    MPI_Send(&u[myrows * np], np, MPI_DOUBLE, myid + 1, iter,
                             MPI_COMM_WORLD);
                }
                // CSG: and upwards for the next round (we receive them
                // at iteration iter + 1).
                MPI_Send(&u[np], np, MPI_DOUBLE, myid - 1, iter,
                         MPI_COMM_WORLD);
                break;
            }

            iter++;

            // CSG: Global Residual only for Jacobi?
            switch  (algorithm)
            {
                case 0: // Jacobi
                case 2: // Gauss
                    MPI_Allreduce(&residual, &global_residual, 1, MPI_DOUBLE,
                                  MPI_SUM, MPI_COMM_WORLD);
                    break;
//                case 2: // Gauss
//                    MPI_Send(&residual, 1, MPI_DOUBLE, 0, myid,
//                             MPI_COMM_WORLD);
//                    MPI_Bcast(&global_residual, 1, MPI_DOUBLE, 0,
//                              MPI_COMM_WORLD);
//                    break;
                default:
                    global_residual = residual;
                    break;
            }

            // solution good enough ?
            if (global_residual < 0.00005) break;

            // max. iteration reached ? (no limit with maxiter=0)
            if (maxiter > 0 && iter >= maxiter) break;
        }

        // Send results back to master
        fprintf(stderr, "Worker %d sent results %i\n", myid, iter);
        MPI_Send(&u[np], myrows * np, MPI_DOUBLE, 0, maxiter + 1,
                 MPI_COMM_WORLD);

        if( u ) free(u); if( uhelp ) free(uhelp);

        fprintf(stdout, "Process %d finished computing %d iterations with residual value = %f\n", myid, iter, residual);

        MPI_Finalize();
        exit(0);
  }
}
