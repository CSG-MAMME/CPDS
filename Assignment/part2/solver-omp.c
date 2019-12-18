#include "heat.h"

#define NB 8

#define min(a,b) ( ((a) < (b)) ? (a) : (b) )

/*
 * Blocked Jacobi solver: one iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
    double diff, sum=0.0;
    int nbx, bx, nby, by;
  
    nbx = NB; // CSG: sizex = nbx = num_threads
    bx = sizex/nbx;
    nby = NB;
    by = sizey/nby;
    // CSG: The first two loops iter block-wise
    // CSG: We will use a parallel for, parallelizing the first loop
    // hence it seems optimal to use OMP_NUMTHREADS = 4
    // CSG: A small note on the variables types:
    //  - diff: we only need it to compute the LOCAL difference, hence
    //          it does not need to be shared (private).
    //  - sum:  variable we do the reduction on.
    #pragma omp parallel for reduction(+:sum) private(diff)
    for (int ii=0; ii<nbx; ii++)
    {
        for (int jj=0; jj<nby; jj++) 
        {
            // CSG: The second two inner loops iter within one block 
            for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) 
                {
                    // CSG: Average of the four neighbours
                    utmp[i*sizey+j] = 0.25 * (u[i*sizey + (j-1)] 
                        +  u[i*sizey + (j+1)]
                        +  u[(i-1)*sizey + j] 
                        +  u[(i+1)*sizey + j]); 
                    diff = utmp[i*sizey+j] - u[i*sizey + j];
                    sum += diff * diff; 
	        }
       }
    }
    // CSG: Note that there is an implicit barrier at the end of the
    // parallel for region.
    return sum;
}

/*
 * Blocked Red-Black solver: one iteration step
 */
double relax_redblack (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;
    int nbx, bx, nby, by;
    int lsw;

    nbx = NB;
    bx = sizex/nbx;
    nby = NB;
    by = sizey/nby;
    // CSG: This is working well as it is, but could we merge both loops?
    // Computing "Red" blocks
    // CSG: one or two parallel regions? 
    // one pragma omp parallel and two fors
    // CSG: nowait or not
    #pragma omp parallel for reduction(+:sum) private(diff, unew, lsw)
    for (int ii=0; ii<nbx; ii++) {
        lsw = ii%2;
        for (int jj=lsw; jj<nby; jj=jj+2) 
            for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	            unew= 0.25 * (    u[ i*sizey	+ (j-1) ]+  // left
				      u[ i*sizey	+ (j+1) ]+  // right
				      u[ (i-1)*sizey	+ j     ]+  // top
				      u[ (i+1)*sizey	+ j     ]); // bottom
	            diff = unew - u[i*sizey+ j];
	            sum += diff * diff; 
	            u[i*sizey+j]=unew;
	        }
    }

    // Computing "Black" blocks
    #pragma omp parallel for reduction(+:sum) private(diff, unew, lsw)
    for (int ii=0; ii<nbx; ii++) {
        lsw = (ii+1)%2;
        for (int jj=lsw; jj<nby; jj=jj+2) 
            for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	            unew= 0.25 * (    u[ i*sizey	+ (j-1) ]+  // left
				      u[ i*sizey	+ (j+1) ]+  // right
				      u[ (i-1)*sizey	+ j     ]+  // top
				      u[ (i+1)*sizey	+ j     ]); // bottom
	            diff = unew - u[i*sizey+ j];
	            sum += diff * diff; 
	            u[i*sizey+j]=unew;
	        }
    }

    return sum;
}

/*
 * Blocked Gauss-Seidel solver: one iteration step
 */
double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;
    int nbx, bx, nby, by;
    // CSG: Define auxiliary block to compute and keep track of dependencies
    int *block;
    block = (int *)calloc(NB * NB, sizeof(int));

    nbx = NB;
    bx = sizex/nbx;
    nby = NB;
    by = sizey/nby;
    // CSG: Each block has a dependency from the two previous ones
    // (left and up) so i-1 and j-1
    // CSG: To keep track of the dependencies, we make use an auxiliary
    // matrix.

    #pragma omp parallel
    #pragma omp single
    for (int ii=0; ii<nbx; ii++)
    {
        for (int jj=0; jj<nby; jj++) 
        {
            if (ii > 0 && jj > 0)
            {
                #pragma omp task private(unew, diff) \
                    depend(in: block[(ii-1)*nbx + jj], block[ii*nbx + (jj-1)]) \
                    depend(out: block[ii*nbx + jj])
                {
                    for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                    {
                        for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) 
                        {
                            {
                                unew= 0.25 * (u[i*sizey + (j-1)] + 
                                      u[i*sizey + (j+1)] +  
                                      u[(i-1)*sizey + j] +
                                      u[(i+1)*sizey + j]);
                                diff = unew - u[i*sizey+ j];
                                sum += diff * diff; 
                                u[i*sizey+j]=unew;
                            }
                        }
                    }
                }
            }
            else if (ii > 0)
            {
                #pragma omp task private(unew, diff) \
                    depend(in: block[(ii-1)*nbx + jj]) \
                    depend(out: block[ii*nbx + jj])
                {
                    for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                    {
                        for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) 
                        {
                            {
                                unew= 0.25 * (u[i*sizey + (j-1)] + 
                                      u[i*sizey + (j+1)] +  
                                      u[(i-1)*sizey + j] +
                                      u[(i+1)*sizey + j]);
                                diff = unew - u[i*sizey+ j];
                                sum += diff * diff; 
                                u[i*sizey+j]=unew;
                            }
                        }
                    }
                }
            }
            else if (jj > 0)
            {
                #pragma omp task private(unew, diff) \
                    depend(in: block[ii*nbx + (jj-1)]) \
                    depend(out: block[ii*nbx + jj])
                {
                    for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                    {
                        for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) 
                        {
                            {
                                unew= 0.25 * (u[i*sizey + (j-1)] + 
                                      u[i*sizey + (j+1)] +  
                                      u[(i-1)*sizey + j] +
                                      u[(i+1)*sizey + j]);
                                diff = unew - u[i*sizey+ j];
                                sum += diff * diff; 
                                u[i*sizey+j]=unew;
                            }
                        }
                    }
                }
            }
            else 
            {
                #pragma omp task private(unew, diff) \
                    depend(out: block[ii*nbx + jj])
                {
                    for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
                    {
                        for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) 
                        {
                            {
                                unew= 0.25 * (u[i*sizey + (j-1)] + 
                                      u[i*sizey + (j+1)] +  
                                      u[(i-1)*sizey + j] +
                                      u[(i+1)*sizey + j]);
                                diff = unew - u[i*sizey+ j];
                                sum += diff * diff; 
                                u[i*sizey+j]=unew;
                            }
                        }
                    }
                }
            }
        }
    }
    //#pragma omp taskwait

    return sum;
}

