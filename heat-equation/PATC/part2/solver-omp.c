#include "heat.h"
#include <omp.h>

#define min(a,b) ( ((a) < (b)) ? (a) : (b) )

/*
 * One Jacobi iteration step
 */
double relax_jacobi (double *u, double *utmp, unsigned sizex, unsigned sizey)
{
    double diff, sum=0.0;
    int nbx, bx, nby, by;
  
    nbx=omp_get_max_threads();
    bx = sizex/nbx;
    nby = nbx;
    by = sizey/nby;
#pragma omp parallel for schedule(static) private(diff) reduction(+:sum)
    for (int ii=0; ii<nbx; ii++)
    for (int jj=0; jj<nby; jj++) 
    for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
        for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	    utmp[i*sizey+j]= 0.25 * (u[ i*sizey     + (j-1) ]+  // left
				     u[ i*sizey     + (j+1) ]+  // right
				     u[ (i-1)*sizey + j     ]+  // top
				     u[ (i+1)*sizey + j     ]); // bottom
	    diff = utmp[i*sizey+j] - u[i*sizey + j];
	    sum += diff * diff; 
	}
  
    return sum;
}

/*
 * One Red-Black iteration step
 */
double relax_redblack (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;
    int nbx, bx, nby, by;
    int lsw;

    nbx=omp_get_max_threads();
    bx = sizex/nbx;
    nby = nbx;
    by = sizey/nby;
    //lsw = 0;
#pragma omp parallel for private(lsw,diff,unew) reduction(+:sum)
    for (int ii=0; ii<nbx; ii++) {
    lsw = ii%2;
    for (int jj=lsw; jj<nby; jj=jj+2) 
    for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
        for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	    unew= 0.25 * (	  u[ i*sizey	+ (j-1) ]+  // left
				  u[ i*sizey	+ (j+1) ]+  // right
				  u[ (i-1)*sizey	+ j     ]+  // top
				  u[ (i+1)*sizey	+ j     ]); // bottom
	    diff = unew - u[i*sizey+ j];
	    sum += diff * diff; 
	    u[i*sizey+j]=unew;
	}
    //lsw = 1 - lsw;
    }

    //lsw = 1;
#pragma omp parallel for private(lsw,diff,unew) reduction(+:sum)
    for (int ii=0; ii<nbx; ii++) {
    lsw = (ii+1)%2;
    for (int jj=lsw; jj<nby; jj=jj+2) 
    for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
        for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	    unew= 0.25 * (	  u[ i*sizey	+ (j-1) ]+  // left
				  u[ i*sizey	+ (j+1) ]+  // right
				  u[ (i-1)*sizey	+ j     ]+  // top
				  u[ (i+1)*sizey	+ j     ]); // bottom
	    diff = unew - u[i*sizey+ j];
	    sum += diff * diff; 
	    u[i*sizey+j]=unew;
	}
    //lsw = 1 - lsw;
    }

    return sum;
}

/*
 * One Gauss-Seidel iteration step
 */

double relax_gauss (double *u, unsigned sizex, unsigned sizey)
{
    double unew, diff, sum=0.0;
    int nbx, bx, nby, by;
    int blockfinished[24]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    nbx=omp_get_max_threads();
    bx = sizex/nbx;
    nby = nbx;
    by = sizey/nby;
#pragma omp parallel for schedule(static,1) private(diff,unew) reduction(+:sum)
    for (int ii=0; ii<nbx; ii++)
    for (int jj=0; jj<nby; jj++) {
    if (ii>0) {
	while(blockfinished[ii-1]<=jj) {
		#pragma omp flush(blockfinished)
	}
    }
    for (int i=1+ii*bx; i<=min((ii+1)*bx, sizex-2); i++) 
        for (int j=1+jj*by; j<=min((jj+1)*by, sizey-2); j++) {
	    unew= 0.25 * (	  u[ i*sizey	+ (j-1) ]+  // left
				  u[ i*sizey	+ (j+1) ]+  // right
				  u[ (i-1)*sizey	+ j     ]+  // top
				  u[ (i+1)*sizey	+ j     ]); // bottom
	    diff = unew - u[i*sizey+ j];
	    sum += diff * diff; 
	    u[i*sizey+j]=unew;
	}
    blockfinished[ii]++;
    #pragma omp flush(blockfinished)
    }

    return sum;
}

