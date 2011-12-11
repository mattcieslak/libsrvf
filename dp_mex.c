#include "mex.h"
#include <math.h>
#include "dp_grid.h"


/* For C89 compilers */
static int _round( double v ){
  double f, r;

  f = floor(v);
  r = v-f;
  
  if ( r <= 0.5 ) return (int)f;
  else            return (int)ceil(v);
}

/* Signature:
 * function [G T] = dp_mex( Q1, T1, Q2, T2, ndp )
 * Arguments are checked in dp.m, not here.  */
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]){
  double *Q1 = 0;
  double *T1 = 0;
  double *Q2 = 0;
  double *T2 = 0;
  int nsamps1;
  int nsamps2;
  double *G = 0;
  double *T = 0;
  int Gsize;
  int ndp = 0;
  int dim = 0;
  double *E = 0; /* E[ndp*i+j] = cost of best path to (i,j) */
  int *P = 0; /* P[ndp*i+j] = predecessor of (i,j) along best path */
  double m, rootm;
  int r1, r2, c1, c2;
  int i, j, k;

  Q1 = mxGetPr( prhs[0] );
  T1 = mxGetPr( prhs[1] );
  Q2 = mxGetPr( prhs[2] );
  T2 = mxGetPr( prhs[3] );
  ndp = (int)mxGetScalar( prhs[4] );
  dim = mxGetM( prhs[0] );
  nsamps1 = mxGetN( prhs[1] ); /* = columns(T1) = columns(Q1)+1 */
  nsamps2 = mxGetN( prhs[3] ); /* = columns(T2) = columns(Q2)+1 */


  if ( !(E=(double*)mxMalloc(ndp*ndp*sizeof(double))) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate E" );
    goto cleanup;
  }
  if ( !(P=(int*)mxCalloc(ndp*ndp,sizeof(int))) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate P" );
    goto cleanup;
  }
  if ( !(plhs[0]=mxCreateDoubleMatrix(1,ndp,mxREAL)) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleMatrix failed" );
    goto cleanup;
  }
  if ( !(plhs[1]=mxCreateDoubleMatrix(1,ndp,mxREAL)) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleMatrix failed" );
    goto cleanup;
  }

  G = mxGetPr( plhs[0] );
  T = mxGetPr( plhs[1] );

  /* Compute cost of best path from (0,0) to every other grid point */
  dp_costs( Q1, T1, nsamps1, Q2, T2, nsamps2, dim, E, P, ndp, ndp );

  /* Reconstruct best path from (0,0) to (1,1) */
  Gsize = dp_build_gamma( P, ndp, ndp, G, T );
  mxSetN( plhs[0], Gsize );
  mxSetN( plhs[1], Gsize );

cleanup:
  if ( E ) mxFree( E );
  if ( P ) mxFree( P );
}
