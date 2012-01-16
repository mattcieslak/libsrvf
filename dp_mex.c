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
 * function [G T] = dp_mex( Q1, T1, Q2, T2, tv1, tv2 )
 * Arguments are checked in dp.m, not here.  */
void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]){
  double *Q1 = 0;
  double *T1 = 0;
  double *Q2 = 0;
  double *T2 = 0;
  int nsamps1;
  int nsamps2;
  double *tv1 = 0;
  double *tv2 = 0;
  int ntv1;
  int ntv2;
  double *G = 0;
  double *T = 0;
  int Gsize;
  int dim = 0;
  double *E = 0; /* E[ntv1*j+i] = cost of best path to (tv1[i],tv2[j]) */
  int *P = 0; /* P[ntv1*j+i] = predecessor of (tv1[i],tv2[j]) along best path */
  double m, rootm;
  int sr, sc; /* source row and column index */
  int tr, tc; /* target row and column index */
  int i, j, k;
  int Galloc_size;
  double *pres;

  Q1 = mxGetPr( prhs[0] );
  T1 = mxGetPr( prhs[1] );
  Q2 = mxGetPr( prhs[2] );
  T2 = mxGetPr( prhs[3] );
  tv1 = mxGetPr( prhs[4] );
  tv2 = mxGetPr( prhs[5] );

  dim = mxGetM( prhs[0] );
  nsamps1 = mxGetN( prhs[1] ); /* = columns(T1) = columns(Q1)+1 */
  nsamps2 = mxGetN( prhs[3] ); /* = columns(T2) = columns(Q2)+1 */
  ntv1 = mxGetN( prhs[4] );
  ntv2 = mxGetN( prhs[5] );
  Galloc_size = ntv1>ntv2 ? ntv1 : ntv2;


  if ( !(E=(double*)mxMalloc(ntv1*ntv2*sizeof(double))) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate E" );
    goto cleanup;
  }
  if ( !(P=(int*)mxCalloc(ntv1*ntv2,sizeof(int))) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate P" );
    goto cleanup;
  }
  if ( !(plhs[0]=mxCreateDoubleMatrix(1,Galloc_size,mxREAL)) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleMatrix failed" );
    goto cleanup;
  }
  if ( !(plhs[1]=mxCreateDoubleMatrix(1,Galloc_size,mxREAL)) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleMatrix failed" );
    goto cleanup;
  }
  if ( !(plhs[2]=mxCreateDoubleScalar(0.0)) )
  { 
    mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleScalar failed" );
    goto cleanup;
  }

  G = mxGetPr( plhs[0] );
  T = mxGetPr( plhs[1] );
  pres = mxGetPr( plhs[2] );

  /* Compute cost of best path from (0,0) to every other grid point */
  *pres = dp_costs( Q1, T1, nsamps1, Q2, T2, nsamps2, 
                    dim, tv1, ntv1, tv2, ntv2, E, P );

  /* Reconstruct best path from (0,0) to (1,1) */
  Gsize = dp_build_gamma( P, tv1, ntv1, tv2, ntv2, G, T );
  mxSetN( plhs[0], Gsize );
  mxSetN( plhs[1], Gsize );

cleanup:
  if ( E ) mxFree( E );
  if ( P ) mxFree( P );
}
