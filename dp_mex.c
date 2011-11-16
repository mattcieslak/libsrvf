#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include "dp_nbhd.h"


/* For C89 compilers */
static int _round( double v ){
  double f, r;

  f = floor(v);
  r = v-f;
  
  if ( r <= 0.5 ) return (int)f;
  else            return (int)ceil(v);
}


static void do_usage(){
  mexErrMsgIdAndTxt( 
    "dp:InvalidArguments", 
    "Usage:  gamma = dp( q1, q2, Ndp ), where \n\t%s, \n\t%s, \n\t%s, and \n\t%s.\n", 
    "q1 represents the reference SRV function (one row for each dimension, one column for each sample point)",
    "q2 represents the target SRV function (must be the same size as q1)",
    "Ndp is the size of the DP grid to be used",
    "gamma represents the optimal reparametrization of q2"
  );
}


void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]){
  double *q1 = NULL;      /* reference curve SRV fcn, evaluated at [0:Nq-1]/(Nq-1) */
  double *q2 = NULL;      /* target curve SRV fcn, evaluated at [0:Nq-1]/(Nq-1) */
  double *gamma = NULL;   /* (result) optimal reparametrization for q2 */
  int     Nq;             /* number of sample points in q1, q2, and gamma */
  int     Ndp;            /* width and height of DP grid */
  int     dim;            /* ambient curve dimension */
  double *E = NULL;       /* E[Ndp*i+j] = cost of best path from (0,0) to (row_i,col_j) */
  int    *P = NULL;       /* P[Ndp*i+j] = predecessor of (i,j) on best path from (0,0) to (i,j) */
  int     r1,r2,c1,c2;
  int     i,j,k;
  int     idx1,idx2;
  double  m, rootm, d, dxi;
  double  s, ds;
  double  curE, integralVal;


  if ( nrhs != 3 ){mexErrMsgIdAndTxt("dp:BadArgs","Wrong number of input arguments"); do_usage(); return;}
  if ( nlhs != 1 ){mexErrMsgIdAndTxt("dp:BadArgs","Wrong number of output arguments"); do_usage(); return;}

  /* Unpack arguments and validate them */
  q1 = mxGetPr( prhs[0] );
  q2 = mxGetPr( prhs[1] );
  Ndp = (int)*mxGetPr( prhs[2] );
  dim  = mxGetM( prhs[0] );
  Nq   = mxGetN( prhs[0] );

  if ( Ndp < 2 ){mexErrMsgIdAndTxt("dp:BadArgs","Ndp < 2"); return;}
  if ( Nq  < 2 ){mexErrMsgIdAndTxt("dp:BadArgs","Nq  < 2"); return;}
  if ( dim < 1 ){mexErrMsgIdAndTxt("dp:BadArgs","dim < 1"); return;}
  if ( mxGetM(prhs[1])!=dim ){mexErrMsgIdAndTxt("dp:BadArgs","q1 and q2 have different dimensions"); return;}
  if ( mxGetN(prhs[1])!=Nq  ){mexErrMsgIdAndTxt("dp:BadArgs","q1 and q2 have different numbers of sample points"); return;}

  /* Scale factor */
  s=(double)(Nq-1) / (double)(Ndp-1);
  /* DP square width */
  ds=1.0 / (double)(Ndp-1);

  if ( !(E=(double*)mxMalloc(Ndp*Ndp*sizeof(double))) )
  { mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate E" ); goto cleanup; }
  if ( !(P=(int*)mxCalloc(Ndp*Ndp,sizeof(int))) )
  { mexErrMsgIdAndTxt( "dp:AllocFailed", "failed to allocate P" ); goto cleanup; }
  if ( !(plhs[0]=mxCreateDoubleMatrix(1,Ndp,mxREAL)) )
  { mexErrMsgIdAndTxt( "dp:AllocFailed", "mxCreateDoubleMatrix failed" ); goto cleanup; }

  gamma = mxGetPr( plhs[0] );
  E[0] = 0.;  
  for ( i=1; i<Ndp; i++ ){
    E[i] = E[i*Ndp] = 1e9;
  }

  for ( r2=1; r2<Ndp; r2++ ){
    for ( c2=1; c2<Ndp; c2++ ){
      E[r2*Ndp+c2] = 1E9;
      
      for ( i=0; i<DP_NBHD_COUNT; i++ ){
        r1 = r2-dp_nbhd[i][0];
        c1 = c2-dp_nbhd[i][1];
        if ( r1<0 || c1<0 ) continue;

        m=((double)(r2-r1))/((double)(c2-c1));
        rootm = sqrt(m);
        curE = E[r1*Ndp+c1];   /* E is a row-major array */
        integralVal = 0.0;

        for ( j=c1; j<=c2; j++ ){
          idx1 = _round(j*s);
          idx2 = _round((m*(j-c1)+r1)*s);

          d = 0.0;
          for ( k=0; k<dim; k++ ){
            /* q1 and q2 are column-major arrays */
            dxi = q1[idx1*dim+k] - q2[idx2*dim+k]*rootm;
            d += dxi*dxi;
          }

          integralVal += d * ds;
        }

        curE += integralVal;
        if ( curE < E[r2*Ndp+c2] ){
          /* E and P are row-major arrays */
          E[r2*Ndp+c2] = curE;
          P[r2*Ndp+c2] = r1*Ndp+c1;
        }
      }
    }
  }

  /* mexPrintf( "DP score: %f\n", E[Ndp*Ndp-1] ); */

  /* Reconstruct best path from (0,0) to (1,1) */
  gamma[0] = 0.;  gamma[Ndp-1] = 1.;
  c2=Ndp-1;  r2=Ndp-1;
  while( c2>0 ){
    c1=P[r2*Ndp+c2]%Ndp;
    r1=P[r2*Ndp+c2]/Ndp;
    m = ((double)(r2-r1))/((double)(c2-c1));
    for ( i=c1; i<c2; i++ ){
      gamma[i] = (m*(i-c1)+r1) * ds;
    }
    c2=c1;  r2=r1;
  }


cleanup:
  if ( E ) mxFree( E );
  if ( P ) mxFree( P );
}
