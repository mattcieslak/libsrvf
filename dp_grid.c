#include <stdio.h>
#include <math.h>
#include "dp_nbhd.h"
#include "dp_grid.h"


void dp_all_edge_weights(
  double *Q1, double *T1, int nsamps1, 
  double *Q2, double *T2, int nsamps2,
  int dim, 
  double *grid, int grid_rows, int grid_cols )
{
  int sr, sc;  /* source row and column */
  int tr, tc;  /* target row and column */
  int l1, l2, l3;  /* for multidimensional array mapping */
  double a, c; /* source coordinates */
  double b, d; /* target coordinates */
  int i;
  
  for ( i=0; i<grid_cols*grid_rows*grid_cols*grid_rows; grid[i++]=1e6 );

  /* grid is a grid_rows x grid_cols x grid_rows x grid_cols array.  
   * Weight of edge from (T1[i],T2[j]) to (T1[k],T2[l]) (Cartesian coordinates) 
   * is in grid(j,i,l,k).
   * Mapping: 
   *  (j,i,l,k) :--> j*grid_cols*grid_rows*grid_cols + 
   *                 i*grid_rows*grid_cols +
   *                 l*grid_cols + 
   *                 k
   */
  l1 = grid_cols * grid_rows * grid_cols;
  l2 = grid_rows * grid_cols;
  l3 = grid_cols;

  for ( tr=1; tr<grid_rows; ++tr )
  {
    d = ((double)tr) / ((double)(grid_rows-1));
    for ( tc=1; tc<grid_cols; ++tc )
    {
      b = ((double)tc) / ((double)(grid_cols-1));
      for ( i=0; i<DP_NBHD_COUNT; ++i )
      {
        sr = tr - dp_nbhd[i][0];
        sc = tc - dp_nbhd[i][1];

        if ( sr < 0 || sc < 0 ) continue;

        a = ((double)sc) / ((double)(grid_cols-1));
        c = ((double)sr) / ((double)(grid_rows-1));

        /* grid(sr,sc,tr,tc) */
        grid[sr*l1+sc*l2+tr*l3+tc] = 
         dp_edge_weight( Q1, T1, nsamps1, Q2, T2, nsamps2, dim, a, b, c, d );
        
        /*
        printf( "(%0.2f,%0.2f) --> (%0.2f,%0.2f) = %0.2f\n", 
          a, c, b, d, grid[sr*l1+sc*l2+tr*l3+tc] );
        */
      }
    }
  }
}


void dp_costs(
  double *Q1, double *T1, int nsamps1, 
  double *Q2, double *T2, int nsamps2,
  int dim, double *E, int *P, int grid_rows, int grid_cols )
{
  int sr, sc;  /* source row and column */
  int tr, tc;  /* target row and column */
  double a, c; /* source coordinates */
  double b, d; /* target coordinates */
  double w, cand_cost;
  int i;
  
  E[0] = 0.0;
  for ( i=1; i<grid_cols; E[i++]=1e9 );
  for ( i=1; i<grid_rows; E[grid_cols*i++]=1e9 );

  for ( tr=1; tr<grid_rows; ++tr )
  {
    d = ((double)tr) / ((double)(grid_rows-1));
    for ( tc=1; tc<grid_cols; ++tc )
    {
      b = ((double)tc) / ((double)(grid_cols-1));
      E[grid_cols*tr + tc] = 1e9;

      for ( i=0; i<DP_NBHD_COUNT; ++i )
      {
        sr = tr - dp_nbhd[i][0];
        sc = tc - dp_nbhd[i][1];

        if ( sr < 0 || sc < 0 ) continue;

        a = ((double)sc) / ((double)(grid_cols-1));
        c = ((double)sr) / ((double)(grid_rows-1));

        w = dp_edge_weight( Q1, T1, nsamps1, Q2, T2, nsamps2, dim, a, b, c, d );
        cand_cost = E[grid_cols*sr+sc] + w;
        if ( cand_cost < E[grid_cols*tr+tc] )
        {
          E[grid_cols*tr+tc] = cand_cost;
          P[grid_cols*tr+tc] = grid_cols*sr + sc;
        }
      }
    }
  }

  tr = grid_rows - 1;
  tc = grid_cols - 1;

  for ( tr=1; tr<grid_rows; ++tr )
  {
    for ( tc=1; tc<grid_cols; ++tc )
    {
      printf( "E[%d,%d]=%0.3f, ", tr, tc, E[grid_cols*tr+tc] );
      printf( "P[%d,%d]=(%d,%d)\n", tr, tc, P[grid_cols*tr+tc]/grid_cols,
                                            P[grid_cols*tr+tc]%grid_cols );
    }
  }
}

double dp_edge_weight(
  double *Q1, double *T1, int nsamps1, 
  double *Q2, double *T2, int nsamps2,
  int dim,
  double a, double b, 
  double c, double d )
{
  double res = 0.0;
  int Q1idx, Q2idx;
  int Q1idxnext, Q2idxnext;
  double t1, t2;
  double t1next, t2next;
  double t1nextcand1, t1nextcand2;
  double slope, rslope;
  double dq, dqi;
  int i;

  Q1idx = param_to_idx( T1, nsamps1, a );
  Q2idx = param_to_idx( T2, nsamps2, c );

  t1 = a;
  t2 = c;

  slope = (d-c) / (b-a);
  rslope = sqrt( slope );

  while( t1 < b && t2 < d )
  {
    /* Find endpoint of current interval */
    t1nextcand1 = T1[Q1idx+1];
    t1nextcand2 = a + (T2[Q2idx+1]-c) / slope; /* preimage of T2[Q2idx+1] */

    if ( t1nextcand1 < t1nextcand2 )
    {
      t1next = t1nextcand1;
      t2next = c + slope * (t1next - a);
      Q1idxnext = Q1idx+1;
      Q2idxnext = Q2idx;
    } else {
      t1next = t1nextcand2;
      t2next = T2[Q2idx+1];
      Q1idxnext = Q1idx;
      Q2idxnext = Q2idx+1;
    }

    if ( t1next > b ) t1next = b;
    if ( t2next > d ) t2next = d;

    /* Get contribution for current interval */
    dq = 0.0;
    for ( i=0; i<dim; ++i )
    {
      dqi = Q1[Q1idx] - rslope * Q2[Q2idx];
      dq += dqi*dqi;
    }
    res += (t1next - t1) * dq;

    t1 = t1next;
    t2 = t2next;
    Q1idx = Q1idxnext;
    Q2idx = Q2idxnext;
  }

  return res;
}


/**
 * Given predecessor matrix P (as created by dp_costs()), build the 
 * corrseponding piecewise-linear reparametrization function gamma.
 *
 * \param P P[i*grid_cols+j] = predecessor of gridpoint at row i, column j.  
 *   If predecessor is at row pr and column pc, then 
 *   P[i*grid_cols+j] = pr*grid_cols + pc.
 * \param grid_rows
 * \param grid_cols
 * \param G pointer to array of length at least \c grid_cols, to receive 
 *   the function values.
 * \param T pointer to array of length at least \c grid_cols, to receive 
 *   the corresponding parameter values.
 */
int dp_build_gamma( int *P, int grid_rows, int grid_cols, double *G, double *T )
{
  int sr, sc;
  int tr, tc;
  int p, i;
  double dx = 1.0 / (grid_cols - 1);
  double dy = 1.0 / (grid_rows - 1);
  int npts;  /* result = length of Tg */

  /* Dry run first, to determine length of Tg */
  npts = 1;
  tr = grid_rows-1;
  tc = grid_cols-1;
  while( tr > 0 && tc > 0 )
  {
    p = P[tr*grid_cols+tc];
    tr = p / grid_cols;
    tc = p % grid_cols;
    ++npts;
  }

  G[npts-1] = 1.0;
  T[npts-1] = 1.0;

  tr = grid_rows-1;
  tc = grid_cols-1;
  i = npts-2;
  while( tr > 0 && tc > 0 )
  {
    p = P[tr*grid_cols+tc];
    sr = p / grid_cols;
    sc = p % grid_cols;
    
    G[i] = dy * sr;
    T[i] = dx * sc;

    tr = sr;
    tc = sc;
    --i;
  }

  return npts;
}


int param_to_idx( double *T, int n, double t )
{
  int l, m, r;

  if ( t < 0.9999 )
  {
    l=0;
    r=n;
    m=(l+r)/2;

    /* TODO: are these comparisons OK??? */
    while( 1 )
    {
      if ( t >= T[m+1] )
        l = m;
      else if ( t < T[m] )
        r = m;
      else
        break;
      
      m = (r+l)/2;
    }

    return m;
  } else {
    return n-2;
  }
}

