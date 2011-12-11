#ifndef DP_GRID_H
#define DP_GRID_H 1

/**
 * Computes weights of all edges in the DP matching graph, 
 * which is needed for the Floyd-Warshall all-pairs shortest-path 
 * algorithm.
 */
void dp_all_edge_weights( 
  double *Q1, double *T1, int nsamps1, 
  double *Q2, double *T2, int nsamps2,
  int dim, double *grid, int grid_rows, int grid_cols );

/**
 * Computes cost of best path from (0,0) to all other gridpoints.
 */
void dp_costs(
  double *Q1, double *T1, int nsamps1, 
  double *Q2, double *T2, int nsamps2,
  int dim, double *E, int *P, int grid_rows, int grid_cols );

/**
 * Computes the weight of the edge from (a,c) to (b,d) int the DP grid.
 */
double dp_edge_weight(
  double *Q1, double *T1, int nsamps1, 
  double *Q2, double *T2, int nsamps2,
  int dim,
  double a, double b, 
  double c, double d );

/**
 * Given predecessor table P, builds the piecewise-linear reparametrization 
 * function gamma.
 */
int dp_build_gamma( int *P, int grid_rows, int grid_cols, 
  double *G, double *T );

/**
 * Given t in [0,1], return the integer i such that t lies in the interval 
 * [T[i],T[i+1]) (or returns n-2 if t==1).
 */
int param_to_idx( double *T, int n, double t );

#endif /* DP_GRID_H */
