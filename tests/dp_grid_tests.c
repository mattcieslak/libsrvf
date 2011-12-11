#include <stdlib.h>
#include <time.h>
#include <CUnit/CUnit.h>

#include "dp_grid.h"


#if RAND_MAX < 100
#error "RAND_MAX < 100!!!"
#endif

static void random_partition( double *T, int n )
{
  double s;
  int i;

  T[0] = 0.0;
  for ( i=1; i<n; ++i )
  {
    T[i] = T[i-1] + (double)(rand() % 100);
  }

  s = T[n-1];
  for ( i=1; i<n; ++i )
  {
    T[i] /= s;
  }
}


static void test_param_to_idx_basic1()
{
  double T[] = { 0.0, 0.1, 0.4, 0.6, 1.0 };
  double params[] = { 
    0.0, 0.05, 0.09, 0.0999, 0.1, 0.11, 0.35, 0.3999, 0.4, 
    0.5, 0.6, 0.7, 0.75, 0.89, 0.9999, 1.0 };
  int indexes[] = {
    0, 0, 0, 0, 1, 1, 1, 1, 2, 
    2, 3, 3, 3, 3, 3, 3 };
  int nsamps, nparams;
  int i, idx;

  nsamps = sizeof(T) / sizeof(double);
  nparams = sizeof(params) / sizeof(double);

  for ( i=0; i<nparams; ++i )
  {
    idx = param_to_idx( T, nsamps, params[i] );
    CU_ASSERT_EQUAL( idx, indexes[i] );

    /*
    if ( idx < nsamps-2 )
      printf( "%f in [%f,%f)\n", params[i], T[idx], T[idx+1] );
    else
      printf( "%f in [%f,%f]\n", params[i], T[idx], T[idx+1] );
    */
  }
}


static void test_param_to_idx_random1()
{
  double T[10];
  int nsamps;
  int ntrials = 1000;
  double t;
  int i, idx;

  srand( time(NULL) );

  nsamps = sizeof(T) / sizeof(double);
  random_partition( T, nsamps );

  for ( i=0; i<ntrials; ++i )
  {
    t = (rand() % 100) / 99.0;
    idx = param_to_idx( T, nsamps, t );
    CU_ASSERT( (t>0.9999 && idx==nsamps-2) || (t >= T[idx] && t < T[idx+1]) );

    /*
    if ( idx < nsamps-2 )
      printf( "%f in [%f,%f)\n", t, T[idx], T[idx+1] );
    else
      printf( "%f in [%f,%f]\n", t, T[idx], T[idx+1] );
    */
  }
}

static void test_param_to_idx_random2()
{
  double T[256];
  int nsamps;
  int ntrials = 10000;
  double t;
  int i, idx;

  srand( time(NULL) );

  nsamps = sizeof(T) / sizeof(double);
  random_partition( T, nsamps );

  for ( i=0; i<ntrials; ++i )
  {
    t = (rand() % 100) / 99.0;
    idx = param_to_idx( T, nsamps, t );
    CU_ASSERT( (t>0.9999 && idx==nsamps-2) || (t >= T[idx] && t < T[idx+1]) );

    /*
    if ( idx < nsamps-2 )
      printf( "%f in [%f,%f)\n", t, T[idx], T[idx+1] );
    else
      printf( "%f in [%f,%f]\n", t, T[idx], T[idx+1] );
    */
  }
}

static void test_dp_edge_weight_basic1()
{
  double T1[] = { 0.0, 1.0/6.0, 0.5, 4.0/6.0, 5.0/6.0, 1.0 };
  double Q1[] = { 1.0, -1.0, 1.0, -1.0, 1.0 };
  double T2[] = { 0.0, 2.0/6.0, 5.0/6.0, 1.0 };
  double Q2[] = { -1.0, 1.0, -1.0 };
  double grid[625];
  double E[25];
  int P[25];
  double G[5];
  double T[5];
  int npts;

  double a[] = { 0.25, 0.75, 0.0, 0.0, 0.25 };
  double b[] = { 0.75, 1.0, 0.75, 0.25, 0.5 };
  double c[] = { 0.25, 0.25, 0.0, 0.0, 0.0 };
  double d[] = { 0.5, 1.0, 0.25, 0.25, 0.25 };
  double expected[] = { 0.51430, 0.90377, 0.90377, 0.66667, 0.0 };

  int nsamps1 = sizeof(T1) / sizeof(double);
  int nsamps2 = sizeof(T2) / sizeof(double);
  int dim = 1;

  double actual;
  int ncases = sizeof(expected) / sizeof(double);
  int i;

  for ( i=0; i<ncases; ++i )
  {
    actual = dp_edge_weight ( 
      Q1, T1, nsamps1, 
      Q2, T2, nsamps2, 
      dim, a[i], b[i], c[i], d[i] );
    CU_ASSERT_DOUBLE_EQUAL( actual, expected[i], 0.0001 );
  }

  dp_all_edge_weights( Q1, T1, 6, Q2, T2, 4, 1, grid, 5, 5 );
  dp_costs( Q1, T1, 6, Q2, T2, 4, 1, E, P, 5, 5 );
  npts = dp_build_gamma( P, 5, 5, G, T );

  printf( "gamma:\n" );
  for ( i=0; i<npts; ++i )
  {
    printf( "(%0.2f,%0.2f)\n", T[i], G[i] );
  }
}

static void test_dp_edge_weight_timing1()
{
  double T1[10];
  double Q1[] = {0.23, -0.2, 0.57, -0.12, 0.332, -0.14, 0.1, -0.3, 0.2, -0.1};
  double T2[8];
  double Q2[] = { -0.1, 0.23, -0.5, 0.6, -0.13, 0.2, -0.012, 0.27 };
  double E[80*80];
  int P[80*80];
  int grid_rows = 80, grid_cols = 80;

  int nsamps1 = sizeof(T1) / sizeof(double);
  int nsamps2 = sizeof(T2) / sizeof(double);

  random_partition( T1, nsamps1 );
  random_partition( T2, nsamps2 );

  dp_costs( Q1, T1, nsamps1, Q2, T2, nsamps2, 1, E, P, grid_rows, grid_cols );
}


CU_ErrorCode grid_tests_suite()
{
  CU_pSuite suite = CU_add_suite( "grid_tests", NULL, NULL );
  if ( suite == NULL )
  {
    return CU_get_error();
  }
  
  if ( CU_add_test( suite, "param_to_idx() basic test 1", 
                    (CU_TestFunc)test_param_to_idx_basic1 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "param_to_idx() random test 1", 
                    (CU_TestFunc)test_param_to_idx_random1 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "param_to_idx() random test 2", 
                    (CU_TestFunc)test_param_to_idx_random2 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_edge_weight() basic test 1", 
                    (CU_TestFunc)test_dp_edge_weight_basic1 ) == NULL )
  {
    return CU_get_error();
  }

  if ( CU_add_test( suite, "dp_edge_weight() timing test 1", 
                    (CU_TestFunc)test_dp_edge_weight_timing1 ) == NULL )
  {
    return CU_get_error();
  }

  return CUE_SUCCESS;
}
