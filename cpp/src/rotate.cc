/*
 * LibSRVF - a shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012  Daniel Robinson
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "matrix.h"
#include "srvf.h"
#include "rotate.h"
#include "exceptions.h"

namespace srvf
{


// Computes the sign of the determinant of A
static int sgndet_(gsl_matrix *A)
{
  int signum;
  int status=0;
  int res;

  gsl_matrix *LU = gsl_matrix_alloc(A->size1, A->size2);
  gsl_permutation *perm = gsl_permutation_alloc(A->size1);
  if (!LU || !perm)
  {
    status=-1;
    goto cleanup;
  }

  gsl_matrix_memcpy(LU, A);
  gsl_linalg_LU_decomp(LU, perm, &signum);
  res = gsl_linalg_LU_sgndet(LU, signum);

cleanup:
  if (LU)   gsl_matrix_free(LU);
  if (perm) gsl_permutation_free(perm);
  if (status==-1)
  {
    throw std::bad_alloc();
  }

  return res;
}


// Returns a dim x dim matrix representing rotation by 
// theta radians about 
static Matrix perturbation_matrix_(double theta, size_t dim)
{
  double sinx, cosx;

  sinx = sin(theta);
  cosx = cos(theta);

  Matrix R = Matrix::identity(dim);
  R(0,0) = cosx; R(0,1) = -sinx;
  R(1,0) = sinx; R(1,1) = cosx;
  return R;
}


// Computes integral of Q1 * Q2^T
static void build_A_(const Srvf &Q1, const Srvf &Q2, gsl_matrix *A)
{
  std::vector<double> tv=srvf::util::unique(Q1.params(), Q2.params());
  size_t dim = Q1.dim();
  size_t npts = tv.size();

  Pointset Q1vals(dim, npts);
  Pointset Q2vals(dim, npts);

  Q1.evaluate(tv, Q1vals);
  Q2.evaluate(tv, Q2vals);

  gsl_matrix_set_zero(A);
  for (size_t i=0; i<npts-1; ++i)
  {
    double dt = tv[i+1] - tv[i];

    for (size_t j=0; j<dim; ++j)
    {
      for (size_t k=0; k<dim; ++k)
      {
        double vcur = gsl_matrix_get(A, j, k);
        vcur += Q1vals[i][j] * Q2vals[i][k] * dt;
        gsl_matrix_set(A, j, k, vcur);
      }
    }
  }
}


/**
 * Computes the optimal rotational alignment of \a Q2 to \a Q1.
 *
 * Returns a \c Matrix representing a rotation \f$ R \f$ which minimizes 
 * the L^2 norm of \f$ Q_1 - R Q_2 \f$.
 *
 * \a Q1 and \a Q2 must have the same dimension.  If the dimension 
 * is 1, optimizing over rotations doesn't make sense, so the result 
 * will be the 1x1 identity matrix.
 */
Matrix optimal_rotation (const Srvf &Q1, const Srvf &Q2)
{
  int status=0;  // success
  int sgndet;
  size_t dim = Q1.dim();
  Matrix R;
  Matrix Q2pert;

  if (Q1.dim() != Q2.dim())
    throw std::invalid_argument("Q1 and Q2 must have same dimension.");

  if (Q1.dim() < 2)
  {
    return Matrix(1,1,1.0);
  }

  gsl_matrix *A    = gsl_matrix_alloc(dim, dim);
  gsl_matrix *V    = gsl_matrix_alloc(dim, dim);
  gsl_vector *S    = gsl_vector_alloc(dim);
  gsl_matrix *res  = gsl_matrix_alloc(dim, dim);
  gsl_vector *work = gsl_vector_alloc(dim);

  if ( !A || !V || !S || !res || !work )
  {
    status=-1; // GSL alloc error, so throw after cleanup
    goto cleanup;
  }

  // Compute A = \int_0^1 q1(t) q2(t)^T dt
  build_A_(Q1, Q2, A);

  // Compute sign of det(A)
  // If A is singular, apply a small rotation to Q2 and try again.
  sgndet = sgndet_(A);
  while(sgndet == 0){
    Matrix Rpert = perturbation_matrix_(0.001, dim);
    Srvf Q2pert(Q2);
    Q2pert.rotate(Rpert);
    build_A_(Q1, Q2pert, A);
    sgndet = sgndet_(A);
  }

  // Compute SVD of A
  if (gsl_linalg_SV_decomp(A, V, S, work) != GSL_SUCCESS)
  { 
    status=-3;  // SVD failed, so throw after cleanup
    goto cleanup;
  }
  
  // If det(A) < 0, change the sign of the last column of V
  if (sgndet < 0){
    for (size_t i=0; i<dim; i++)
    {
      double x = gsl_matrix_get(V, i, dim-1);
      gsl_matrix_set(V, i, dim-1, -x);
    }
  }

  // Rotation matrix is AV^T
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, V, 0.0, res);

  // Convert result to an srvf::Matrix.
  // gsl_matrix and srvf::Matrix both use row-major order.
  R=Matrix(dim, dim, res->data, Matrix::ROW_MAJOR);

cleanup:
  if ( A )    gsl_matrix_free(A);
  if ( V )    gsl_matrix_free(V);
  if ( S )    gsl_vector_free(S);
  if ( work ) gsl_vector_free(work);
  if ( res )  gsl_matrix_free(res);

  // status!=0 indicates an error
  if (status==-1)
    throw std::bad_alloc();
  else if (status==-2)
    throw AlgorithmFailure("A is singular.");
  else if (status==-3)
    throw AlgorithmFailure("SVD failed.");

  return R;
}

} // namespace srvf
