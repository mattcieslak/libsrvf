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
#ifndef SRVF_INTERP_H
#define SRVF_INTERP_H 1

#include "matrix.h"

namespace srvf
{

namespace interp
{


/**
 * 1-D table lookup.
 *
 * Given a table of numbers \c table[0]<=table[1]<=...<=table[n-1] and 
 * a number t, returns the smallest integer \c idx such that 
 * \c table[idx]<=t<table[idx+1].  If \c t<table[0], then the result is 
 * \c -1, and if \c t>=table[n-1], then the result is \c n-1.
 *
 * \param table an array sorted in non-decreasing order
 * \param n the number of elements in \a table
 * \param t the query
 * \result an index between \c -1 and \c n-1, inclusive
 */
inline int lookup(const double *table, int n, double t)
{
  if (n>1)
  {
    if (t<table[0]) return -1;
    else if (t>=table[n-1]) return n-1;
    else
    {
      int i1=0, i3=n-2;
      int i2=(i1+i3)/2;
      while(i1<i3)
      {
        if (t<table[i2]) i3=i2;
        else if (t>=table[i2+1]) i1=i2+1;
        else break;

        i2=(i1+i3)/2;
      }
      return i2;
    }
  }
  else if (n==1)
  {
    return (t<table[0] ? -1 : 0);
  }
  else
    return -1;
}

/**
 * Computes a weighted sum of a column of \a A with a column of \a B.
 *
 * The result is stored in column \a idx3 of the matrix \a R, which 
 * must have the same number of rows as \a A and \a B.
 *
 * \param A the first matrix
 * \param B the second matrix
 * \param idx1 column index for \a A
 * \param idx2 column index for \a B
 * \param w1 first weight
 * \param w2 second weight
 * \param R a matrix to receive the result
 * \param idx3 column index for \a R
 */
inline void weighted_column_sum(const Matrix &A, const Matrix &B, 
                                int idx1, int idx2, 
                                double w1, double w2, 
                                Matrix &R, int idx3)
{
  for (int i=0; i<A.rows(); ++i)
  {
    R(i,idx3) = w1*A(i,idx1) + w2*B(i,idx2);
  }
}

/**
 * Linear interpolation.
 *
 * If \c params(i)==params(i+1) for some \c i, the interpolant has a jump 
 * discontinuity.  In this case, we take the function to be right-continuous 
 * at that point.
 *
 * \param samps a \c Matrix containing the sample points, one point per column
 * \param params the parameter values corresponding to \a samps
 * \param tv parameters at which to interpolated.  Must be non-decreasing.
 * \param result [output] a \c Matrix to hold the result
 */
inline void interp_linear(const Matrix &samps, const Matrix &params, 
                          const Matrix &tv, Matrix &result)
{
  int idx=0;
  if (params.cols()>1 && tv(0)>params(1))
  {
    idx=lookup(params.data(), params.cols(), tv(0));
    if (idx<0) idx=0;
  }

  for (int i=0; i<tv.cols(); ++i)
  {
    double tvi=tv(i);
    while (idx<params.cols()-1 && tvi>=params(idx+1)) 
    {
      ++idx;
    }
    if (tvi<params(idx))
    {
      tvi=params(idx);
    }

    if (idx<params.cols()-1)
    {
      double w1 = params(idx+1)-tvi;
      double w2 = tvi-params(idx);
      double w = w1+w2;

      if (w>1e-6)
      {
        w1 /= w;
        w2 /= w;
      }
      else
      {
        w1 = 0.0;
        w2 = 1.0;
      }

      weighted_column_sum(samps,samps,idx,idx+1,w1,w2,result,i);
    }
    else
    {
      weighted_column_sum(samps,samps,idx,idx,1.0,0.0,result,i);
    }
  }
}

} // namespace srvf::interp

} // namespace srvf

#endif // SRVF_INTERP_H
