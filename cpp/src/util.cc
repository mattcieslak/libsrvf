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

#include <algorithm>
#include <stdexcept>

#include "util.h"

namespace srvf
{

namespace util
{

Matrix linspace(double a, double b, int n)
{
  if (n<0) std::invalid_argument("n<0");

  Matrix R(1,n);
  if (n>1)
  {
    R(0)=a;
    double dt=(b-a)/(double)(n-1);
    for (int i=1; i<n; ++i)
    {
      R(i)=R(i-1)+dt;
    }
  }
  else if (n>0)
  {
    R(0)=a;
  }
  return R;
}

/**
 * Returns a new \c Matrix containing unique elements of \a v1 and \a v2, 
 * sorted in ascending order.
 *
 * The result will be a \c 1xN matrix, where \c N is the number of elements 
 * in the set-theoretic union of \c v1 and \c v2.
 *
 * For the purposes of this function, two \c doubles \c a and \c b are 
 * considered unique if \c fabs(a-b)>=1e-6 .
 *
 * \param v1
 * \param v2
 */
Matrix unique(Matrix v1, Matrix v2)
{
  int nv1=v1.size();
  int nv2=v2.size();
  double *buf = new double[nv1+nv2];
  int bufsize=0;

  std::sort(&(v1.data()[0]), &(v1.data()[nv1]));
  std::sort(&(v2.data()[0]), &(v2.data()[nv2]));
  
  int i1=0;
  int i2=0;
  while (i1<nv1 && i2<nv2)
  {
    if (v1(i1) < v2(i2)-1e-9)
    {
      buf[bufsize++]=v1(i1++);
    }
    else if (v2(i2) < v1(i1)-1e-9)
    {
      buf[bufsize++]=v2(i2++);
    }
    else
    {
      buf[bufsize++]=v1(i1++);
      ++i2;
    }
  }
  while (i1<nv1) buf[bufsize++]=v1(i1++);
  while (i2<nv2) buf[bufsize++]=v2(i2++);

  Matrix result(1,bufsize,buf);
  delete[] buf;
  return result;
}

/**
 * Compute the first-order differences of the columns of \a X.
 *
 * \param X
 * \return a new \c Matrix containing the first differences of the columns 
 *   of \a X.
 */
Matrix diff(const Matrix &X)
{
  Matrix res(X.rows(),X.cols()-1);
  for (int i=0; i<res.cols(); ++i)
  {
    for (int j=0; j<res.rows(); ++j)
    {
      res(j,i)=X(j,i+1)-X(j,i);
    }
  }
  return res;
}

/**
 * Forward difference approximation of first derivative.
 *
 * Equivalent to \c diff(X)./diff(tv)
 *
 * \param X
 * \param tv abscissae corresponding to the columns of \a X
 * \return a new \c Matrix containing the forward difference approximation
 */
Matrix diff(const Matrix &X, const Matrix &tv)
{
  if (X.cols()!=tv.cols())
    std::invalid_argument("X and tv must have the same length.");
  if (tv.rows()!=1)
    std::invalid_argument("tv must have 1 row");

  Matrix res(X.rows(),X.cols()-1);
  for (int i=0; i<res.cols(); ++i)
  {
    double dt=tv(i+1)-tv(i);
    for (int j=0; j<res.rows(); ++j)
    {
      res(j,i)=(X(j,i+1)-X(j,i)) / dt;
    }
  }
  return res;
}

} // namespace srvf::util

} // namespace srvf
