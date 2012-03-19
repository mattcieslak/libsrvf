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

#include <cmath>
#include <stdexcept>

#include "plf.h"
#include "matrix.h"
#include "util.h"
#include "interp.h"

namespace srvf
{

/**
 * Creates a new \c Plf with the given sample points.
 *
 * Corresponding parameter values will be uniformly spaced from 0 to 1.
 *
 * \param samples
 */
Plf::Plf(const Matrix &samples)
  : samps_(samples), params_()
{
  params_=srvf::util::linspace(0.0,1.0,samples.cols());
}

/**
 * Creates a new \c Plf with the given sample points and parameters.
 *
 * \param samples
 * \param parameters
 */
Plf::Plf(const Matrix &samples, const Matrix &parameters)
  : samps_(samples), params_(parameters)
{
  if (params_.rows()!=1)
    std::invalid_argument("parameters must have 1 row");
  if (samps_.cols()!=params_.cols())
    std::invalid_argument("samples.cols()!=parameters.cols()");
}

/**
 * Evaluate the PLF at the given parameter value.
 *
 * The result is stored in the first column of \a result.
 *
 * \param t the parameter value
 * \param result a \c Matrix to receive the result.  Must have at least 
 *   as many rows as the dimension of this function.
 */
void Plf::evaluate(double t, Matrix &result)
{
  Matrix tv(1,1,&t);
  evaluate(tv,result);
}

/**
 * Evaluate the PLF at several parameter values.
 *
 * The result will be stored in \a result, with one point per column.
 * The caller is responsible for making sure that \a result is large enough 
 * to hold the result.
 *
 * \param tv the parameter values at which the function will be evaluated
 * \param result a \c Matrix to receive the result
 */
void Plf::evaluate(const Matrix &tv, Matrix &result)
{
  srvf::interp::interp_linear(samps_,params_,tv,result);
}

/**
 * Computes the preimages of the numbers in \a tv.
 *
 * This function is only supported if the \c Plf is an increasing 
 * 1-D function (i.e. \c samples().rows()==1 and samples(i)<=samples(i+1) for 
 * \c i=0,...,samples.cols()-2 ).
 *
 * The numbers in \a tv must be sorted in non-decreasing order, and 
 * must lie between the first and last elements of \c samples(), provided 
 * that this \c Plf is non-empty.  If this \c Plf is empty (i.e. 
 * samples().size()==0), then this routine will return immediately and 
 * \a result will be left unchanged.
 *
 * The matrix \a result must be the same size as \a tv.
 *
 * \param tv
 * \param result
 */
void Plf::preimages(const Matrix &tv, Matrix &result)
{
  if (samps_.rows()>1)
    std::logic_error("preimages() only supported for 1-D functions");
  if (tv.rows()!=1) 
    std::invalid_argument("tv must have 1 row");
  if (result.rows()!=1) 
    std::invalid_argument("result must have 1 row");
  if (tv.cols()!=result.cols()) 
    std::invalid_argument("tv and result must have the same size");
  
  // Do nothing if this PLF is the empty map, or if tv is empty
  if (samps_.size()==0) return;
  if (tv.size()==0) return ;

  srvf::interp::interp_linear(params_, samps_, tv, result);
}

/**
 * Computes the arc length of the \c Plf.
 *
 * Since the function is piecewise-linear, the arc length is just the sum 
 * of the lengths of the linear segments.
 */
double Plf::arc_length() const
{
  double res=0.0;
  for (int i=0; i<samps_.cols()-1; ++i)
  {
    double cur=0.0;
    for (int j=0; j<samps_.rows(); ++j)
    {
      double dxi=samps_(j,i+1)-samps_(j,i);
      cur += dxi*dxi;
    }
    res += sqrt(cur);
  }
  return res;
}

/**
 * Apply a translation to this \c Plf.
 *
 * Adds \a v to each of this function's sample points.  \a v must have 
 * the same number of rows as \c this->samples().
 *
 * \param v the vector by which to translate
 */
void Plf::translate(const Matrix &v)
{
  if (v.rows()!=samples().rows())
    std::invalid_argument("v has incorrect length");
  if (v.cols()<1)
    std::invalid_argument("v is empty");

  for (int i=0; i<samples().cols(); ++i)
  {
    for (int j=0; j<samples().rows(); ++j)
    {
      samples()(j,i)+=v(j,0);
    }
  }
}

/**
 * Apply a rotation to the curve.
 *
 * \param R a square rotation matrix with dimension \c samples().rows()
 */
void Plf::rotate(const Matrix &R)
{
  samps_=product(R,samps_);
}

/**
 * Apply a uniform scaling to the curve.
 *
 * \param s the scale factor.
 */
void Plf::scale(double s)
{
  samps_*=s;
}

void Plf::linear_combine(const Plf &F, double w1, double w2)
{

}

void Plf::compose(const Plf &F)
{

}

void Plf::invert()
{

}


Srvf Plf::to_srvf() const
{

}


Plf translation(const Plf &F, const Matrix &v)
{

}

Plf rotation(const Plf &F, const Matrix &R)
{

}

Plf scaling(const Plf &F, double s)
{

}

Plf linear_combination(const Plf &F1, const Plf &F2, 
             double w1, double w2)
{

}

Plf composition(const Plf &F1, const Plf &F2)
{

}

Plf inverse(const Plf &F)
{

}


}
