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
#ifndef SRVF_POINTSET_H
#define SRVF_POINTSET_H 1

#include <cmath>
#include <vector>
#include <stdexcept>

#include "matrix.h"

namespace srvf
{

/**
 * Represents an ordered set of points in \f$ R^n \f$.
 */
class Pointset
{
public:
  
  /**
   * Indicates how a matrix containing sample points should be interpreted.
   */
  enum PackingMethod {
    POINT_PER_ROW=0,
    POINT_PER_COLUMN=1
  };

  /**
   * Creates an empty pointset.
   */
  Pointset () { }

  /**
   * Creates a new \c Pointset with the given size.
   *
   * The points will be uninitialized
   *
   * \param dim the dimension of the ambient space
   * \param npts the number of points
   */
  Pointset (size_t dim, size_t npts)
   : data_(npts,dim)
  { }

  /**
   * Creates a new \c Pointset with the given size.
   *
   * The points will be initialized so that all of their components are 
   * set equal to \a v.
   *
   * \param dim the dimension of the ambient space
   * \param npts the number of points
   * \param v all components of all points will be set to this value
   */
  Pointset (size_t dim, size_t npts, double v)
   : data_(npts,dim,v)
  { }

  /**
   * Creates a new \c Pointset representing the given points.
   *
   * \param points a \c Matrix containing the points
   * \param packing indicates how the points are stored in \a data
   */
  Pointset (const Matrix &points, PackingMethod packing=POINT_PER_ROW)
  { 
    // Internally, Pointsets always use point per row
    if (packing==POINT_PER_ROW)
    {
      data_=points;
    }
    else
    {
      data_=transpose(points);
    }
  }

  /** 
   * Creates a new \c Pointset initialized from the given \c vector.
   *
   * \param data
   * \param dim the dimension of the points.  Must divide \c data.size().
   *        Default is \c dim=1.
   * \param packing POINT_PER_ROW if each point is stored contiguously 
   *        in \a data, POINT_PER_COLUMN if each component is stored 
   *        contiguously.  Default is \c packing=POINT_PER_ROW.
   */
  Pointset (size_t dim, size_t npts, 
            const std::vector<double> &data, 
            PackingMethod packing=POINT_PER_ROW)
  {
    if (dim==0 || npts==0)
    {
      return;
    }
    if (npts != data.size()/dim)
    {
      throw std::invalid_argument("data.size() != dim*npts");
    }

    Matrix::Majorness majorness = 
      (packing==POINT_PER_ROW ? Matrix::ROW_MAJOR : Matrix::COLUMN_MAJOR);

    data_ = Matrix(npts,dim,data,majorness);
  }

  /** 
   * Creates a new \c Pointset initialized from the given data.
   *
   * \param data
   * \param dim the dimension of the points
   * \param npts the number of points
   * \param packing POINT_PER_ROW if each point is stored contiguously 
   *        in \a data, POINT_PER_COLUMN if each component is stored 
   *        contiguously.  Default is \c packing=POINT_PER_ROW.
   */
  Pointset (size_t dim, size_t npts, 
            const double *data, 
            PackingMethod packing=POINT_PER_ROW)
  {
    if (dim==0 || npts==0)
      return;
    if (!data)
      throw std::invalid_argument("data is NULL");

    Matrix::Majorness majorness = 
      (packing==POINT_PER_ROW ? Matrix::ROW_MAJOR : Matrix::COLUMN_MAJOR);

    data_ = Matrix(npts,dim,data,majorness);
  }

  /** Return the dimension of the space containing the points. */
  size_t dim() const { return data_.cols(); }

  /** Returns the number of points in the set. */
  size_t npts() const { return data_.rows(); }

  /** Returns component \a comp_idx of point \a point_idx. */
  double& operator() (size_t point_idx, size_t comp_idx)
  { return data_(point_idx,comp_idx); }

  /** Returns component \a comp_idx of point \a point_idx. */
  const double& operator() (size_t point_idx, size_t comp_idx) const
  { return data_(point_idx,comp_idx); }

  /** 
   * Computes the Euclidean distance between point \a i1 and point \a i2. 
   * 
   * \param i1 index of the first point
   * \param i2 index of the second point
   * \return the Euclidean distance between the two points
   */
  double distance(size_t i1, size_t i2) const
  {
    double d=0.0;
    for (size_t i=0; i<dim(); ++i)
    {
      double dxi=data_(i1,i)-data_(i2,i);
      d += dxi*dxi;
    }
    return sqrt(d);
  }

  /** 
   * Computes the dot product of point \a i1 and point \a i2. 
   * 
   * \param i1 index of the first point
   * \param i2 index of the second point
   * \return the dot product of the two points
   */
  double dot_product(size_t i1, size_t i2) const
  {
    double ip=0.0;
    for (size_t i=0; i<dim(); ++i)
    {
      ip += data_(i1,i)*data_(i2,i);
    }
    return ip;
  }

  /**
   * Computes the Euclidean norm of point \a idx.
   *
   * \param idx an index between 0 and \c npts()-1, inclusive
   * \result the norm of the point with index \a idx
   */
  double norm(size_t idx) const
  {
    double nrm=0.0;
    for (size_t i=0; i<dim(); ++i)
    {
      nrm += data_(idx,i)*data_(idx,i);
    }
    return sqrt(nrm);
  }

  /**
   * Computes the centroid of this \c Pointset.
   */
  std::vector<double> centroid() const
  {
    std::vector<double> res(dim(),0.0);
    if (npts() < 1) return res;

    for (size_t i=0; i<npts(); ++i)
    {
      for (size_t j=0; j<dim(); ++j)
      {
        res[j] += data_(i,j);
      }
    }

    for (size_t j=0; j<dim(); ++j)
    {
      res[j] /= (double)npts();
    }
    return res;
  }

  /**
   * Returns a \c vector containing this \c Pointset's data.
   *
   * \param packing if \c POINT_PER_ROW, then each point will be stored 
   *        contiguously in the result; otherwise, each component will 
   *        be stored contiguously.
   * \return a new \c vector containing the data.
   */
  std::vector<double> to_vector(PackingMethod packing=POINT_PER_ROW) const
  {
    std::vector<double> res(dim()*npts());
    for (size_t i=0; i<npts(); ++i)
    {
      for (size_t j=0; j<dim(); ++j)
      {
        if (packing==POINT_PER_ROW)
          res[i*dim()+j]=data_(i,j);
        else
          res[j*npts()+i]=data_(i,j);
      }
    }
    return res;
  }

  /**
   * Apply a translation to a single point in this \c Pointset.
   *
   * \param idx the index of the point to be translated
   * \param v a \c vector with size equal to the dimension of this \c Pointset
   */
  void translate (size_t idx, const std::vector<double> &v)
  {
    for (size_t i=0; i<data_.cols(); ++i)
    {
      data_(idx,i) += v[i];
    }
  }

  /**
   * Apply a translation to each point in this \c Pointset.
   *
   * \param v a \c vector with size equal to the dimension of this \c Pointset
   */
  void translate (const std::vector<double> &v)
  {
    for (size_t i=0; i<npts(); ++i)
    {
      translate (i, v);
    }
  }

  /**
   * Apply a rotation to a single point in this \c Pointset.
   *
   * \param idx the index of the point to be rotated
   * \param R a rotation matrix
   */
  void rotate (int idx, const Matrix &R)
  {
    // TODO
    throw std::logic_error("not implemented");
  }

  /**
   * Apply a rotation to each point in this \c Pointset.
   *
   * \param R a rotation matrix
   */
  void rotate (const Matrix &R)
  {
    // Each point is a row in the matrix, so multiply on right by 
    // the inverse (transpose) of R.
    // TODO: do the multiplication in place like rather than explicitly 
    // creating Rt
    Matrix Rt=transpose(R);
    data_ = product(data_,Rt);
  }

  /**
   * Apply a scaling to a single point in this \c Pointset.
   *
   * \param idx the index of the point to be scaled
   * \parm sf the scale factor
   */
  void scale (size_t idx, double sf)
  {
    for (size_t i=0; i<data_.cols(); ++i)
    {
      data_(idx,i) *= sf;
    }
  }

  /**
   * Apply a scaling to each point in this \c Pointset.
   *
   * \parm sf the scale factor
   */
  void scale (double sf)
  {
    data_ *= sf;
  }

private:
  
  Matrix data_;
};

/**
 * Returns the Euclidean distance between a pair of points in two 
 * \c Pointset instances.
 *
 * \param S1 the \c Pointset containing the first point
 * \param S2 the \c Pointset containing the second point
 * \param i1 the index of the first point
 * \param i2 the index of the first point
 * \return the distance between point \a i1 in \a S1 and point 
 *         \a i2 in \a S2.
 */
inline double distance(const Pointset &S1, const Pointset &S2, int i1, int i2)
{
  size_t dim=S1.dim();
  double d=0.0;

  for (size_t i=0; i<dim; ++i)
  {
    double dxi=S1(i1,i)-S1(i2,i);
    d += dxi*dxi;
  }
  return sqrt(d);
}

/**
 * Returns the dot product of a pair of points in two \c Pointset instances.
 *
 * \param S1 the \c Pointset containing the first point
 * \param S2 the \c Pointset containing the second point
 * \param i1 the index of the first point
 * \param i2 the index of the first point
 * \return the distance between point \a i1 in \a S1 and point 
 *         \a i2 in \a S2.
 */
inline double dot_product (
  const Pointset &S1, const Pointset &S2, 
  size_t i1, size_t i2)
{
  size_t dim=S1.dim();
  double ip=0.0;

  for (size_t i=0; i<dim; ++i)
  {
    ip += S1(i1,i)*S2(i2,i);
  }
  return ip;
}

/**
 * Computes a weighted sum of a point in \a S1 with a point in \a S2, and 
 * stores the result in \a S3.
 *
 * \param S1 \c Pointset containing first point
 * \param S2 \c Pointset containing second point
 * \param i1 index of first point
 * \param i2 index of second point
 * \param w1 weight for first point
 * \param w2 weight for second point
 * \param S3 \c Pointset to receive the result
 * \param i3 index of point in \a S3 that will receive the result
 */
inline void weighted_sum(
  const Pointset &S1, const Pointset &S2, 
  size_t i1, size_t i2, 
  double w1, double w2,
  Pointset &S3, size_t i3)
{
  if (S1.dim()!=S2.dim() || S2.dim()!=S3.dim())
    throw std::invalid_argument("Dimension mismatch");
  
  for (size_t i=0; i<S1.dim(); ++i)
  {
    S3(i3,i) = w1*S1(i1,i) + w2*S2(i2,i);
  }
}

} // namespace srvf

#endif // SRVF_POINTSET_H
