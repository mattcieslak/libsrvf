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
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include "point.h"
#include "matrix.h"

namespace srvf
{

/**
 * Represents an ordered set of points in \f$ R^n \f$.
 */
class Pointset
{
public:

  typedef Point value_type;
  typedef Point* pointer;
  typedef const Point* const_pointer;
  typedef Point& reference;
  typedef const Point& const_reference;
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  
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
  Pointset() { }

  /**
   * Creates a new \c Pointset with the given size.
   *
   * The points will be uninitialized
   *
   * \param dim the dimension of the ambient space
   * \param npts the number of points
   */
  Pointset (size_t dim, size_t npts)
   : data_(npts, Point(dim, 0.0))
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
   : data_(npts, Point(dim, v))
  { }

  /**
   * Creates a new \c Pointset representing the given points.
   *
   * \param points a \c Matrix containing the points
   * \param packing indicates how the points are stored in \a data
   */
  Pointset (const Matrix &points, PackingMethod packing=POINT_PER_ROW)
   : data_( (packing==POINT_PER_ROW ? points.rows() : points.cols()), 
       Point(packing==POINT_PER_ROW ? points.cols() : points.rows()) )
  { 
    if (packing==POINT_PER_ROW)
    {
      for (size_t i=0; i<points.rows(); ++i)
        for (size_t j=0; j<points.cols(); ++j)
          data_[i][j] = points(i, j);
    }
    else
    {
      for (size_t i=0; i<points.cols(); ++i)
        for (size_t j=0; j<points.rows(); ++j)
          data_[i][j] = points(j, i);
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
   : data_(npts, Point(dim))
  {
    if (dim==0 || npts==0) return;
    if (npts != data.size()/dim)
      throw std::invalid_argument("data.size() != dim*npts");
    
    if (packing==POINT_PER_ROW)
    {
      for (size_t i=0; i<npts; ++i)
        for (size_t j=0; j<dim; ++j)
          data_[i][j] = data[i*dim + j];
    }
    else
    {
      for (size_t i=0; i<npts; ++i)
        for (size_t j=0; j<dim; ++j)
          data_[i][j] = data[j*npts + i];
    }
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
   : data_(npts, Point(dim))
  {
    if (dim==0 || npts==0) return;
    if (!data) throw std::invalid_argument("data is NULL");

    if (packing==POINT_PER_ROW)
    {
      for (size_t i=0; i<npts; ++i)
        for (size_t j=0; j<dim; ++j)
          data_[i][j] = data[i*dim + j];
    }
    else
    {
      for (size_t i=0; i<npts; ++i)
        for (size_t j=0; j<dim; ++j)
          data_[i][j] = data[j*npts + i];
    }
  }

  /** Return the dimension of the space containing the points. */
  inline size_t dim() const { return (!data_.empty() ? data_[0].dim() : 0); }

  /** Returns the number of points in the set. */
  inline size_t npts() const { return data_.size(); }

  /** Indicates whether or not this Pointset contains any points. */
  inline bool empty() const { return (dim() > 0 && npts() > 0); }

  /** Returns a reference to the specified point. */
  inline Point &operator[] (size_t idx)
  {
    if (idx >= npts()) throw std::out_of_range("Point index out of range.");
    return data_[idx];
  }

  /** Returns a \c const reference to the specified point. */
  inline const Point &operator[] (size_t idx) const
  {
    if (idx >= npts()) throw std::out_of_range("Point index out of range.");
    return data_[idx];
  }

  /** Returns a reference to the first point in the set. */
  inline Point &front()
  {
    if (empty()) throw std::out_of_range("front() called on empty Pointset.");
    return data_[0];
  }

  /** Returns a \c const reference to the first point in the set. */
  inline const Point &front() const
  {
    if (empty()) throw std::out_of_range("front() called on empty Pointset.");
    return data_[0];
  }

  /** Returns a reference to the last point in the set. */
  inline Point &back()
  {
    if (empty()) throw std::out_of_range("back() called on empty Pointset.");
    return data_[npts()-1];
  }

  /** Returns a \c const reference to the last point in the set. */
  inline const Point &back() const
  {
    if (empty()) throw std::out_of_range("back() called on empty Pointset.");
    return data_[npts()-1];
  }

  /** Append a new point to the end of this pointset. */
  inline void push_back(const Point &P)
  {
    if (dim()>0 && P.dim() != dim())
      throw std::invalid_argument("New point has wrong dimension.");

    data_.push_back(P);
  }

  /** Remove the last point from this pointset. */
  inline void pop()
  {
    if (!data_.empty()) data_.pop_back();
  }

  /** 
   * Computes the Euclidean distance between point \a i1 and point \a i2. 
   * 
   * \param i1 index of the first point
   * \param i2 index of the second point
   * \return the Euclidean distance between the two points
   */
  inline double distance(size_t i1, size_t i2) const
  {
    if (i1 >= npts()) throw std::out_of_range("i1 out of range");
    if (i2 >= npts()) throw std::out_of_range("i2 out of range");

    if (i1 == i2) return 0.0;
    else return data_[i1].distance_to(data_[i2]);
  }

  /** 
   * Computes the dot product of point \a i1 and point \a i2. 
   * 
   * \param i1 index of the first point
   * \param i2 index of the second point
   * \return the dot product of the two points
   */
  inline double dot_product(size_t i1, size_t i2) const
  {
    if (i1 >= npts()) throw std::out_of_range("i1 out of range");
    if (i2 >= npts()) throw std::out_of_range("i2 out of range");

    return data_[i1].dot_product(data_[i2]);
  }

  /**
   * Computes the Euclidean norm of point \a idx.
   *
   * \param idx an index between 0 and \c npts()-1, inclusive
   * \result the norm of the point with index \a idx
   */
  inline double norm(size_t idx) const
  {
    if (idx >= npts()) throw std::out_of_range("idx out of range");
    return data_[idx].norm();
  }

  /**
   * Computes the centroid of this \c Pointset.
   */
  Point centroid() const
  {
    Point res(dim(), 0.0);
    if (npts() < 1) return res;

    for (size_t i=0; i<npts(); ++i)
      res += data_[i];

    res /= (double)npts();
    return res;
  }

  /**
   * Determines whether or not point \a i1 is equal to the product 
   * of point \a i2 and a non-negative scalar.
   *
   * The result is true if and only if one of these holds:
   *  - point \a i1 is zero (norm < \a tol), or
   *  - point \a i2 is zero (norm < \a tol), or
   *  - after scaling both points to unit norm, their l1 distance 
   *    is less than \a tol
   */
  inline bool on_same_ray(size_t i1, size_t i2, double tol=1e-4) const
  {
    if (i1 >= npts()) throw std::out_of_range("i1 out of range");
    if (i2 >= npts()) throw std::out_of_range("i2 out of range");

    return data_[i1].on_same_ray(data_[i2], tol);
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
          res[i*dim()+j]=data_[i][j];
        else
          res[j*npts()+i]=data_[i][j];
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
  void translate (size_t idx, const Point &v)
  {
    if (idx >= npts()) throw std::out_of_range("idx out of range");

    data_[idx] += v;
  }

  /**
   * Apply a translation to each point in this \c Pointset.
   *
   * \param v a \c vector with size equal to the dimension of this \c Pointset
   */
  void translate (const Point &v)
  {
    for (size_t i=0; i<npts(); ++i)
      translate(i, v);
  }

  /**
   * Apply a rotation to a single point in this \c Pointset.
   *
   * \param idx the index of the point to be rotated
   * \param R a rotation matrix
   */
  void rotate (int idx, const Matrix &R)
  {
    if (idx >= npts()) throw std::out_of_range("idx out of range");
    if (R.cols() != dim()) throw std::invalid_argument("size mismatch: R");
    
    data_[idx] = product(R, data_[idx]);
  }

  /**
   * Apply a rotation to each point in this \c Pointset.
   *
   * \param R a rotation matrix having the same dimension as the 
   *   points in this \c Pointset
   */
  void rotate (const Matrix &R)
  {
    if (R.cols() != dim()) throw std::invalid_argument("size mismatch: R");

    for (size_t i=0; i<npts(); ++i)
      data_[i] = product(R, data_[i]);
  }

  /**
   * Apply a scaling to a single point in this \c Pointset.
   *
   * \param idx the index of the point to be scaled
   * \parm sf the scale factor
   */
  void scale (size_t idx, double sf)
  {
    if (idx >= npts()) throw std::out_of_range("idx out of range");

    data_[idx] *= sf;
  }

  /**
   * Apply a scaling to each point in this \c Pointset.
   *
   * \parm sf the scale factor
   */
  void scale (double sf)
  {
    for (size_t i=0; i<npts(); ++i)
      data_[i] *= sf;
  }


  /**
   * An iterator that iterates over a single component of this \c Pointset.
   */
  class component_iterator
   : public std::iterator<std::random_access_iterator_tag, double>
  {
    friend class Pointset;

  public:

    bool operator== (const component_iterator &A) const
    { 
      return ( ((S_) == (A.S_)) && (point_idx_ == A.point_idx_) && 
               (comp_idx_ == A.comp_idx_) ); 
    }
    bool operator!= (const component_iterator &A) const
    { return !( (*this) == A ); }

    bool operator< (const component_iterator &A) const
    { return (((S_) == (A.S_)) && (point_idx_ < A.point_idx_)
      && (comp_idx_ == A.comp_idx_)); }
    bool operator<= (const component_iterator &A) const
    { return (((S_) == (A.S_)) && (point_idx_ <= A.point_idx_)
      && (comp_idx_ == A.comp_idx_)); }

    bool operator> (const component_iterator &A) const
    { return (((S_) == (A.S_)) && (point_idx_ > A.point_idx_)
      && (comp_idx_ == A.comp_idx_)); }
    bool operator>= (const component_iterator &A) const
    { return (((S_) == (A.S_)) && (point_idx_ >= A.point_idx_)
      && (comp_idx_ == A.comp_idx_)); }

    component_iterator &operator++ ()
    { point_idx_ = std::min( point_idx_ + 1, (ptrdiff_t)S_->npts() ); }
    component_iterator operator++ (int)
    { 
      component_iterator tmp(*this); 
      point_idx_ = std::min( point_idx_ + 1, (ptrdiff_t)S_->npts() ); 
      return tmp;
    }

    component_iterator &operator-- ()
    { point_idx_ = std::max(point_idx_-1, -1); return *this; }
    component_iterator operator-- (int)
    { component_iterator tmp(*this);  
      point_idx_ = std::max(point_idx_-1, -1);  return tmp; }

    component_iterator operator+ (ptrdiff_t n) const
    { return component_iterator(
        S_, comp_idx_, std::min(point_idx_+n, (ptrdiff_t)S_->npts())); }
    component_iterator &operator+= (ptrdiff_t n)
    { point_idx_ = std::min(point_idx_+n, (ptrdiff_t)S_->npts());  
      return *this; }

    component_iterator operator- (ptrdiff_t n) const
    { return component_iterator(
        S_, comp_idx_, std::max(point_idx_-n, -1)); }
    component_iterator &operator-= (ptrdiff_t n)
    { point_idx_ = std::max(point_idx_-n, -1); return *this; }

    ptrdiff_t operator- (const component_iterator &A)
    { return (point_idx_ - A.point_idx_); }

    double& operator* ()
    {
      // Bounds checking done in Pointset, exceptions propagate to caller
      return (*S_)[point_idx_][comp_idx_];
    }

    double& operator[] (ptrdiff_t n)
    {
      // Bounds checking done in Pointset, exceptions propagate to caller
      return (*S_)[point_idx_+n][comp_idx_];
    }

  private:
    
    component_iterator(Pointset *S, ptrdiff_t comp_idx, ptrdiff_t point_idx)
     : S_(S), comp_idx_(comp_idx), point_idx_(point_idx)
    { }

    Pointset *S_;
    ptrdiff_t comp_idx_;
    ptrdiff_t point_idx_;
  };


  /**
   * Returns an iterator to iterate over the points in this \c Pointset.
   */
  std::vector<Point>::iterator begin_points()
  {
    return data_.begin();
  }

  /**
   * Returns a past-the-end value for the points iterator.
   */
  std::vector<Point>::iterator end_points()
  {
    return data_.end();
  }

  /**
   * Returns an iterator for a single component of this \c Pointset.
   */
  component_iterator begin_component(size_t comp_idx)
  {
    return component_iterator(this, comp_idx, 0);
  }

  /**
   * Returns a past-the-end iterator for the specified component.
   */
  component_iterator end_component(size_t comp_idx)
  {
    return component_iterator(this, comp_idx, npts());
  }

private:
  
  std::vector<Point> data_;
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
  return S1[i1].distance_to(S2[i2]);
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
inline double 
dot_product ( const Pointset &S1, const Pointset &S2, size_t i1, size_t i2)
{
  return S1[i1].dot_product(S2[i2]);
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
  
  S3[i3] = (S1[i1]*w1) + (S2[i2]*w2);
}

} // namespace srvf

#endif // SRVF_POINTSET_H
