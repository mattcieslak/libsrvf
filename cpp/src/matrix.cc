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
#include <stdexcept>

#include "matrix.h"

namespace srvf{

#define __ELEMENTWISE_OP(A,B,C,op) \
  int nc=(A).size(); \
  for (int i=0; i<nc; ++i) { \
    (A)(i) = (B)(i) op (C)(i); \
  }

#define __SCALAR_OP(A,B,v,op) \
  int nc=(A).size(); \
  for (int i=0; i<nc; ++i) { \
    (A)(i) = (B)(i) op (v); \
  }
  
#define __MATRIX_MUL(A,B,C) \
  for (int i=0; i<(B).rows(); ++i) { \
    for (int j=0; j<(C).cols(); ++j) { \
      (A)(i,j)=0.0; \
      for (int k=0; k<(B).cols(); ++k) { \
        (A)(i,j) += (B)(i,k)*(C)(k,j); \
  } } } \

/**
 * Create a new \c Matrix with the specified size.
 *
 * \param rows
 * \param cols
 */
Matrix::Matrix (int rows, int cols)
  : data_(0), rows_(rows), cols_(cols)
{
  if (rows<0) throw std::invalid_argument("rows < 0");
  if (cols<0) throw std::invalid_argument("cols < 0");
  if (rows*cols<0) throw std::overflow_error("rows*cols<0");

  data_ = new double[rows*cols];
}

/**
 * Create a new \c Matrix with the specified size, with all elements 
 * set to the given value.
 *
 * \param rows
 * \param cols
 * \param val
 */
Matrix::Matrix (int rows, int cols, double val)
  : data_(0), rows_(rows), cols_(cols)
{
  if (rows<0) throw std::invalid_argument("rows < 0");
  if (cols<0) throw std::invalid_argument("cols < 0");
  if (rows*cols<0) throw std::overflow_error("rows*cols<0");

  data_ = new double[rows*cols];
  for (int i=0; i<rows*cols; ++i)
    data_[i]=val;
}

/**
 * Create a new \c Matrix initialized with the given data.
 *
 * A deep copy of the data is made.
 *
 * \param rows
 * \param cols
 * \param data pointer to a \c (rows)x(cols) array of \c doubles
 */
Matrix::Matrix (int rows, int cols, double *data)
  : data_((double*)0), rows_(rows), cols_(cols)
{
  if (rows<0) throw std::invalid_argument("rows < 0");
  if (cols<0) throw std::invalid_argument("cols < 0");
  int nc=rows*cols;
  if (nc<0) throw std::overflow_error("rows*cols<0");
  if (!data) throw std::invalid_argument("data is null");
  
  data_=new double[nc];
  for (int i=0; i<nc; ++i)
  {
    data_[i]=data[i];
  }
}

/**
 * Copy constructor.
 *
 * Creates a deep copy of the given \c Matrix.
 *
 * \param A
 */
Matrix::Matrix (const Matrix &A)
  : data_(0), rows_(A.rows_), cols_(A.cols_)
{
  data_ = new double[rows_*cols_];
  for (int i=0; i<rows_*cols_; ++i)
    data_[i]=A.data_[i];
}

/**
 * Assignment operator.
 *
 * Sets the current \c Matrix to a deep copy of the given \c Matrix.
 *
 * \param A
 */
Matrix& Matrix::operator= (const Matrix &A)
{
  if (this != &A)
  {
    rows_ = A.rows_;
    cols_ = A.cols_;
    data_ = new double[rows_*cols_];
    for (int i=0; i<rows_*cols_; ++i)
      data_[i]=A.data_[i];
  }
  return *this;
}

/**
 * Sets this \c Matrix to the empty matrix, freeing its resources.
 */
void Matrix::clear()
{
  delete[] data_;
  rows_=0;
  cols_=0;
  data_=(double*)0;
}

/**
 * Resizes this \c Matrix to have the specified size.
 *
 * If the matrix is non-empty, then any entries in the new matrix that 
 * were present in the old matrix will still have the same value as before 
 * the call to \c resize().  New entries are not initialized.
 *
 * If either of \a rows or \a cols is zero, the Matrix is cleared.
 *
 * \param rows the new number of rows
 * \param cols the new number of columns
 */
void Matrix::resize(int rows, int cols)
{
  if (rows<0) throw std::invalid_argument("rows<0");
  if (cols<0) throw std::invalid_argument("cols<0");
  int new_size=rows*cols;
  if (new_size<0) throw std::overflow_error("rows*cols<0");
  if (new_size==0)
  { clear(); return; }

  if (rows!=rows_ || cols!=cols_)
  {
    double *old_data=data_;
    int old_cols=cols_;
    int copy_rows=(rows<rows_ ? rows : rows_);
    int copy_cols=(cols<cols_ ? cols : cols_);

    data_=new double[new_size];
    rows_=rows;
    cols_=cols;
    
    for (int i=0; i<copy_rows; ++i)
    {
      for (int j=0; j<copy_cols; ++j)
      {
        data_[i*cols_+j]=old_data[i*old_cols+j];
      }
    }
    delete[] old_data;
  }
  else
  {
    data_=new double[new_size];
    rows_=rows;
    cols_=cols;
  }
}

/**
 * Add \a A to this \c Matrix.
 *
 * \param A
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator+= (const Matrix &A)
{
  if (rows_!=A.rows_ || cols_!=A.cols_)
    std::invalid_argument("size mismatch");

  __ELEMENTWISE_OP(*this,*this,A,+);
  return *this;
}

/**
 * Subtract \a A from this \c Matrix.
 *
 * \param A
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator-= (const Matrix &A)
{
  if (rows_!=A.rows_ || cols_!=A.cols_)
    std::invalid_argument("size mismatch");

  __ELEMENTWISE_OP(*this,*this,A,-);
  return *this;
}

/**
 * Multiply this \c Matrix elementwise by \a A.
 *
 * \param A
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator*= (const Matrix &A)
{
  if (rows_!=A.rows_ || cols_!=A.cols_)
    std::invalid_argument("size mismatch");

  __ELEMENTWISE_OP(*this,*this,A,*);
  return *this;
}

/**
 * Divide this \c Matrix elementwise by \a A.
 *
 * \param A
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator/= (const Matrix &A)
{
  if (rows_!=A.rows_ || cols_!=A.cols_)
    std::invalid_argument("size mismatch");

  __ELEMENTWISE_OP(*this,*this,A,/);
  return *this;
}

/**
 * Add \a v to all elements of this \c Matrix.
 *
 * \param v
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator+= (double v)
{
  __SCALAR_OP(*this,*this,v,+);
}

/**
 * Subtract \a v from all elements of this \c Matrix.
 *
 * \param v
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator-= (double v)
{
  __SCALAR_OP(*this,*this,v,-);
}

/**
 * Multiply all elements of this \c Matrix by \a v.
 *
 * \param v
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator*= (double v)
{
  __SCALAR_OP(*this,*this,v,*);
}

/**
 * Divide all elements of this \c Matrix by \a v.
 *
 * \param v
 * \return a reference to this \c Matrix
 */
Matrix& Matrix::operator/= (double v)
{
  __SCALAR_OP(*this,*this,v,/);
}

/**
 * Elementwise sum of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing \c A+B
 */
Matrix operator+ (const Matrix &A, const Matrix &B)
{
  if (A.rows()!=B.rows() || A.cols()!=B.cols())
    std::invalid_argument("size mismatch");

  Matrix R(A.rows(),A.cols());
  __ELEMENTWISE_OP(R,A,B,+);
  return R;
}

/**
 * Elementwise difference of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing \c A-B
 */
Matrix operator- (const Matrix &A, const Matrix &B)
{
  if (A.rows()!=B.rows() || A.cols()!=B.cols())
    std::invalid_argument("size mismatch");

  Matrix R(A.rows(),A.cols());
  __ELEMENTWISE_OP(R,A,B,-);
  return R;
}

/**
 * Elementwise product of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing \c A.*B
 */
Matrix operator* (const Matrix &A, const Matrix &B)
{
  if (A.rows()!=B.rows() || A.cols()!=B.cols())
    std::invalid_argument("size mismatch");

  Matrix R(A.rows(),A.cols());
  __ELEMENTWISE_OP(R,A,B,*);
  return R;
}

/**
 * Elementwise division of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing \c A.*B
 */
Matrix operator/ (const Matrix &A, const Matrix &B)
{
  if (A.rows()!=B.rows() || A.cols()!=B.cols())
    std::invalid_argument("size mismatch");

  Matrix R(A.rows(),A.cols());
  __ELEMENTWISE_OP(R,A,B,/);
  return R;
}

/**
 * Scalar addition.
 *
 * \param A
 * \param v
 */
Matrix operator+ (const Matrix &A, double v)
{
  Matrix R(A);
  R+=v;
  return R;
}

/**
 * Scalar subtraction.
 *
 * \param A
 * \param v
 */
Matrix operator- (const Matrix &A, double v)
{
  return A+(-v);
}

/**
 * Scalar multiplication.
 *
 * \param A
 * \param v
 */
Matrix operator* (const Matrix &A, double v)
{
  Matrix R(A);
  R*=v;
  return R;
}

/**
 * Scalar division.
 *
 * \param A
 * \param v
 */
Matrix operator/ (const Matrix &A, double v)
{
  Matrix R(A);
  R/=v;
  return R;
}

/**
 * Matrix product of two matrices.
 *
 * \param A
 * \param B
 * \return a new \c Matrix representing the matrix product \c AB
 */
Matrix product(const Matrix &A, const Matrix &B)
{
  if (A.cols()!=B.rows()) std::invalid_argument("size mismatch");

  Matrix R(A.rows(),B.cols());
  __MATRIX_MUL(R,A,B);
  return R;
}

/**
 * Returns a new \c Matrix representing the transpose of the given \c Matrix.
 *
 * \param A
 * \return the transpose of A
 */
Matrix transpose(const Matrix &A)
{
  Matrix At(A.cols(),A.rows());
  for (int i=0; i<A.rows(); ++i)
  {
    for (int j=0; j<A.cols(); ++j)
    {
      At(j,i)=A(i,j);
    }
  }
  return At;
}

} // namespace srvf
