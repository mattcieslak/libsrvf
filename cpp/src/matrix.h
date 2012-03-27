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
#ifndef SRVF_MATRIX_H
#define SRVF_MATRIX_H 1

#include <vector>

namespace srvf
{

/**
 * A basic matrix class.
 */
class Matrix
{
public:

  /**
   * Indicates either row-major or column-major ordering.
   * This is used when converting between \c Matrix objects and linear 
   * representations.  In row-major ordering, rows are stored contiguously 
   * in memory, and in column-major ordering, columns are stored 
   * contiguously.
   */
  enum Majorness
  {
    ROW_MAJOR=0,
    COLUMN_MAJOR=1
  };
  
  Matrix ();
  Matrix (size_t rows, size_t cols);
  Matrix (size_t rows, size_t cols, double val);
  Matrix (size_t rows, size_t cols, 
          const double *data, Majorness layout=ROW_MAJOR);
  Matrix (size_t rows, size_t cols, 
          const std::vector<double> &data, Majorness layout=ROW_MAJOR);
  Matrix (const Matrix &A);
  Matrix &operator= (const Matrix &A);

  /** Destructor. */
  ~Matrix () { delete[] data_; }

  /** Returns the number of rows in this \c Matrix. */
  size_t rows() const 
  { return rows_; }

  /** Returns the number of columns in this \c Matrix. */
  size_t cols() const 
  { return cols_; }

  /** Returns the total number of entries in this \c Matrix. */
  size_t size() const 
  { return rows_*cols_; }

  void clear();
  void resize(size_t rows, size_t cols);
  
  /** 
   * Returns the \f$ n^{th} \f$ entry in the matrix.
   * Row-major ordering is used.  For example, in a \c 2x3 \c Matrix, 
   * calling this method with \c n=3 references the entry in row 1, column 0. 
   */ 
  double& operator() (size_t n)
  { return data_[n]; }

  /** 
   * Returns the \f$ n^{th} \f$ entry in the matrix.
   * Row-major ordering is used.  For example, in a \c 2x3 \c Matrix, 
   * calling this method with \c n=3 references the entry in row 1, column 0. 
   */ 
  const double& operator() (size_t n) const
  { return data_[n]; };

  /** Returns the specified entry. */
  double& operator() (size_t r, size_t c)
  { return data_[r*cols_+c]; }

  /** Returns the specified entry. */
  const double& operator() (size_t r, size_t c) const
  { return data_[r*cols_+c]; }

  /** Returns a pointer to the raw data for this \c Matrix. */
  double* data()
  { return data_; }

  /** Returns a \c const pointer to the raw data for this \c Matrix. */
  const double* data() const 
  { return data_; }

  // In-place elementwise operations
  Matrix& operator+= (const Matrix &A);
  Matrix& operator-= (const Matrix &A);
  Matrix& operator*= (const Matrix &A);
  Matrix& operator/= (const Matrix &A);

  Matrix& operator+= (double v);
  Matrix& operator-= (double v);
  Matrix& operator*= (double v);
  Matrix& operator/= (double v);

  // Elementwise operations
  friend Matrix operator+ (const Matrix &A, const Matrix &B);
  friend Matrix operator- (const Matrix &A, const Matrix &B);
  friend Matrix operator* (const Matrix &A, const Matrix &B);
  friend Matrix operator/ (const Matrix &A, const Matrix &B);
  friend Matrix operator+ (const Matrix &A, double v);
  friend Matrix operator- (const Matrix &A, double v);
  friend Matrix operator* (const Matrix &A, double v);
  friend Matrix operator/ (const Matrix &A, double v);

  // Matrix multiplication
  friend Matrix product(const Matrix &A, const Matrix &B);

  // Related matrices
  friend Matrix transpose(const Matrix &A);

private:
  double *data_;
  size_t  rows_;
  size_t  cols_;
};

Matrix operator+ (const Matrix &A, const Matrix &B);
Matrix operator- (const Matrix &A, const Matrix &B);
Matrix operator* (const Matrix &A, const Matrix &B);
Matrix operator/ (const Matrix &A, const Matrix &B);
Matrix operator+ (const Matrix &A, double v);
Matrix operator- (const Matrix &A, double v);
Matrix operator* (const Matrix &A, double v);
Matrix operator/ (const Matrix &A, double v);
Matrix product (const Matrix &A, const Matrix &B);
Matrix transpose(const Matrix &A);

} // namespace srvf


#endif // SRVF_MATRIX_H
