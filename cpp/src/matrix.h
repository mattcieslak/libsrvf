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

namespace srvf
{

class Matrix
{
public:
  
  Matrix () : data_(0), rows_(0), cols_(0) { }
  Matrix (int rows, int cols);
  Matrix (int rows, int cols, double val);
  Matrix (int rows, int cols, double *data);
  Matrix (const Matrix &A);
  Matrix &operator= (const Matrix &A);
  ~Matrix ();

  int rows() const 
  { return rows_; }

  int cols() const 
  { return cols_; }

  int size() const 
  { return rows_*cols_; }

  void clear();
  void resize(int rows, int cols);
  
  double& operator() (int n)
  { return data_[n]; }

  const double& operator() (int n) const
  { return data_[n]; };

  double& operator() (int r, int c)
  { return data_[r*cols_+c]; }

  const double& operator() (int r, int c) const
  { return data_[r*cols_+c]; }

  double* data()
  { return data_; }

  const double* data() const 
  { return (const double*)data_; }

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
  double   *data_;
  int       rows_;
  int       cols_;
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
