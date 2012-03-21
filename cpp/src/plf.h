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
#ifndef SRVF_PLF_H
#define SRVF_PLF_H 1

#include "matrix.h"

namespace srvf {

class Plf
{
public:
  
  Plf() { }
  Plf(const Matrix &samples);
  Plf(const Matrix &samples, const Matrix &parameters);

  Matrix &samps() { return samps_; }
  const Matrix &samps() const { return samps_; }

  Matrix &params() { return params_; }
  const Matrix &params() const { return params_; }

  int dim() const { return samps_.rows(); }
  int nsamps() const { return samps_.cols(); }
 
  void evaluate(double t, Matrix &result) const;
  void evaluate(const Matrix &tv, Matrix &result) const;
  void preimages(const Matrix &tv, Matrix &result) const;
  double arc_length() const;

  void translate(const Matrix &v);
  void rotate(const Matrix &R);
  void scale(double s);

  friend Plf linear_combination(const Plf &F1, const Plf &F2, 
                                double w1, double w2);
  friend Plf composition(const Plf &F1, const Plf &F2);
  friend Plf inverse(const Plf &F);

private:
  Matrix samps_;
  Matrix params_;
};

Plf  linear_combination(const Plf &F1, const Plf &F2, double w1, double w2);
Plf  composition(const Plf &F1, const Plf &F2);
Plf  inverse(const Plf &F);

} // namespace srvf

#endif // SRVF_PLF_H
