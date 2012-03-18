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
#include "srvf.h"

namespace srvf {

class Srvf;

class Plf
{
public:
  
  Plf() { }
  Plf(const Matrix &samples);
  Plf(const Matrix &samples, const Matrix &parameters);
 
  void evaluate(double t, Matrix &result);
  void evaluate(const Matrix &tv, Matrix &result);
  void preimages(const Matrix &tv, Matrix &result);
  double arc_length() const;

  void translate(const Matrix &v);
  void rotate(const Matrix &R);
  void scale(double s);
  void linear_combine(const Plf &F, double w1, double w2);
  void compose(const Plf &F);
  void invert();

  Srvf to_srvf() const;

  friend Plf translation(const Plf &F, const Matrix &v);
  friend Plf rotation(const Plf &F, const Matrix &R);
  friend Plf scaling(const Plf &F, double s);
  friend Plf linear_combination(const Plf &F1, const Plf &F2, 
               double w1, double w2);
  friend Plf composition(const Plf &F1, const Plf &F2);
  friend Plf inverse(const Plf &F);

private:
  Matrix samps_;
  Matrix params_;
};

} // namespace srvf

#endif // SRVF_PLF_H
