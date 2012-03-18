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
#ifndef SRVF_SRVF_H
#define SRVF_SRVF_H 1

#include <vector>

#include "matrix.h"
#include "plf.h"

namespace srvf {

class Plf;
  
class Srvf
{
public:

  void evaluate(double t, Matrix &result);
  void evaluate(const Matrix &tv, Matrix &result);

  double l2_norm();
  double l2_product(const Srvf &Q);
  double l2_distance(const Srvf &Q);
  double sphere_distance(const Srvf &Q);

  void linear_combine(const Srvf &Q, double w1, double w2);
  void refine(const Matrix &tv);
  void gamma_act(const Plf &gamma);

  Plf to_plf();

  friend double l2_product(const Srvf &Q1, const Srvf &Q2);
  friend double l2_distance(const Srvf &Q1, const Srvf &Q2);
  friend double sphere_distance(const Srvf &Q1, const Srvf &Q2);
  friend Srvf   linear_combination(const Srvf &Q1, const Srvf &Q2, 
                  double w1, double w2);
  friend Srvf   refinement(const Srvf &Q, const Matrix &tv);
  friend Srvf   gamma_action(const Srvf &Q, const Plf &gamma);

private:
  Matrix samps_;
  Matrix params_;
};

} // namespace srvf

#endif // SRVF_SRVF_H
