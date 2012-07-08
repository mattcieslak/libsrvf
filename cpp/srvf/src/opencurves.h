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
#ifndef SRVF_OPENCURVES_H
#define SRVF_OPENCURVES_H 1

#include "srvf.h"

#include <cstddef>
#include <vector>


namespace srvf
{

class Srvf;

namespace opencurves
{

Srvf karcher_mean(const std::vector<Srvf> &Qs, 
                  double tol=1e-3, 
                  size_t max_iters=0, 
                  bool optimize_rots=true, 
                  bool optimize_reparams=true);

Srvf shooting_vector(const Srvf &Q1, const Srvf &Q2);

} //namespace opencurves

} // namespace srvf

#endif // SRVF_OPENCURVES_H
