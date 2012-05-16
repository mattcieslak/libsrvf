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
#include "functions.h"

#include "srvf.h"
#include "plf.h"

#include <vector>


namespace srvf
{

namespace functions
{


std::vector<Plf> optimal_reparam(const Srvf &Q1, const Srvf &Q2)
{
  // TODO
  std::vector<Plf> res;
  return res;
}


Srvf karcher_mean(const std::vector<Srvf> &Qs, 
                  double tol, size_t max_iters)
{
  // TODO
  Srvf res;
  return res;
}

} // namespace srvf::functions

} // namespace srvf
