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
#ifndef PMATCH_UTIL_H
#define PMATCH_UTIL_H 1

#include <srvf.h>

#include <vector>
#include <map>
#include <cstddef>


namespace srvf
{
namespace pmatch
{


/**
 * Build a map that sends each index \c i of \a tv to the index of the 
 * parameter subinterval of \a Q which contains \c tv[i].
 *
 * These maps are used in \c calculate_edge_weights().
 */
std::map<size_t,size_t> build_tv_idx_to_Q_idx_map (
  const std::vector<double> &tv, const srvf::Srvf &Q );
  

} // namespace pmatch
} // namespace srvf

#endif // PMATCH_UTIL_H
