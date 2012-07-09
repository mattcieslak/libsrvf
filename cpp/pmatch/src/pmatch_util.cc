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
#include "pmatch_util.h"


namespace srvf
{
namespace pmatch
{


std::map<size_t,size_t> build_tv_idx_to_Q_idx_map (
  const std::vector<double> &tv, const srvf::Srvf &Q )
{
  std::map<size_t,size_t> result;

  for (size_t tvi=0, Qi=0; tvi+1<tv.size(); ++tvi)
  {
    result[tvi] = Qi;
    while (Qi+2 < Q.ncp() && tv[tvi+1] > (Q.params()[Qi+1]-1e-5))
      ++Qi;
  }
  result[tv.size()-1] = Q.ncp()-2;

  return result;
}

} // namespace pmatch
} // namespace srvf
