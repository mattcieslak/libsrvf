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
#include "rotselect.h"
#include "pmatch_util.h"
#include "kdtree.h"
#include "rotate.h"
#include "dpnbhd.h"

#include <vector>
#include <map>

namespace srvf
{
namespace pmatch
{


std::vector<Matrix> select_rotations (
  const Srvf &Q1, const Srvf &Q2, 
  const std::vector<double> &tv1, const std::vector<double> &tv2 )
{
  std::map<size_t,size_t> tv1_idx_to_Q1_idx = 
    build_tv_idx_to_Q_idx_map(tv1,Q1);
  std::map<size_t,size_t> tv2_idx_to_Q2_idx = 
    build_tv_idx_to_Q_idx_map(tv2,Q2);

  size_t grid_width = tv1.size();
  size_t grid_height = tv2.size();

  KdTree<Matrix> rot_tree;
  for (size_t ct=1; ct<grid_width; ++ct)
  {
  for (size_t rt=1; rt<grid_height; ++rt)
  {

    for (size_t i=0; i<DP_NBHD_SIZE; ++i)
    {
      size_t cs = ct - dp_nbhd[i][0];
      size_t rs = rt - dp_nbhd[i][1];

      if (cs < ct && rs < rt)
      {
        rot_tree.insert_cond (
          srvf::optimal_rotation (
            Q1, Q2, tv1[cs], tv1[ct], 
            tv2[rs], tv2[rt], 
            tv1_idx_to_Q1_idx[cs], tv2_idx_to_Q2_idx[rs] ),
          0.2 );
      }
    }
  
  }
  }

  return rot_tree.to_vector();
}


} // namespace pmatch
} // namespace srvf
