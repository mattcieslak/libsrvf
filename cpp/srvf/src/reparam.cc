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
#include "reparam.h"
#include "dpnbhd.h"
#include "pointset.h"
#include "matrix.h"

#include <limits>
#include <cmath>
#include <stack>


namespace srvf
{

namespace opencurves
{


/**
 * Computes an optimal reparametrization of \a Q2 relative to \a Q1.
 */
Plf optimal_reparam (const Srvf &Q1, const Srvf &Q2)
{
  size_t npts1 = Q1.ncp();
  size_t npts2 = Q2.ncp();

  // DP grid.  Columns correspond to parameters of Q1, and 
  // rows correspond to parameters of Q2.
  Matrix grid(npts2, npts1, std::numeric_limits<double>::infinity());

  // Predecessor matrix.  V=preds[r*npts1+c] holds predecessor of gridpoint 
  // at row r and column c, encoded so that V/npts1 gives the row of the 
  // predecessor and V%npts1 gives the column.
  std::vector<int> preds(npts1*npts2);
  

  // Main DP loop
  grid(0,0) = 0.0;
  for (size_t r2=1; r2<npts2; ++r2)
  {
    for (size_t c2=1; c2<npts1; ++c2)
    {
      for (size_t k=0; k<DP_NBHD_SIZE; ++k)
      {
        size_t r1 = r2 - dp_nbhd[k][0];
        size_t c1 = c2 - dp_nbhd[k][1];

        // Check for underflow
        if (r1>r2 || c1>c2) continue;
        
        double w = match_cost(Q1, c1, c2, Q2, r1, r2);
        double cand_cost = grid(r1,c1) + w;
        if (cand_cost < grid(r2,c2))
        {
          grid(r2,c2) = cand_cost;
          preds[r2*npts1+c2] = r1*npts1+c1;
        }
      }
    }
  }

  // Reconstruct the path by following the trail of predecessors 
  // from the top-right corner of the grid to the bottom-left corner.
  std::stack<size_t> path;
  path.push((npts2-1)*npts1 + (npts1-1));
  while (path.top() != 0)
  {
    path.push(preds[path.top()]);
  }

  // Reconstruct the piecewise-linear reparametrization function.
  std::vector<double> params(path.size());
  Pointset samps(1, path.size());
  for (size_t i=0; !path.empty(); ++i)
  {
    size_t i1 = path.top() % npts1;
    size_t i2 = path.top() / npts1;
    path.pop();
    params[i] = Q1.params()[i1];
    samps[i][0] = Q2.params()[i2];
  }

  return Plf(samps, params);
}


} // namespace opencurves

} // namespace srvf
