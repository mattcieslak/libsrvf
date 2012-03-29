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
#include <limits>
#include <cmath>
#include <vector>
#include <stack>

#include "reparam.h"
#include "plf.h"
#include "srvf.h"
#include "pointset.h"
#include "matrix.h"


#ifndef DP_NBHD_DIM
#define DP_NBHD_DIM 7 
#endif

#if DP_NBHD_DIM == 17
size_t dp_nbhd[][2] = {
{  1,  1 }, {  1,  2 }, {  1,  3 }, {  1,  4 }, {  1,  5 }, {  1,  6 }, {  1,  7 }, {  1,  8 }, {  1,  9 }, {  1, 10 }, 
{  1, 11 }, {  1, 12 }, {  1, 13 }, {  1, 14 }, {  1, 15 }, {  1, 16 }, {  1, 17 }, {  2,  1 }, {  2,  3 }, {  2,  5 }, 
{  2,  7 }, {  2,  9 }, {  2, 11 }, {  2, 13 }, {  2, 15 }, {  2, 17 }, {  3,  1 }, {  3,  2 }, {  3,  4 }, {  3,  5 }, 
{  3,  7 }, {  3,  8 }, {  3, 10 }, {  3, 11 }, {  3, 13 }, {  3, 14 }, {  3, 16 }, {  3, 17 }, {  4,  1 }, {  4,  3 }, 
{  4,  5 }, {  4,  7 }, {  4,  9 }, {  4, 11 }, {  4, 13 }, {  4, 15 }, {  4, 17 }, {  5,  1 }, {  5,  2 }, {  5,  3 }, 
{  5,  4 }, {  5,  6 }, {  5,  7 }, {  5,  8 }, {  5,  9 }, {  5, 11 }, {  5, 12 }, {  5, 13 }, {  5, 14 }, {  5, 16 }, 
{  5, 17 }, {  6,  1 }, {  6,  5 }, {  6,  7 }, {  6, 11 }, {  6, 13 }, {  6, 17 }, {  7,  1 }, {  7,  2 }, {  7,  3 }, 
{  7,  4 }, {  7,  5 }, {  7,  6 }, {  7,  8 }, {  7,  9 }, {  7, 10 }, {  7, 11 }, {  7, 12 }, {  7, 13 }, {  7, 15 }, 
{  7, 16 }, {  7, 17 }, {  8,  1 }, {  8,  3 }, {  8,  5 }, {  8,  7 }, {  8,  9 }, {  8, 11 }, {  8, 13 }, {  8, 15 }, 
{  8, 17 }, {  9,  1 }, {  9,  2 }, {  9,  4 }, {  9,  5 }, {  9,  7 }, {  9,  8 }, {  9, 10 }, {  9, 11 }, {  9, 13 }, 
{  9, 14 }, {  9, 16 }, {  9, 17 }, { 10,  1 }, { 10,  3 }, { 10,  7 }, { 10,  9 }, { 10, 11 }, { 10, 13 }, { 10, 17 }, 
{ 11,  1 }, { 11,  2 }, { 11,  3 }, { 11,  4 }, { 11,  5 }, { 11,  6 }, { 11,  7 }, { 11,  8 }, { 11,  9 }, { 11, 10 }, 
{ 11, 12 }, { 11, 13 }, { 11, 14 }, { 11, 15 }, { 11, 16 }, { 11, 17 }, { 12,  1 }, { 12,  5 }, { 12,  7 }, { 12, 11 }, 
{ 12, 13 }, { 12, 17 }, { 13,  1 }, { 13,  2 }, { 13,  3 }, { 13,  4 }, { 13,  5 }, { 13,  6 }, { 13,  7 }, { 13,  8 }, 
{ 13,  9 }, { 13, 10 }, { 13, 11 }, { 13, 12 }, { 13, 14 }, { 13, 15 }, { 13, 16 }, { 13, 17 }, { 14,  1 }, { 14,  3 }, 
{ 14,  5 }, { 14,  9 }, { 14, 11 }, { 14, 13 }, { 14, 15 }, { 14, 17 }, { 15,  1 }, { 15,  2 }, { 15,  4 }, { 15,  7 }, 
{ 15,  8 }, { 15, 11 }, { 15, 13 }, { 15, 14 }, { 15, 16 }, { 15, 17 }, { 16,  1 }, { 16,  3 }, { 16,  5 }, { 16,  7 }, 
{ 16,  9 }, { 16, 11 }, { 16, 13 }, { 16, 15 }, { 16, 17 }, { 17,  1 }, { 17,  2 }, { 17,  3 }, { 17,  4 }, { 17,  5 }, 
{ 17,  6 }, { 17,  7 }, { 17,  8 }, { 17,  9 }, { 17, 10 }, { 17, 11 }, { 17, 12 }, { 17, 13 }, { 17, 14 }, { 17, 15 }, 
{ 17, 16 } };
#define DP_NBHD_COUNT 191

#elif DP_NBHD_DIM == 12
size_t dp_nbhd[][2] = {
{  1,  1 }, {  1,  2 }, {  1,  3 }, {  1,  4 }, {  1,  5 }, {  1,  6 }, {  1,  7 }, {  1,  8 }, {  1,  9 }, {  1, 10 }, 
{  1, 11 }, {  1, 12 }, {  2,  1 }, {  2,  3 }, {  2,  5 }, {  2,  7 }, {  2,  9 }, {  2, 11 }, {  3,  1 }, {  3,  2 }, 
{  3,  4 }, {  3,  5 }, {  3,  7 }, {  3,  8 }, {  3, 10 }, {  3, 11 }, {  4,  1 }, {  4,  3 }, {  4,  5 }, {  4,  7 }, 
{  4,  9 }, {  4, 11 }, {  5,  1 }, {  5,  2 }, {  5,  3 }, {  5,  4 }, {  5,  6 }, {  5,  7 }, {  5,  8 }, {  5,  9 }, 
{  5, 11 }, {  5, 12 }, {  6,  1 }, {  6,  5 }, {  6,  7 }, {  6, 11 }, {  7,  1 }, {  7,  2 }, {  7,  3 }, {  7,  4 }, 
{  7,  5 }, {  7,  6 }, {  7,  8 }, {  7,  9 }, {  7, 10 }, {  7, 11 }, {  7, 12 }, {  8,  1 }, {  8,  3 }, {  8,  5 }, 
{  8,  7 }, {  8,  9 }, {  8, 11 }, {  9,  1 }, {  9,  2 }, {  9,  4 }, {  9,  5 }, {  9,  7 }, {  9,  8 }, {  9, 10 }, 
{  9, 11 }, { 10,  1 }, { 10,  3 }, { 10,  7 }, { 10,  9 }, { 10, 11 }, { 11,  1 }, { 11,  2 }, { 11,  3 }, { 11,  4 }, 
{ 11,  5 }, { 11,  6 }, { 11,  7 }, { 11,  8 }, { 11,  9 }, { 11, 10 }, { 11, 12 }, { 12,  1 }, { 12,  5 }, { 12,  7 }, 
{ 12, 11 } };
#define DP_NBHD_COUNT 91

#elif DP_NBHD_DIM == 10

size_t dp_nbhd[][2] = {
{  1,  1 }, {  1,  2 }, {  1,  3 }, {  1,  4 }, {  1,  5 }, {  1,  6 }, {  1,  7 }, {  1,  8 }, {  1,  9 }, {  1, 10 }, 
{  2,  1 }, {  2,  3 }, {  2,  5 }, {  2,  7 }, {  2,  9 }, {  3,  1 }, {  3,  2 }, {  3,  4 }, {  3,  5 }, {  3,  7 }, 
{  3,  8 }, {  3, 10 }, {  4,  1 }, {  4,  3 }, {  4,  5 }, {  4,  7 }, {  4,  9 }, {  5,  1 }, {  5,  2 }, {  5,  3 }, 
{  5,  4 }, {  5,  6 }, {  5,  7 }, {  5,  8 }, {  5,  9 }, {  6,  1 }, {  6,  5 }, {  6,  7 }, {  7,  1 }, {  7,  2 }, 
{  7,  3 }, {  7,  4 }, {  7,  5 }, {  7,  6 }, {  7,  8 }, {  7,  9 }, {  7, 10 }, {  8,  1 }, {  8,  3 }, {  8,  5 }, 
{  8,  7 }, {  8,  9 }, {  9,  1 }, {  9,  2 }, {  9,  4 }, {  9,  5 }, {  9,  7 }, {  9,  8 }, {  9, 10 }, { 10,  1 }, 
{ 10,  3 }, { 10,  7 }, { 10,  9 } };
#define DP_NBHD_COUNT 63

#elif DP_NBHD_DIM == 7

size_t dp_nbhd[][2] = {
{  1,  1 }, {  1,  2 }, {  1,  3 }, {  1,  4 }, {  1,  5 }, {  1,  6 }, {  1,  7 }, {  2,  1 }, {  2,  3 }, {  2,  5 }, 
{  2,  7 }, {  3,  1 }, {  3,  2 }, {  3,  4 }, {  3,  5 }, {  3,  7 }, {  4,  1 }, {  4,  3 }, {  4,  5 }, {  4,  7 }, 
{  5,  1 }, {  5,  2 }, {  5,  3 }, {  5,  4 }, {  5,  6 }, {  5,  7 }, {  6,  1 }, {  6,  5 }, {  6,  7 }, {  7,  1 }, 
{  7,  2 }, {  7,  3 }, {  7,  4 }, {  7,  5 }, {  7,  6 } };
#define DP_NBHD_COUNT 35

#else // DP_NBHD_DIM = 6
size_t dp_nbhd[][2] = {
{  1,  1 }, {  1,  2 }, {  1,  3 }, {  1,  4 }, {  1,  5 }, {  1,  6 }, {  2,  1 }, {  2,  3 }, {  2,  5 }, {  3,  1 }, 
{  3,  2 }, {  3,  4 }, {  3,  5 }, {  4,  1 }, {  4,  3 }, {  4,  5 }, {  5,  1 }, {  5,  2 }, {  5,  3 }, {  5,  4 }, 
{  5,  6 }, {  6,  1 }, {  6,  5 } };
#define DP_NBHD_COUNT 23

#endif  /* DP_NBHD_DIM */



namespace srvf
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
      for (size_t k=0; k<DP_NBHD_COUNT; ++k)
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
    samps(i,0) = Q2.params()[i2];
  }

  return Plf(samps, params);
}

} // namespace srvf

