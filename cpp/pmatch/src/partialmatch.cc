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

#include "partialmatch.h"
#include "paretoset.h"

#include <srvf.h>
#include <dpnbhd.h>

#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>

namespace srvf
{

namespace pmatch
{


/**
 * Computes the weight of an edge in the matching graph.
 */
double edge_weight (
  const srvf::Srvf &Q1, const srvf::Srvf &Q2, 
  const std::vector<double> &tv1, const std::vector<double> &tv2, 
  size_t tv1_i1, size_t tv1_i2, 
  size_t tv2_i1, size_t tv2_i2, 
  size_t Q1_idx_start, size_t Q2_idx_start )
{
  double a = tv1[tv1_i1], b = tv1[tv1_i2];
  double c = tv2[tv2_i1], d = tv2[tv2_i2];

  double m = (d-c) / (b-a);
  double rm = sqrt(m);
  
  size_t Q1_idx = Q1_idx_start;
  size_t Q2_idx = Q2_idx_start;
  double t1 = a;
  double t2 = c;
  double result = 0.0;


  while (t1 < (b-1e-5) && t2 < (d-1e-5))
  {
    double dx1 = std::min(b, Q1.params()[Q1_idx+1]) - t1;
    double dy1 = m * dx1;
    double dy2 = std::min(d, Q2.params()[Q2_idx+1]) - t2;
    double dx2 = dy2 / m;
    
    double dQi = 0.0;
    for (size_t j=0; j<Q1.dim(); ++j)
    {
      double dQij = Q1.samps()[Q1_idx][j] - rm * Q2.samps()[Q2_idx][j];
      dQi += dQij * dQij;
    }

    if ( fabs(dx1 - dx2) < 1e-5 )
    {
      result += dx1 * dQi;
      t1 += dx1;
      t2 += dy1;
      ++Q1_idx;
      ++Q2_idx;
    }
    else if (dx1 < dx2)
    {
      result += dx1 * dQi;
      t1 += dx1;
      t2 += dy1;
      ++Q1_idx;
    }
    else
    {
      result += dx2 * dQi;
      t1 += dx2;
      t2 += dy2;
      ++Q2_idx;
    }
  }

  return result;
}


/**
 * Build a map that sends each index \c i of \a tv to the index of the 
 * parameter subinterval of \a Q which contains \c tv[i].
 *
 * These maps are used in \c calculate_edge_weights().
 */
static std::map<size_t,size_t> build_tv_idx_to_Q_idx_map_ (
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
  

/**
 * Calculates the weights of all edges in the matching graph.
 */
MatchingGraph calculate_edge_weights (
  const srvf::Srvf &Q1, const srvf::Srvf &Q2, 
  const std::vector<double> &tv1, const std::vector<double> &tv2 )
{
  std::map<size_t,size_t> tv1_idx_to_Q1_idx = 
    build_tv_idx_to_Q_idx_map_(tv1,Q1);
  std::map<size_t,size_t> tv2_idx_to_Q2_idx = 
    build_tv_idx_to_Q_idx_map_(tv2,Q2);

  MatchingGraph result(tv1.size(), tv2.size());

  for (size_t ct=1; ct<tv1.size(); ++ct)
  {
    for (size_t rt=1; rt<tv2.size(); ++rt)
    {
      for (size_t i=0; i<DP_NBHD_SIZE; ++i)
      {
        size_t cs = ct - srvf::dp_nbhd[i][0];
        size_t rs = rt - srvf::dp_nbhd[i][1];
        if (cs < tv1.size() && rs < tv2.size()){
          double w = edge_weight (
            Q1, Q2, tv1, tv2, 
            cs, ct, rs, rt, 
            tv1_idx_to_Q1_idx[cs], 
            tv2_idx_to_Q2_idx[rs] );

          if (w < 0.0)
            w = 0.0;

          result(cs, rs, ct, rt) = w;
        }
      }
    }
  }

  return result;
}


/**
 * Use the Floyd-Warshall algorithm to solve the all-pairs-shortest 
 * path problem in \a G.
 */
void calculate_match_scores (MatchingGraph &G)
{
  for (size_t cm=1; cm+1<G.width(); ++cm)
  {
  for (size_t rm=1; rm+1<G.height(); ++rm)
  {

    for (size_t ct=cm+1; ct<G.width(); ++ct)
    {
    for (size_t rt=rm+1; rt<G.height(); ++rt)
    {

      for (size_t cs=0; cs<cm; ++cs)
      {
      for (size_t rs=0; rs<rm; ++rs)
      {
        double cand_w = G(cs, rs, cm, rm) + G(cm, rm, ct, rt);
        if (cand_w < G(cs, rs, ct, rt))
          G(cs, rs, ct, rt) = cand_w;
      }
      }
    
    }
    }

  }
  }
}


ParetoSet find_matches (
  const srvf::Srvf &Q1, const srvf::Srvf &Q2, 
  bool do_rots, size_t grid_width, size_t grid_height, 
  size_t nbuckets, double bucket_thresh )
{
  std::vector<double> tv1;
  std::vector<double> tv2;

  if (grid_width == 0)
    tv1 = Q1.params();
  else
    tv1 = srvf::util::linspace(Q1.domain_lb(), Q1.domain_ub(), grid_width);

  if (grid_height == 0)
    tv2 = Q2.params();
  else
    tv2 = srvf::util::linspace(Q2.domain_lb(), Q2.domain_ub(), grid_height);

  if (nbuckets == 0)
  {
    // Default: nbuckets = number of possible match lengths.
    size_t min_chunks = 2;
    size_t max_chunks = (tv1.size()-1) + (tv2.size()-1);
    nbuckets = max_chunks - min_chunks + 1;
  }

  MatchingGraph G = calculate_edge_weights(Q1, Q2, tv1, tv2);
  calculate_match_scores(G);

  ParetoSet S(nbuckets, bucket_thresh);
  for (size_t ct=1; ct<grid_width; ++ct)
  {
  for (size_t rt=1; rt<grid_height; ++rt)
  {

    for (size_t cs=0; cs<ct; ++cs)
    {
    for (size_t rs=0; rs<rt; ++rs)
    {
      S.insert( PartialMatch( 
        tv1[cs], tv1[ct], 
        tv2[rs], tv2[rt],
        G(cs, rs, ct, rt) ) );
    }
    }
  
  }
  }

  return S;
}


} // namespace pmatch
} // namespace srvf
