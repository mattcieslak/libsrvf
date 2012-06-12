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
#include <utility>
#include <stack>
#include <map>
#include <algorithm>
#include <iterator>
#include <iostream>


namespace srvf
{

namespace functions
{

static inline double
time_up_(const Srvf &Q, size_t i1, size_t i2)
{
  double v=0.0;
  for (size_t i=i1; i<i2; i+=2)
  {
    v += Q.params()[i+1] - Q.params()[i];
  }
  return v;
}


static inline bool
is_peak_(const Srvf &Q, size_t idx)
{
  return ( (idx>0 && Q.samps()[idx-1][0] > 0.0) ||
           (idx+1<Q.ncp() && Q.samps()[idx][0] < 0.0) );
}


static inline double
edge_variation_(const Srvf &Q, size_t i1, size_t i2)
{
  double v=0.0;
  for (size_t i=i1; i<i2; i+=2)
  {
    double Qi = Q.samps()[i][0];
    v += Qi * fabs(Qi) * (Q.params()[i+1]-Q.params()[i]);
  }
  return v;
}


static inline double 
edge_score_(const Srvf &Q1, const Srvf &Q2, 
  size_t sc, size_t sr, size_t tc, size_t tr)
{
  //double v1 = edge_variation_(Q1,sc,tc);
  //double v2 = edge_variation_(Q2,sr,tr);
  double v1 = time_up_(Q1,sc,tc);
  double v2 = time_up_(Q2,sr,tr);
  return sqrt(v1*v2);
}


// Builds the pieces of gamma1 and gamma2 corresponding to 
// the edge (sc,sr)-->(tc,tr) in the matching graph.  The 
// output is appended to G1samps and G2samps.
static inline void 
build_gamma_segment_(const Srvf &Q1, const Srvf &Q2, 
  size_t sc, size_t sr, size_t tc, size_t tr, 
  std::vector<double> &G1samps, std::vector<double> &G2samps)
{
  double v1 = time_up_(Q1, sc, tc);
  double v2 = time_up_(Q2, sr, tr);
  double R = v2 / v1;
  double u1 = 0.0;
  double u2 = 0.0;
  double t1 = Q1.params()[sc];
  double t2 = Q2.params()[sr];
  size_t c = sc;
  size_t r = sr;

  while (c<tc && r<tr)
  {
    if (c+1 < tc || r+1 < tr)
    {
      double p1 = (u1 + (Q1.params()[c+1] - t1)) / v1;
      double p2 = (u2 + (Q2.params()[r+1] - t2)) / v2;

      if (fabs(p1-p2) > 1e-4)
      {
        if (p1 < p2)
        {
          double lambda = t2 + R*(Q1.params()[c+1]-t1);
          G1samps.push_back(Q1.params()[c+1]);
          G1samps.push_back(Q1.params()[c+2]);
          G2samps.push_back(lambda);
          G2samps.push_back(lambda);
          u1 += (Q1.params()[c+1]-t1);
          u2 += (lambda-t2);
          t1 = Q1.params()[c+2];
          t2 = lambda;
          c += 2;
        }
        else
        {
          double lambda = t1 + (Q2.params()[r+1]-t2) / R;
          G1samps.push_back(lambda);
          G1samps.push_back(lambda);
          G2samps.push_back(Q2.params()[r+1]);
          G2samps.push_back(Q2.params()[r+2]);
          u1 += (lambda-t1);
          u2 += (Q2.params()[r+1]-t2);
          t1 = lambda;
          t2 = Q2.params()[r+2];
          r += 2;
        }
      }
      else
      {
        G1samps.push_back(Q1.params()[c+1]);
        G1samps.push_back(Q1.params()[c+2]);
        G2samps.push_back(Q2.params()[r+1]);
        G2samps.push_back(Q2.params()[r+2]);
        u1 += (Q1.params()[c+1]-t1);
        u2 += (Q2.params()[r+1]-t2);
        t1 = Q1.params()[c+2];
        t2 = Q2.params()[r+2];
        c += 2;
        r += 2;
      }
    }
    else
    {
      G1samps.push_back(Q1.params()[tc]);
      G2samps.push_back(Q2.params()[tr]);
      break;
    } 
  }
}


static void calculate_scores_(const Srvf &Q1, const Srvf &Q2, 
  std::map<match_vertex_t,double> &score, 
  std::map<match_vertex_t,match_vertex_t> &pred)
{
  size_t n1 = Q1.ncp();
  size_t n2 = Q2.ncp();

  // Dynamic programming step: flood-fill the score matrix.
  for (size_t tr=1; tr<n2; ++tr)
  {
    size_t tc=1;
    if (is_peak_(Q1,tc) != is_peak_(Q2,tr)) ++tc;

    for (/* NOOP */; tc<n1; tc+=2)
    {
      match_vertex_t tv(tc,tr);
      double best_score = -1.0;
      match_vertex_t best_pred;
      
      for (size_t sr=1-(tr%2); sr<tr; sr+=2)
      {
        for (size_t sc=1-(tc%2); sc<tc; sc+=2)
        {
          match_vertex_t sv(sc,sr);
          double w = edge_score_(Q1,Q2,sc,sr,tc,tr);
          double s = score[sv] + w;
          if (s > best_score)
          {
            best_score = s;
            best_pred = sv;
          }
        }
      }

      if (best_score > -1.0)
      {
        score[tv] = best_score;
        pred[tv] = best_pred;
      }
    }
  }

  //for (size_t r=n2-1; r<n2; --r)
  //{
  //  for (size_t c=0; c<n1; ++c)
  //  {
  //    if (score.find(match_vertex_t(c,r)) != score.end())
  //    {
  //      std::cout << "[" << score[match_vertex_t(c,r)] << " : " << 
  //        "(" << pred[match_vertex_t(c,r)].first << "," <<
  //        pred[match_vertex_t(c,r)].second << ")-->" << 
  //        "(" << c << "," << r << ")] ";
  //    }
  //    else
  //    {
  //      std::cout << "[(" << c << "," << r << ") xxxxxxxxxxxxxxxxxx] ";
  //    }
  //  }
  //  std::cout << std::endl;
  //}
  //std::cout << std::endl;
}


bool fuzzy_less_(double a, double b)
{
  return a < b-1e-6;
}

// Computes reparametrizations eta, gamma_1, gamma_2, ..., gamma_k such 
// that Mu*eta is optimally matched to Qs[i]*gamma_i for i=1,2,...k.
// The result is a vector of Plf's which contains eta, then gamma_1, 
// gamma_2, ..., and gamma_k.
// 
// IMPORTANT:  Mu and all elements of Qs must already be normalized, so 
// that they only take on values in {1,-1}.
static std::vector<Plf> 
groupwise_match_(const Srvf &Mu, const std::vector<Srvf> &Qs)
{
  size_t nfuncs = Qs.size();
  std::vector<double> X;
  std::vector<std::vector<double> > Xs(nfuncs);
  std::vector<std::vector<double> > XCs(nfuncs);
  std::vector<std::vector<double> > XCPs(nfuncs);
  std::vector<Pointset> Ys;
  std::vector<Plf> Es;
  std::vector<Plf> Gs;

  // Compute the individual matchings between Mu and each Qs[i]
  for (size_t i=0; i<nfuncs; ++i)
  {
    std::vector<Plf> EGcur = optimal_reparam(Mu, Qs[i]);
    Es.push_back(EGcur[0]);
    Gs.push_back(EGcur[1]);
  }

  // Find the extra points on Mu resulting from each of these matchings
  for (size_t i=0; i<nfuncs; ++i)
  {
    std::vector<double> Eivals = Es[i].samps().to_vector();
    std::set_difference(Eivals.begin(), Eivals.end(), 
      Mu.params().begin(), Mu.params().end(), 
      std::back_insert_iterator<std::vector<double> > (Xs[i]), 
      fuzzy_less_);
    X.insert(X.end(), Xs[i].begin(), Xs[i].end());
  }

  // The common parameters for all reparametrizations in the result
  std::vector<double> uparams = srvf::util::linspace(0.0, 1.0, 
    Mu.ncp() + X.size() );

  for (size_t i=0; i<nfuncs; ++i)
  {
    for (size_t j=0; j<nfuncs; ++j)
    {
      if (i==j) continue;
      XCs[j].insert(XCs[j].end(), Xs[i].begin(), Xs[i].end());
    }
  }

  for (size_t i=0; i<nfuncs; ++i)
  {
    //std::sort(XCs[i].begin(), XCs[i].end());
    //XCPs[i] = Es[i].preimages(XCs[i]);
    //Ys.push_back(Gs[i].evaluate(XCPs[i]));
  }
}


std::vector<Plf> optimal_reparam(const Srvf &Q1, const Srvf &Q2)
{
  std::map<match_vertex_t,double> score;
  std::map<match_vertex_t,match_vertex_t> pred;

  size_t n1 = Q1.ncp();
  size_t n2 = Q2.ncp();

  // All paths start at start_vertex.  If both functions start with a peak, 
  // or both functions start with a valley, then the start vertex is (0,0).  
  // Otherwise, we make a dummy vertex, and add zero-score edges from it 
  // to (1,0) and (0,1).
  match_vertex_t start_vertex;
  if (is_peak_(Q1,0) == is_peak_(Q2,0))
  {
    start_vertex = match_vertex_t(0,0);
    score[start_vertex] = 0.0;
  }
  else
  {
    start_vertex = match_vertex_t(n1,n2);  // dummy vertex: invalid indices
    score[start_vertex] = 0.0;
    score[match_vertex_t(0,1)] = 0.0;
    pred[match_vertex_t(0,1)] = start_vertex;
    score[match_vertex_t(1,0)] = 0.0;
    pred[match_vertex_t(1,0)] = start_vertex;
  }

  // Call the main DP routine
  calculate_scores_(Q1, Q2, score, pred);

  // Select the end vertex.  If both functions end with a peak or both 
  // functions end with a valley, the end vertex is (n1-1,n2-1).  Otherwise, 
  // the end vertex is (n1-2,n2-1) or (n1-1,n2-2), depending on which one 
  // has the higher score.
  match_vertex_t end_vertex(n1-1,n2-1);
  if (is_peak_(Q1,n1-1) != is_peak_(Q2,n2-1))
  {
    match_vertex_t end_cand1(n1-2,n2-1);
    match_vertex_t end_cand2(n1-1,n2-2);

    if (score[end_cand1] > score[end_cand2])
      end_vertex = end_cand1;
    else
      end_vertex = end_cand2;
  }

  // Reconstruct the optimal path from start_vertex to end_vertex.
  std::stack<match_vertex_t> path;
  path.push(end_vertex);
  while (path.top() != start_vertex) path.push(pred[path.top()]);

  // Discard the first vertex.  It will either be a dummy vertex, 
  // or it will be (0,0).  In either case, we don't need it any more.
  path.pop();

  // Translate the path into a pair of piecewise-linear reparametrization
  // functions gamma1 and gamma2.
  std::vector<double> G1samps(1,0.0);
  std::vector<double> G2samps(1,0.0);
  size_t sr = path.top().first; 
  size_t sc = path.top().second;

  // Handle case where one function begins with a peak and the other 
  // begins with a valley.
  if (sr == 0)
  {
    G1samps.push_back(Q1.params()[sc]);
    G2samps.push_back(Q2.params()[0]);
  }
  else if (sc == 0)
  {
    G1samps.push_back(Q1.params()[0]);
    G2samps.push_back(Q2.params()[sr]);
  }

  while (!path.empty())
  {
    size_t tc = path.top().first;
    size_t tr = path.top().second;

    build_gamma_segment_(Q1,Q2,sc,sr,tc,tr,G1samps,G2samps);

    sr = tr;
    sc = tc;
    path.pop();
  }

  // Handle case where one function ends with a peak and the other ends 
  // with a valley.
  if (sc < n1-1 || sr < n2-1)
  {
    G1samps.push_back(Q1.domain_ub());
    G2samps.push_back(Q2.domain_ub());
  }

  // The parameters for gamma1 and gamma2 are arbitrary, so we take 
  // uniformly-spaced parameters in [0,1].
  std::vector<double> G1params = 
    util::linspace(0.0, 1.0, G1samps.size());
  std::vector<double> G2params = 
    util::linspace(0.0, 1.0, G2samps.size());

  std::vector<Plf> res;
  res.push_back(Plf(Pointset(1,G1samps.size(),G1samps), G1params));
  res.push_back(Plf(Pointset(1,G2samps.size(),G2samps), G2params));
  return res;
}


Srvf karcher_mean(const std::vector<Srvf> &Qs, 
                  double tol, size_t max_iters)
{
  // TODO
  Srvf res;
  return res;
}



double TestAccess::edge_variation(const Srvf &Q, size_t i1, size_t i2)
{
  return edge_variation_(Q, i1, i2);
}

double TestAccess::edge_score(const Srvf &Q1, const Srvf &Q2, 
  size_t sc, size_t sr, size_t tc, size_t tr)
{
  return edge_score_(Q1, Q2, sc, sr, tc, tr);
}

void TestAccess::calculate_scores(const Srvf &Q1, const Srvf &Q2, 
  std::map<match_vertex_t,double> &score, 
  std::map<match_vertex_t,match_vertex_t> &pred)
{
  calculate_scores_(Q1, Q2, score, pred);
}

void TestAccess::build_gamma_segment(const Srvf &Q1, const Srvf &Q2, 
  size_t sc, size_t sr, size_t tc, size_t tr, 
  std::vector<double> &G1samps, std::vector<double> &G2samps)
{
  build_gamma_segment_(Q1, Q2, sc, sr, tc, tr, G1samps, G2samps);
}

} // namespace srvf::functions

} // namespace srvf
