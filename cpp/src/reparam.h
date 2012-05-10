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
#ifndef SRVF_REPARAM_H
#define SRVF_REPARAM_H 1

#include <vector>

#include "plf.h"
#include "srvf.h"

namespace srvf
{

namespace opencurves
{

Plf optimal_reparam(const Srvf &Q1, const Srvf &Q2);

//TODO: support arbitrary grids
//Plf optimal_reparam (const Srvf &Q1, const Srvf &Q2
//                 std::vector<double> &gv1, std::vector<double> &gv2);


/**
 * Computes a partial matching cost.
 *
 * Computes the cost of matching part of \a Q1 to part of \a Q2.  
 * Let \c a=Q1.params()[Q1i1], \c b=Q1.params()[Q1i2], 
 * \c c=Q2.params()[Q2i1], and \c d=Q2.params()[Q2i2].
 * Then the matching cost is
 *
 * \f[ \int_a^b|Q_1(t)-\sqrt{\dot{\eta}(t)}Q_2(\eta(t))|^2 dt \f]
 *
 * where \f$ \eta:[a,b]\to[c,d] \f$ is the linear map with 
 * \f$ \eta(a)=c \f$ and \f$ \eta(b)=d \f$.
 *
 * \param Q1 
 * \param Q1i1
 * \param Q1i2
 * \param Q2
 * \param Q2i1
 * \param Q2i2
 * \return the matching cost
 */
inline double
match_cost (const Srvf &Q1, size_t Q1i1, size_t Q1i2, 
            const Srvf &Q2, size_t Q2i1, size_t Q2i2)
{
  double slope = (double)(Q2.params()[Q2i2]-Q2.params()[Q2i1]) / 
                 (double)(Q1.params()[Q1i2]-Q1.params()[Q1i1]);
  double rslope = sqrt(slope);
  size_t dim = Q1.dim();
  double res = 0.0;

  size_t i1 = Q1i1;
  size_t i2 = Q2i1;
  double t1prev = Q1.params()[i1];
  double t2prev = Q2.params()[i2];
  
  while (i1<Q1i2 && i2<Q2i2)
  {
    double dt1 = Q1.params()[i1+1] - t1prev;
    double dt2 = Q2.params()[i2+1] - t2prev;

    double t1next_cand1 = Q1.params()[i1+1];
    double t1next_cand2 = t1prev + dt2 / slope;
    double t1next, t2next;
    double i1next, i2next;
    
    if (t1next_cand1 < t1next_cand2-1e-6)
    {
      t1next = t1next_cand1;
      t2next = t2prev + dt1 * slope;
      i1next = i1+1;
      i2next = i2;
    }
    else if (t1next_cand2 < t1next_cand1-1e-6)
    {
      t1next = t1next_cand2;
      t2next = Q2.params()[i2+1];
      i1next = i1;
      i2next = i2+1;
    }
    else
    {
      t1next = Q1.params()[i1+1];
      t2next = Q2.params()[i2+1];
      i1next = i1+1;
      i2next = i2+1;
    }

    double s = 0.0;
    for (size_t j=0; j<dim; ++j)
    {
      double dqi = Q1.samps()(i1,j) - rslope * Q2.samps()(i2,j);
      s += dqi * dqi;
    }
    double dt = t1next-t1prev;
    res += s * dt;

    t1prev = t1next;
    t2prev = t2next;
    i1 = i1next;
    i2 = i2next;
  }

  return res;
}

} // namespace opencurves

} // namespace srvf

#endif // SRVF_REPARAM_H
