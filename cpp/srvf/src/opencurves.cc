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
#include "opencurves.h"
#include "plf.h"
#include "rotate.h"
#include "reparam.h"

#include <cmath>
#include <iostream>


namespace srvf
{

namespace opencurves
{

/**
 * Computes the Karcher mean of a collection of points on the unit 
 * sphere in \f$ L^2(I,R^n) \f$.
 */
Srvf karcher_mean(const std::vector<Srvf> &Qs, double tol, size_t max_iters, 
                  bool optimize_rots, bool optimize_reparams)
{
  // Corner case:  if Qs is empty, return an empty Srvf
  if (Qs.empty())
  {
    Srvf empty;
    return empty;
  }

  double radius = srvf::l2_norm(Qs[0]);
  double stepsize = 0.25;

  // Initialize mean to one of the Qs.  This is arbitrary, and not necessarily 
  // the best solution.  But (hopefully) good enough.
  Srvf Mu(Qs[0]);
  
  for (size_t iter=1; (max_iters==0 || iter<=max_iters); ++iter)
  {
    // The update direction is the average of the shooting directions.
    Srvf update_dir(Mu.domain_lb(), Mu.domain_ub(), 
                    std::vector<double>(Mu.dim(),0.0));

#pragma omp parallel
{
    #pragma omp for schedule (static) nowait
    for (size_t i=0; i<Qs.size(); ++i)
    {
      Srvf Svi;

      if (optimize_rots && optimize_reparams)
      {
        Srvf Qsi(Qs[i]);  // mutable local copy

        // Rotate
        Matrix Ri = optimal_rotation(Mu, Qsi);
        Qsi.rotate(Ri);

        // Reparametrize
        Plf Gi = optimal_reparam(Mu, Qsi);
        Srvf Qsir = gamma_action(Qsi, Gi);

        // And rotate again
        Ri = optimal_rotation(Mu, Qsir);
        Qsir.rotate(Ri);

        // Compute current shooting vector
        Svi = shooting_vector(Mu, Qsir);
      }
      else if (optimize_rots)
      {
        Srvf Qsi(Qs[i]);  // mutable local copy

        Matrix Ri = optimal_rotation(Mu, Qsi);
        Qsi.rotate(Ri);
        Svi = shooting_vector(Mu, Qsi);
      } 
      else if (optimize_reparams)
      {
        Srvf Qsi(Qs[i]);  // mutable local copy

        Plf Gi = optimal_reparam(Mu, Qsi);
        Srvf Qsir = gamma_action(Qsi, Gi);
        Svi = shooting_vector(Mu, Qsir);
      }
      else
      {
        Svi = shooting_vector(Mu, Qs[i]);
      }

      #pragma omp critical
      update_dir = srvf::linear_combination(update_dir, Svi, 1.0, 1.0);
    }
}

    update_dir.scale(1.0 / (double)Qs.size());
    double udr = l2_norm(update_dir);
    std::cout << iter << ") karcher_mean(): gradient norm = " << udr 
              << std::endl;
    if (udr < tol)
    {
      // Norm of gradient is small enough to stop here
      break;
    }
    else
    {
      // Update mean estimate and keep going
      Mu = srvf::linear_combination(Mu, update_dir, 1.0, stepsize);
      Mu.scale(radius / l2_norm(Mu));
    }
  }

  return Mu;
}


/**
 * Returns the shooting vector from \a Q1 to \a Q2.
 *
 * \a Q1 and \a Q2 must have the same \f$ L^2 \f$ norm.
 * 
 * The shooting vector from \a Q1 to \a Q2 is the tangent vector at \a Q1 
 * of a spherical geodesic from \a Q1 to \a Q2.  Its length is equal to the 
 * geodesic distance between the two points.
 */
Srvf shooting_vector(const Srvf &Q1, const Srvf &Q2)
{
  double radius = srvf::l2_norm(Q1);
  double radius2 = radius * radius;

  // Check for case where Q1 and Q2 are close to zero
  if (radius2 < 1e-8)
  {
    return Srvf(Q1);
  }

  // Angle between Q1 and Q2
  double cos_theta = srvf::l2_product(Q1, Q2) / radius2;
  if (cos_theta >  1.0) cos_theta =  1.0;
  if (cos_theta < -1.0) cos_theta = -1.0;
  double theta = acos(cos_theta);

  // Shooting vector
  Srvf sv = srvf::linear_combination(Q1, Q2, -cos_theta, 1.0);

  // Scale shooting vector to have length equal to great-circle distance
  double norm_sv = srvf::l2_norm(sv);
  if (norm_sv > 1e-6)
  {
    sv.scale(radius * theta / norm_sv);
  }

  return sv;
}

} //namespace opencurves
} // namespace srvf
