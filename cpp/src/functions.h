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
#ifndef SRVF_FUNCTIONS_H
#define SRVF_FUNCTIONS_H 1

#include "srvf.h"
#include "plf.h"

#include <vector>


namespace srvf
{

namespace functions
{

/**
 * Returns reparametrizations which optimally register \a Q1 and \a Q2.
 *
 * Unlike \c srvf::opencurves::optimal_reparam(), which uses the older 
 * dynamic programming algorithm, this function implements the new 
 * specialized 1-D function matching.
 *
 * \param Q1 the SRVF of the first function
 * \param Q2 the SRVF of the second function
 * \return a \c vector containing two \c Plf's representing non-decreasing 
 *   functions \f$ \gamma_i:[0,1]\to[0,1] \f$.  When the first is applied 
 *   to \c Q1 and the second is applied to \c Q2 (using \c gamma_action()), 
 *   the result is a pair of representatives of the closed-up orbits of 
 *   \c Q1 and \c Q2 which are at minimum distance.
 */
std::vector<Plf> optimal_reparam(const Srvf &Q1, const Srvf &Q2);


/**
 * Computes the Karcher mean of a collecton of 1-D SRVFs.
 *
 * Unlike \c srvf::opencurves::karcher_mean(), this function uses the 
 * new 1-D function matching algorithm.
 */
Srvf karcher_mean(const std::vector<Srvf> &Qs, 
                  double tol=1e-3, size_t max_iters=0);

  
} // namespace srvf::functions

} // namespace srvf

#endif // SRVF_FUNCTIONS_H
