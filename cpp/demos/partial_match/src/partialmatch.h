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
#ifndef PARTIAL_MATCH_H
#define PARTIAL_MATCH_H 1

#include "paretoset.h"
#include <srvf.h>

/**
 * Returns the set of Pareto-optimal partial matches between two curves.
 *
 * \param Q1 the SRVF of the first curve
 * \param Q2 the SRVF of the second curve
 * \param do_rots Optimize over rotations?  Default is \c true.
 * \param grid_width partial matching grid width.  If unspecified, the 
 *        matching grid will have a column for every changepoint parameter 
 *        of \a Q1.
 * \param grid_height partial matching grid height.  If unspecified, the 
 *        matching grid will have a column for every changepoint parameter 
 *        of \a Q2.
 * \param nbuckets number of buckets to use in the Pareto set.  Each 
 *        bucket corresponds to a range of match lengths.  If unspecified, 
 *        the Pareto set will have a bucket for every possible match length.
 * \param bucket_thresh matches having a shape distance within 
 *        \a bucket_thresh of the minimum shape distance for their length 
 *        range are considered Pareto optimal.
 */
ParetoSet find_matches (
  const srvf::Srvf &Q1, const srvf::Srvf &Q2, 
  bool do_rots=true, size_t grid_width=0, size_t grid_height=0, 
  size_t nbuckets=0, double bucket_thresh=0.01 );


#endif // PARTIAL_MATCH_H
