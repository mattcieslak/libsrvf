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
#ifndef SRVF_PLF_H
#define SRVF_PLF_H 1

#include <stdexcept>
#include <vector>

#include "pointset.h"
#include "matrix.h"
#include "util.h"

namespace srvf {


/**
 * Represents a piecewise-linear function.
 */
class Plf
{
public:
  
  /** Default constructor. */
  Plf() { }

  /**
   * Creates a new \c Plf with the given sample points.
   *
   * Corresponding parameter values will be uniformly spaced from 0 to 1.
   *
   * \param samples
   */
  Plf(const Pointset &samples)
    : samps_(samples), params_()
  {
    params_=srvf::util::linspace(0.0,1.0,samples.npts());
  }

  /**
   * Creates a new \c Plf with the given sample points and parameters.
   *
   * \param samples
   * \param parameters
   */
  Plf(const Pointset &samples, const std::vector<double> &parameters)
    : samps_(samples), params_(parameters)
  {
    if (samples.npts()!=parameters.size())
      throw std::invalid_argument("number of samples != number of parameters");
  }

  /**
   * Copy constructor.
   *
   * Creates a deep copy of \a F
   * \param F an existing \c Plf
   */
  Plf(const Plf &F)
    : samps_(F.samps_), params_(F.params_)
  { }

  /**
   * Assignment operator.
   *
   * Sets the current \c Plf to a deep copy of \a F.
   * \param F an existing \c Plf
   * \return a reference to the current \c Plf
   */
  Plf &operator= (const Plf &F)
  {
    samps_ = F.samps_;    // deep copy
    params_ = F.params_;  // deep copy
    return *this;
  }

  /** Returns the sample points. */
  Pointset &samps() { return samps_; }

  /** Returns the sample points. */
  const Pointset &samps() const { return samps_; }

  /** Returns the changepoint parameters. */
  std::vector<double> &params() { return params_; }

  /** Returns the changepoint parameters. */
  const std::vector<double> &params() const { return params_; }

  /** Returns the dimension of the ambient space. */
  size_t dim() const { return samps_.dim(); }

  /** Returns the number of changepoints. */
  size_t ncp() const { return samps_.npts(); }

  /** Does this \c Plf represent the empty map? */
  bool is_empty() const { return (samps_.npts() == 0); }
 
  void evaluate(double t, Pointset &result) const;
  void evaluate(const std::vector<double> &tv, Pointset &result) const;
  void preimages(const std::vector<double> &tv, 
                 std::vector<double> &result) const;

  Pointset evaluate(double t) const;
  Pointset evaluate(const std::vector<double> &tv) const;
  std::vector<double> preimages(const std::vector<double> &tv) const;
  double arc_length() const;

  void translate(const std::vector<double> &v);
  void rotate(const Matrix &R);
  void scale(double s);

  friend Plf linear_combination(const Plf &F1, const Plf &F2, 
                                double w1, double w2);
  friend Plf composition(const Plf &F1, const Plf &F2);
  friend Plf inverse(const Plf &F);

private:
  Pointset samps_;
  std::vector<double> params_;
};

Plf linear_combination(const Plf &F1, const Plf &F2, double w1, double w2);
Plf composition(const Plf &F1, const Plf &F2);
Plf inverse(const Plf &F);

} // namespace srvf

#endif // SRVF_PLF_H
