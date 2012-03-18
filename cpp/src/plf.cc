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

#include <stdexcept>

#include "plf.h"
#include "matrix.h"
#include "util.h"

namespace srvf
{

/**
 * Creates a new \c Plf with the given sample points.
 *
 * Corresponding parameter values will be uniformly spaced from 0 to 1.
 *
 * \param samples
 */
Plf::Plf(const Matrix &samples)
  : samps_(samples), params_()
{
  params_=srvf::linspace(0.0,1.0,samples.cols());
}

Plf::Plf(const Matrix &samples, const Matrix &parameters)
  : samps_(samples), params_(parameters)
{
  if (params_.rows()!=1)
    std::invalid_argument("parameters must have 1 row");
  if (samps_.cols()!=params_.cols())
    std::invalid_argument("samples.cols()!=parameters.cols()");
}

void Plf::evaluate(double t, Matrix &result)
{

}

void Plf::evaluate(const Matrix &tv, Matrix &result)
{

}

void Plf::preimages(const Matrix &tv, Matrix &result)
{

}

double Plf::arc_length() const
{

}


void Plf::translate(const Matrix &v)
{

}

void Plf::rotate(const Matrix &R)
{

}

void Plf::scale(double s)
{

}

void Plf::linear_combine(const Plf &F, double w1, double w2)
{

}

void Plf::compose(const Plf &F)
{

}

void Plf::invert()
{

}


Srvf Plf::to_srvf() const
{

}


Plf translation(const Plf &F, const Matrix &v)
{

}

Plf rotation(const Plf &F, const Matrix &R)
{

}

Plf scaling(const Plf &F, double s)
{

}

Plf linear_combination(const Plf &F1, const Plf &F2, 
             double w1, double w2)
{

}

Plf composition(const Plf &F1, const Plf &F2)
{

}

Plf inverse(const Plf &F)
{

}


}
