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
#ifndef SRVF_UTIL_H
#define SRVF_UTIL_H 1

#include "matrix.h"

namespace srvf
{

namespace util
{

Matrix linspace(double a, double b, int n);
Matrix unique(Matrix v1, Matrix v2);
Matrix diff(const Matrix &X);
Matrix diff(const Matrix &X, const Matrix &tv);

} // namespace srvf::util

} // namespace srvf

#endif // SRVF_UTIL_H
