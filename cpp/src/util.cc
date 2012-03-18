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

#include "util.h"

namespace srvf
{

Matrix linspace(double a, double b, int n)
{
  if (n<0) std::invalid_argument("n<0");

  Matrix R(1,n);
  if (n>1)
  {
    R(0)=a;
    double dt=(b-a)/(double)(n-1);
    for (int i=1; i<n; ++i)
    {
      R(i)=R(i-1)+dt;
    }
  }
  else if (n>0)
  {
    R(0)=a;
  }
  return R;
}

}
