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
#include <cmath>
#include "qmap.h"
#include "util.h"
#include "plf.h"
#include "srvf.h"

namespace srvf
{

/**
 * Returns a new \c Srvf representing the square-root velocity function of 
 * the given \c Plf.
 *
 * \param F a \c Plf
 * \return a new \c Srvf representing the SRVF of \a F
 */
Srvf plf_to_srvf(const Plf &F)
{
  Matrix dF=srvf::util::diff(F.samps(),F.params());
  for (int i=0; i<dF.cols(); ++i)
  {
    double nqi=0.0;
    for (int j=0; j<dF.rows(); ++j)
    {
      nqi+=dF(j,i)*dF(j,i);
    }
    nqi=sqrt(sqrt(nqi));
    if (nqi>1e-6)
    {
      for (int j=0; j<dF.rows(); ++j)
      {
        dF(j,i)/=nqi;
      }
    }
    else
    {
      for (int j=0; j<dF.rows(); ++j)
      {
        dF(j,i)=0.0;
      }
    }
  }

  return Srvf(dF,F.params());
}


Plf srvf_to_plf(const Srvf &Q)
{
  Plf dummy;
  return dummy;
}

} // namespace srvf
