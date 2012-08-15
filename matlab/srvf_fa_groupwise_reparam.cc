/*
 * libsrvf
 * =======
 *
 * A shape analysis library using the square root velocity framework.
 *
 * Copyright (C) 2012   FSU Statistical Shape Analysis and Modeling Group
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
#include <srvf/srvf.h>
#include <srvf/functions.h>

#include <mex.h>
#include <vector>


void do_usage()
{
  mexPrintf(
    "USAGE: [Gm,TGm,Gs,TGs] = %s(Qm,Tm,Qs,Ts)\n"
    "Inputs:\n"
    "\tQm, Tm : sample points and parameters of the mean SRVF\n"
    "\tQs, Ts : sample points and parameters of the other SRVFs\n"
    "Outputs:\n"
    "\tGm,TGm = the reparametrization for Qm\n",
    "\tGs,TGs = the reparametrizations for the Qs\n",
    mexFunctionName()
  );
}


extern "C"
{
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  size_t nfuncs;

  std::vector<srvf::Srvf> Qs;
  srvf::Srvf Qm;
  mxArray *sampsi_data;
  mxArray *paramsi_data;
  srvf::Pointset sampsi;
  std::vector<double> paramsi;

  std::vector<srvf::Plf> Gs;
  mxArray *Gi_matrix, *Ti_matrix;
  double *Gi_data, *Ti_data;
  mwSize Gs_dim;


  // Check arguments
  if (nrhs != 4 || nlhs != 4)
  {
    mexPrintf("error: incorrect number of arguments.\n");
    do_usage();
    return;
  }
  if (!mxIsCell(prhs[2]) || !mxIsCell(prhs[3]))
  {
    mexPrintf("error: Qs and Ts must be cell arrays.\n");
    do_usage();
    return;
  }
  nfuncs = mxGetNumberOfElements(prhs[2]);
  if (mxGetNumberOfElements(prhs[3]) != nfuncs)
  {
    mexPrintf("error: Qs and Ts must have the same number of elements.\n");
    do_usage();
    return;
  }


  // Create the SRVFs
  for (size_t i=0; i<nfuncs; ++i)
  {
    sampsi_data = mxGetCell(prhs[2],i);
    paramsi_data = mxGetCell(prhs[3],i);

    // More input checking
    if (mxGetM(sampsi_data) != 1 || mxGetM(paramsi_data) != 1)
    {
      mexPrintf("error: elements of Qs and Ts must have 1 row.\n");
      return;
    }
    if (mxGetN(paramsi_data) != mxGetN(sampsi_data)+1)
    {
      mexPrintf("error: Ts(%d) must have length size(Qs(%d),2)+1.\n", i, i);
      return;
    }

    sampsi = srvf::Pointset(1, mxGetN(sampsi_data), mxGetPr(sampsi_data));
    paramsi = std::vector<double>(mxGetPr(paramsi_data), 
      mxGetPr(paramsi_data)+mxGetN(paramsi_data));

    Qs.push_back(srvf::Srvf(sampsi, paramsi));
  }
  if (mxGetM(prhs[0]) != 1 || mxGetM(prhs[1]))
  {
    mexPrintf("error: Qm and Tm must have 1 row\n");
    return;
  }
  if (mxGetN(prhs[0])+1 != mxGetN(prhs[1]))
  {
    mexPrintf("error: Tm must have length size(Qm,2)+1.\n");
    return;
  }
  sampsi = srvf::Pointset(1,mxGetN(prhs[0]),mxGetPr(prhs[0]));
  paramsi = std::vector<double> (mxGetPr(prhs[1]), 
    mxGetPr(prhs[1])+mxGetN(prhs[1]));
  Qm = srvf::Srvf(sampsi,paramsi);


  // All SRVFs must be unit-norm, constant-speed SRVFs
  for (size_t i=0; i<nfuncs; ++i)
  {
    Qs[i].scale_to_unit_norm();
    Qs[i] = srvf::constant_speed_param(Qs[i]);
  }
  Qm.scale_to_unit_norm();
  Qm = srvf::constant_speed_param(Qm);


  // Compute the groupwise alignment using libsrvf
  Gs = srvf::functions::groupwise_optimal_reparam(Qm,Qs);


  // Allocate output variables
  Gs_dim = (mwSize)nfuncs;
  plhs[2] = mxCreateCellArray(1,&Gs_dim);
  plhs[3] = mxCreateCellArray(1,&Gs_dim);
  if (!plhs[2] || !plhs[3])
  {
    mexPrintf("error: mxCreateCellArray() failed.\n");
    if (plhs[2]) mxDestroyArray(plhs[2]);
    if (plhs[3]) mxDestroyArray(plhs[3]);
    return;
  }

  plhs[0] = mxCreateDoubleMatrix(1, Gs.back().ncp(), mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, Gs.back().ncp(), mxREAL);
  if (!plhs[0] || !plhs[1])
  {
    mexPrintf("error: mxCreateDoubleMatrix() failed.\n");
    mxDestroyArray(plhs[2]);
    mxDestroyArray(plhs[3]);
    if (plhs[0]) mxDestroyArray(plhs[0]);
    if (plhs[1]) mxDestroyArray(plhs[1]);
    return;
  }
  for (size_t i=0; i<nfuncs; ++i)
  {
    Gi_matrix = mxCreateDoubleMatrix(1, Gs[i].ncp(), mxREAL);
    Ti_matrix = mxCreateDoubleMatrix(1, Gs[i].ncp(), mxREAL);
    if (!Gi_matrix || !Ti_matrix)
    {
      if (Gi_matrix) mxDestroyArray(Gi_matrix);
      if (Ti_matrix) mxDestroyArray(Ti_matrix);
      for (size_t j=0; j<i; ++j)
      {
        mxDestroyArray(mxGetCell(plhs[2],j));
        mxDestroyArray(mxGetCell(plhs[3],j));
      }
      mxDestroyArray(plhs[0]);
      mxDestroyArray(plhs[1]);
      mxDestroyArray(plhs[2]);
      mxDestroyArray(plhs[3]);

      mexPrintf("error: mxCreateDoubleMatrix() failed.\n");
      return;
    }
    mxSetCell(plhs[2], i, Gi_matrix);
    mxSetCell(plhs[3], i, Ti_matrix);
  }

  
  // Copy into plhs
  Gi_data = mxGetPr(plhs[0]);
  Ti_data = mxGetPr(plhs[1]);
  for (size_t i=0; i<Gs.back().ncp(); ++i)
  {
    Gi_data[i] = Gs.back().samps()[i][0];
    Ti_data[i] = Gs.back().params()[i];
  }
  for (size_t i=0; i<nfuncs; ++i)
  {
    Gi_data = mxGetPr(mxGetCell(plhs[2],i));
    Ti_data = mxGetPr(mxGetCell(plhs[3],i));

    for (size_t j=0; j<Gs[i].ncp(); ++j)
    {
      Gi_data[j] = Gs[i].samps()[j][0];
      Ti_data[j] = Gs[i].params()[j];
    }
  }

}
} // extern "C"
