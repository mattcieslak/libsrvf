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
    "USAGE: [G1,T1,G2,T2] = %s(Q1,T1,Q2,T2)\n"
    "Inputs:\n"
    "\tQ1, T1 : sample points and parameters of the first SRVF\n"
    "\tQ2, T2 : sample points and parameters of the second SRVF\n"
    "Outputs:\n"
    "\tG1,T1 : the reparametrization for Q1\n"
    "\tG2,T2 : the reparametrization for Q2\n",
    mexFunctionName()
  );
}


extern "C"
{
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  srvf::Srvf Q1, Q2;
  srvf::Pointset Q1samps, Q2samps;
  std::vector<double> Q1params, Q2params;
  size_t dim;
  size_t ncp1, ncp2;

  std::vector<srvf::Plf> Gs;
  double *G1_data, *T1_data;
  double *G2_data, *T2_data;


  // Check arguments
  if (nrhs != 4 || nlhs != 4)
  {
    mexPrintf("error: incorrect number of arguments.\n");
    do_usage();
    return;
  }


  // Get dimension and number of points
  dim = mxGetM(prhs[0]);
  if (mxGetM(prhs[2]) != dim)
  {
    mexPrintf("error: Q1 and Q2 have different dimensions.\n");
    do_usage();
    return;
  }
  ncp1 = mxGetN(prhs[1]);
  ncp2 = mxGetN(prhs[3]);
  if (mxGetN(prhs[0])+1 != ncp1)
  {
    mexPrintf("error: bad length on T1 (must be size(Q1,2) + 1).\n");
    do_usage();
    return;
  }
  if (mxGetN(prhs[2])+1 != ncp2)
  {
    mexPrintf("error: bad length on T2 (must be size(Q2,2) + 1).\n");
    do_usage();
    return;
  }


  // Create the SRVFs
  Q1samps = srvf::Pointset(dim, ncp1-1, mxGetPr(prhs[0]));
  Q1params = std::vector<double> (mxGetPr(prhs[1]), mxGetPr(prhs[1])+ncp1);
  Q1 = srvf::Srvf(Q1samps, Q1params);

  Q2samps = srvf::Pointset(dim, ncp2-1, mxGetPr(prhs[2]));
  Q2params = std::vector<double> (mxGetPr(prhs[3]), mxGetPr(prhs[3])+ncp2);
  Q2 = srvf::Srvf(Q2samps, Q2params);


  // Q1 and Q2 need to be unit-norm, constant-speed SRVFs
  Q1.scale_to_unit_norm();
  Q1 = srvf::constant_speed_param(Q1);
  Q2.scale_to_unit_norm();
  Q2 = srvf::constant_speed_param(Q2);


  // Get the reparametrizations for Q1 and Q2 using libsrvf
  Gs = srvf::functions::optimal_reparam(Q1, Q2);


  // Allocate output variables and copy Gs[0] and Gs[1] into plhs
  plhs[0] = mxCreateDoubleMatrix(1,Gs[0].ncp(), mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,Gs[0].ncp(), mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1,Gs[1].ncp(), mxREAL);
  plhs[3] = mxCreateDoubleMatrix(1,Gs[1].ncp(), mxREAL);
  if (!plhs[0] || !plhs[1] || !plhs[2] || !plhs[3])
  {
    if (plhs[0]) mxDestroyArray(plhs[0]);
    if (plhs[1]) mxDestroyArray(plhs[1]);
    if (plhs[2]) mxDestroyArray(plhs[2]);
    if (plhs[3]) mxDestroyArray(plhs[3]);

    mexPrintf("error: mxCreateDoubleMatrix() failed.\n");
    return;
  }

  G1_data = mxGetPr(plhs[0]); T1_data = mxGetPr(plhs[1]);
  G2_data = mxGetPr(plhs[2]); T2_data = mxGetPr(plhs[3]);

  for (size_t i=0; i<Gs[0].ncp(); ++i)
  {
    G1_data[i] = Gs[0].samps()[i][0];
    T1_data[i] = Gs[0].params()[i];
  }
  for (size_t i=0; i<Gs[1].ncp(); ++i)
  {
    G2_data[i] = Gs[1].samps()[i][0];
    T2_data[i] = Gs[1].params()[i];
  }

}
} // extern "C"
