% Computes an SRVF which is a linear combination of two given SRVF's with 
% the specified weights.  The SRVFs must have the same dimension and must 
% be defined on the same interval.
%
% Inputs
%  Q1,T1 : the first SRVF
%  Q2,T2 : the second SRVF
%  w1,w2 : the weights
%
% Outputs
%  Q,T : an SRVF representing w1*Q1 + w2*Q2
% ---------------------------------------------------------
function [Q T] = srvf_linear_combination( Q1, T1, Q2, T2, w1, w2 )
  T = unique( [T1 T2] );

  Q1r = srvf_refine( Q1, T1, T );
  Q2r = srvf_refine( Q2, T2, T );

  Q = w1*Q1r + w2*Q2r;
end


%!test
%! Q1=[1 -1/2];
%! T1=[0 1/2 1];
%! Q2=[1 -1 1];
%! T2=[0 1/4 3/4 1];
%! w1=e;
%! w2=pi;
%! Qexp=[e+pi e-pi -e/2-pi -e/2+pi];
%! Texp=[0 1/4 1/2 3/4 1];
%! [Q T]=srvf_linear_combination(Q1,T1,Q2,T2,w1,w2);
%! assert(Q,Qexp,1e-5);
%! assert(T,Texp,1e-5);
