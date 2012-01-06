% Computes the square of the L^2 distance between two SRVF's.
% The SRVFs must be defined on the same interval.
%
% Inputs
%  Q1,T1 : the first SRVF
%  Q2,T2 : the second SRVF
%
% Outputs
%  d : the square of the L^2 distance between the SRVFs
% --------------------------------------------------------------------------
function d = srvf_squared_l2distance( Q1, T1, Q2, T2 )
  [Qd Td] = srvf_linear_combination( Q1, T1, Q2, T2, 1, -1 );
  d = srvf_squared_l2norm( Qd, Td );
end


%!test
%! Q1=[1 -1];
%! T1=[0 1/4 1];
%! Q2=[-1 1 2];
%! T2=[0 1/2 3/4 1];
%! assert(srvf_squared_l2distance(Q1,T1,Q2,T2),4.25,1e-4);
