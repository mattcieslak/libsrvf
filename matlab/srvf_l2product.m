% Computes the L^2 inner product of two SRVFs.  The SRVFs must have the 
% same dimension and must be defined on the same interval.
%
% Inputs
%  Q1,T1 : the first SRVF
%  Q2,T2 : the second SRVF
%
% Outputs
%  ip : the L^2 inner product of the two SRVFs
% ---------------------------------------------------------
function ip = srvf_l2product( Q1, T1, Q2, T2 )
  assert( size(Q1,1) == size(Q2,1) );
  assert( size(T1,1) == 1 );
  assert( size(T2,1) == 1 );
  assert( size(Q1,2) == size(T1,2)-1 );
  assert( size(Q2,2) == size(T2,2)-1 );

  Tr = unique( [T1 T2] );

  Q1r = srvf_refine( Q1, T1, Tr );
  Q2r = srvf_refine( Q2, T2, Tr );

  ip = sum( diff(Tr) .* sum(Q1r.* Q2r,1) );
end


%!test
%! Q1=[0 0 0 0];
%! T1=[0 1/8 1/2 3/4 1];
%! Q2=[1];
%! T2=[0 1];
%! assert( srvf_l2product(Q1,T1,Q2,T2), 0, 1e-6 );
%!
%!test
%! Q1=[1 -1 1 -1 1];
%! T1=linspace(0,1,6);
%! Q2=[1];
%! T2=[0 1];
%! assert( srvf_l2product(Q1,T1,Q2,T2), 0.2, 1e-6 );
