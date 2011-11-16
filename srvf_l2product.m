% Compute the L^2 inner product of two SRVFs.
function ip = srvf_l2product( Q1, T1, Q2, T2, closed )
  if ( nargin < 5 )
    closed = 0;
  end

  % Is the sort() necessary?
  Tr = sort( union( T1, T2 ) );

  Q1r = srvf_refine( Q1, T1, Tr );
  Q2r = srvf_refine( Q2, T2, Tr );

  ip = sum( diff(Tr) .* sum(Q1r.* Q2r,1) );
end
