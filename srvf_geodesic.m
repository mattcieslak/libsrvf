% Computes several steps along the minimum-length great circle path between 
% the given SRVF's.
function [P T] = srvf_geodesic( Q1, T1, Q2, T2, Nsteps )
  [dim nsegs] = size(Q1);

  T = unique( [T1 T2] );
  Q1r = srvf_refine( Q1, T1, T );
  Q2r = srvf_refine( Q2, T2, T );

  P = zeros( dim, nsegs, Nsteps );
  Ptv = linspace( 0, 1, Nsteps );

  ip = srvf_l2product( Q1r, T, Q2r, T );
  theta = acos( ip );

  for i=1:Nsteps
    P(:,:,i)=(sin(theta*(1-Ptv(i)))*Q1r + sin(theta*Ptv(i))*Q2r) ./ sin(theta);
  end

end
