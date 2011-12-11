function [P T] = srvf_geodesic( Q1, T1, Q2, T2, Nsteps )
  T = union( T1, T2 );
  Q1r = srvf_refine( Q1, T1, T );
  Q2r = srvf_refine( Q2, T2, T );

  P = zeros( Nsteps, length(T)-1 );
  Ptv = linspace( 0, 1, Nsteps );

  ip = srvf_l2product( Q1r, T, Q2r, T );
  theta = acos( ip );

  for i=1:Nsteps
    P(i,:)=(sin(theta*(1-Ptv(i))).*Q1r + sin(theta*Ptv(i)).*Q2r) ./ sin(theta);
  end

end
