% Computes the great circle distance between two SRVFs.
function d = srvf_preshape_distance( Q1, T1, Q2, T2 )
  ip = srvf_l2product( Q1, T1, Q2, T2 );
  n1 = sqrt( srvf_squared_l2norm( Q1, T1 ) );
  n2 = sqrt( srvf_squared_l2norm( Q2, T2 ) );

  % Preshape distance isn't defined if SRVF's have different norms.
  assert( abs(n1-n2) < max(0.01*max(n1,n2), 1e-3) );
  
  if ( n1*n2 > 1e-4 )
    d = acos( ip / n1 / n2 );
  else
    d = 0;
  end

end
