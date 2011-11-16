function d = srvf_squared_l2distance( Q1, T1, Q2, T2 )

  [Qd Td] = srvf_linear_combination( Q1, T1, Q2, T2, 1, -1 );
  d = srvf_squared_l2norm( Qd, Td );
end
