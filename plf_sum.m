function [F T] = plf_sum( F1, T1, F2, T2 )
  T = union( T1, T2 );
  F = interp1( T1, F1, T ) + interp1( T2, F2, T );
end
