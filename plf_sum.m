% Returns the sum of two PLFs.
function [F T] = plf_sum( F1, T1, F2, T2 )
  T = union( T1, T2 );
  F = plf_evaluate( T1, F1, T ) + plf_eval( T2, F2, T );
end
