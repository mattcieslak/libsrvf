function [Q T] = srvf_linear_combination( Q1, T1, Q2, T2, w1, w2 )
  % Is the sort() necessary?
  T = sort( union( T1, T2 ) );

  Q1r = srvf_refine( Q1, T1, T );
  Q2r = srvf_refine( Q2, T2, T );

  Q = w1*Q1r + w2*Q2r;
end
