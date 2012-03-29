% Returns the sum of two PLFs.
function [F T] = plf_sum( F1, T1, F2, T2 )
  % TODO: why is this union() and not unique()?
  T = union( T1, T2 );
  F = plf_evaluate( F1, T1, T ) + plf_evaluate( F2, T2, T );
end


%!test
%! F1 = [0 1 2 3 4];
%! T1 = [0 1/4 1/2 3/4 1];
%! F2 = [5 6 7 8 9];
%! T2 = [0 1/8 3/8 5/8 1];
%! Texp = [0 1/8 1/4 3/8 1/2 5/8 3/4 1];
%! Fexp = [5.0 6.5 7.5 8.5 9.5 10.5 11.3333 13.0];
%! [F T] = plf_sum(F1,T1,F2,T2);
%! assert( F, Fexp, 1e-4 );
%! tv = linspace(0,1,100);
%! F1tv = interp1(T1,F1,tv);
%! F2tv = interp1(T2,F2,tv);
%! Ftv = interp1(T,F,tv);
%! assert( Ftv, F1tv+F2tv, 1e-4 );
