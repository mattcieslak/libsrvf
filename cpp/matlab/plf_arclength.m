% For closed curves, F(:,1)==F(:,end) (i.e. the seedpoint is duplicated), 
% so computing the arclength is the same as for open curves.
function [l S] = plf_arclength( F, T )
  [dim npts] = size( F );

  if ( dim > 1 )
    dF = diff(F,1,2);

    S = [0 cumsum( sqrt(sum( dF .* dF, 1 )), 2 )];
    l = S(end);
  else
    S = [0 cumsum( abs( diff(F) ) )];
    l = S(end);
  end
end

%!test
%! F=[0 1/4 0 1/4 0];
%! T=linspace(0,1,5);
%! assert(plf_arclength(F,T),1,1e-5);
%!
%!test
%! F=[4 -1 -3 5 7 1 0; 3 7 -3 -7 2 1 0; 1 2 3 4 5 6 7];
%! T=0:6;
%! assert(plf_arclength(F,T),42.898,1e-3);

