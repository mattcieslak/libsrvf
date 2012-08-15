% Returns the arclength parameter values for a polygon P.
%
% Inputs:
%  P - the polygon.  An n x k matrix, where n is the dimension of the 
%      ambient space, and k is the number of vertices.
%
% Outputs:
%  clp - clp(i) = sum of lengths of the first (i-1) segments of P
% --------------------------------------------------------------------------
function clp = poly_arclength_params(P)
  dP = diff(P,1,2);
  clp = [0 cumsum(sqrt(sum(dP.*dP,1)))];
end

%!test
%! P=[0 1 1 0; 0 0 1 1];
%! clp=poly_arclength_params(P);
%! assert(clp,[0 1 2 3]);
