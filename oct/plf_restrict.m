% Returns a new PLF representing the restriction of F,T to [a,b].
% [a,b] should be contained in [T(1),T(end)] (the original domain of the PLF).
% --------------------------------------------------------------------------
function [Fr Tr] = plf_restrict(F,T,a,b)
  Tr = unique([T a b]);
  Fr = plf_refine(F,T,Tr);
  aidx = find(Tr==a);
  bidx = find(Tr==b);
  Tr = Tr(aidx:bidx);
  Fr = Fr(:,aidx:bidx);
end

%!test
%! F=[0 1 0 1 0; 1 2 1 2 1];
%! T=[0 1/4 1/2 3/4 1];
%! a=1/8;
%! b=7/8;
%! Frx=[1/2 1 0 1 1/2; 3/2 2 1 2 3/2];
%! Trx=[1/8 1/4 1/2 3/4 7/8];
%! [Fr Tr]=plf_restrict(F,T,a,b);
%! assert(Tr,Trx);
%! assert(Fr,Frx);
