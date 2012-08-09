% Returns a new SRVF representing the restriction of Q,T to the interval 
% [a,b].  [a,b] should be a subset of [T(1),T(end)], the domain of the SRVF.
% --------------------------------------------------------------------------
function [Qr Tr] = srvf_restrict(Q,T,a,b)
  Tr = unique([T a b]);
  Qr = srvf_refine(Q,T,Tr);
  aidx = find(Tr==a);
  bidx = find(Tr==b);
  Qr = Qr(:,aidx:bidx-1);
  Tr = Tr(aidx:bidx);
end

%!test
%! T=[0 1/4 1/2 3/4 1];
%! Q=[0 1 0 1];
%! a=1/8;
%! b=5/8;
%! Trx=[1/8 1/4 1/2 5/8];
%! Qrx=[0 1 0];
%! [Qr Tr]=srvf_restrict(Q,T,a,b);
%! assert(Qr,Qrx);
%! assert(Tr,Trx);
