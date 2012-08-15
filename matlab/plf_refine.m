% Returns a PLF which is equivalent to F,T, but has changepoint parameters 
% Tr.  Tr should be a refinement of T (i.e. T(1)=Tr(1), T(end)=Tr(end), 
% and T is a subset of Tr).
% --------------------------------------------------------------------------
function Fr = plf_refine(F,T,Tr)
  Fr = plf_evaluate(F,T,Tr);
end

%!test
%! F=[0 1 0 1 0];
%! T=[0 1/4 1/2 3/4 1];
%! Tr=[0 1/8 1/4 3/8 1/2 5/8 3/4 7/8 1];
%! Fr=plf_refine(F,T,Tr);
%! assert(Fr,[0 1/2 1 1/2 0 1/2 1 1/2 0]);
