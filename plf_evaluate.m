% Evaluate a PLF at the given parameter values.  The given PLF will be 
% left-continuous at jump discontinuities, if any.
function f = plf_evaluate( F, T, t )
  if isscalar(t)
    % This 'feature' has unintended consequences and should be removed.
    %t = linspace(0,1,t);
  end

  f = zeros(rows(F),length(t));
  for i=1:rows(F)
    f(i,:) = interp1(T,F(i,:),t,'extrap');
  end

  % left-continuous at jump discontinuities
  jumps = find( diff(T) == 0 );
  for j=1:length(jumps)
    jump_tidx = find( t == T(jumps(j)) );
    for i=1:rows(F)
      f(i,jump_tidx) = F(i,jumps(j));
    end
  end
end


%!test
%! F=[0 1 2 3 4 5 6;11 10 9 8 7 6 5];
%! T=[0 0 1/4 1/4 1/2 1 1];
%! tv=[0 1/8 1/4 3/8 1/2 3/4 1];
%! X=[0.0000   1.5000   2.0000   3.5000   4.0000   4.5000   5.0000;\
%!    11.0000    9.5000    9.0000    7.5000    7.0000    6.5000    6.0000];
%! A=plf_evaluate(F,T,tv);
%! assert(A,X);
