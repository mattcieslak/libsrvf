% Evaluates an SRVF at the given parameter values.
% t must be non-decreasing.
function q = srvf_evaluate( Q, T, t )
  assert( size(T,1) == 1 );
  assert( size(Q,2) == size(T,2)-1 );
  assert( min(diff(T)) >= 0);
  assert( min(diff(t)) >= 0 );
  epsval = 0.01*(T(end)-T(1));
  assert( t(1) > T(1)-epsval && t(end) < T(end)+epsval );

  q = zeros(size(Q,1),length(t));
  Qidx = 1;

  for qidx = 1:length(t)
    while ( Qidx < length(T)-1 && t(qidx) > T(Qidx+1) ) 
      Qidx = Qidx + 1; 
    end

    q(:,qidx) = Q(:,Qidx);
  end
end


%!test
%! Q=[0 1/2 1 -1];
%! T=[0 1/4 1/2 3/4 1];
%! t=[0 0.249 0.25 0.251 0.499 0.5 0.501 0.749 0.75 0.751 0.9 1];
%! qexp=[0 0 0 1/2 1/2 1/2 1 1 1 -1 -1 -1];
%! q=srvf_evaluate(Q,T,t);
%! assert(q,qexp,1e-3);
%!
%!#error
%! Q=[0 1/2 1 -1];
%! T=[0 1/4 1/2 3/4 1];
%! t=[-eps];
%! q=srvf_evaluate(Q,T,t);
%!
%!#error
%! Q=[0 1/2 1 -1];
%! T=[0 1/4 1/2 3/4 1];
%! t=[1+eps];
%! q=srvf_evaluate(Q,T,t);
%!
%!#error
%! Q=[0 1/2 1 -1];
%! T=[0 1/4 1/2 3/4 1];
%! t=[0 1/2 1/4 1];
%! q=srvf_evaluate(Q,T,t);
