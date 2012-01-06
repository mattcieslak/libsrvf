% Compute the SRVF for a piecewise-linear function.  Since SRVF's are 
% only defined for absolutely-continuous functions, T must be strictly 
% increasing.
%
% Inputs:
%  F - function values
%  T - corresponding parameter values
%
% Outputs:
%  Q - SRVF values.  Q(i) = function value on interval [T(i), T(i+1)]
function Q = plf_to_srvf( F, T )
  if ( nargin < 3 ) 
    closed = 0; 
  end
  assert( min(diff(T)) > 0 );

  [dim, ncp] = size(F);

  if ( dim > 1 )
    V = diff(F,1,2);

    for i=1:dim
      V(i,:) = V(i,:) ./ diff( T );
    end

    Vrmag = sqrt( sqrt( sum( V .* V, 1 ) ) );
    zidx = find( Vrmag < 1e-4 );

    Q = zeros( dim, ncp-1 );
    for i=1:dim
      Q(i,:) = V(i,:) ./ Vrmag;
      Q(i,zidx) = 0;
    end

  else
    m = diff(F) ./ diff(T);
    Q = sign(m) .* sqrt(abs(m));
  end
end

%!test
%! F=[0 1 1 0; 
%!    0 0 1 1];
%! T=[0 1/3 2/3 1];
%! Qexp=[sqrt(3)  0.00000 -sqrt(3);
%!       0.00000  sqrt(3)  0.00000];
%! Q=plf_to_srvf(F,T);
%! assert(Q,Qexp,1e-3);
%!
%!test
%! F=[0 1 0 -1;
%!    0 1 2  1];
%! T=[0 1 2 3];
%! x=2**(-1/4);
%! Qexp=[x -x -x;
%!       x  x -x];
%! Q=plf_to_srvf(F,T);
%! assert(Q,Qexp,1e-3);
