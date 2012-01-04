% Compute the SRVF for a piecewise-linear function.  Since SRVF's are 
% only defined for absolutely-continuous functions, T must be strictly 
% increasing.
%
% Inputs:
%  F - function values
%  T - corresponding parameter values
%  closed - nonzero iff closed curve (default 0).
%
% Outputs:
%  Q - SRVF values.  Q(i) = function value on interval [T(i), T(i+1)]
function Q = plf_to_srvf( F, T, closed )
  if ( nargin < 3 ) 
    closed = 0; 
  end
  assert( ismonotone(T,-1e-6), 1 );  % tol < 0 ==> need strictly increasing

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
