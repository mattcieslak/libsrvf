% Compute the SRVF for a piecewise-linear function.
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

  [dim, nsamps] = size(F);

  if ( dim > 1 )
    V = diff(F,1,2);

    for i=1:dim
      V(i,:) = V(i,:) ./ diff( T );
    end

    Vmag = sqrt( sqrt( sum( V .* V ) ) );
    Q = zeros( dim, nsamps-1 );
    for i=1:length(Q)
      if ( Vmag(i) > 1e-4 )
        Q(:,i) = V(:,i) / Vmag(i);
      end
    end
  else
    m = diff(F) ./ diff(T);
    Q = sign(m) .* sqrt(abs(m));
  end
end
