% Useful only for 1-D data (functions).
function [F T] = samps_to_plf( X, TX )
  [dim nsamps] = size(X);

  if ( nargin < 2 )
    TX = linspace(0,1,nsamps);
  end

  if ( dim > 1 )
    % TODO: try to reduce the number of segments for curves?
    T = TX;
    F = X;
  else
    s = updown_sequence( X );
    F = [0 cumsum(s)];
    T = [0 cumsum(abs(s))];
  end
end
