% Useful only for 1-D data (functions).
function [F T S] = samps_to_plf( X, TX )
  [dim nsamps] = size(X);

  if ( nargin < 2 )
    TX = linspace(0,1,nsamps);
  end

  if ( dim > 1 )
    % TODO: try to reduce the number of segments for curves?
    T = TX;
    F = X;
    if ( nargout > 2 )
      error "S is only useful for 1-D curves";
    end
  else
    S = updown_sequence( X );
    F = [0 cumsum(S)];
    T = [0 cumsum(abs(S))];
  end
end
