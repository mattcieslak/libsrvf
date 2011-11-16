function Xn = preprocess( X, T, closed=0 )
  [dim nsamps] = size( X );

  if ( nargin < 2 )
    if ( closed )
      T = linspace(0, 1, nsamps+1)(1:end-1);
    else
      T = linspace(0, 1, nsamps);
    end
  end

  if ( dim > nsamps )
    X = X';
    [dim nsamps] = size( X );
  end

  l = plf_arclength( X, T, closed );
  if ( l > 1e-5 )
    Xn = X / l;
  else
    Xn = zeros(dim, nsamps);
  end
end
