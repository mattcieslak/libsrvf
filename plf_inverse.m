function [Fi Ti] = plf_inverse( F, T, min_diff )
  if ( nargin < 3 )
    min_diff = 1e-8;
  end

  Ti = make_monotone( F, min_diff );
  Fi = make_monotone( T, min_diff );

  Ti = Ti / sum(diff(Ti));
end

function v = make_monotone( v, min_diff )
  d = guess_direction( v );
  assert( abs(d) > 0.5 );
  
  for i=2:length( v )
    if ( d > 0 ) 
      v(i) = max( v(i), v(i-1) + min_diff );
    else
      v(i) = min( v(i), v(i-1) - min_diff );
    end
  end
end

function d = guess_direction( v )
  d = sign( sum( v ) );
end
