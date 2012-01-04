% Returns the PLF corrresponding the the given SRVF.
function F = srvf_to_plf( Q, T )
  [dim nsamps] = size(Q);

  if ( dim > 1 )
    normQ = sqrt( sum( Q.*Q, 1 ) );
    V = Q .* repmat( normQ, dim, 1 ) .* repmat( diff(T), dim, 1 );
    F = [zeros(dim,1) cumsum( V, 2 )];
  else
    v = Q .* abs(Q);
    F = [0 cumsum(v .* diff(T))];
  end
end
