function r = srvf_squared_l2norm( Q, T )
  r = sum( sum( Q.*Q, 1 ) .* diff( T ) );
end
