% Returns the PLF corrresponding the the given SRVF.  The output PLF 
% will have the same change point vector as the input SRVF, so no change 
% point vector is returned.
%
% Inputs
%  Q,T : the SRVF
%
% Outputs
%  F : the PLF function values
% --------------------------------------------------------------------------
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


%!assert(srvf_to_plf([0],[0 1]),[0 0],1e-6);
%!assert(srvf_to_plf([0;0],[0 1]),[0 0;0 0],1e-6);
%!assert(srvf_to_plf([1],[0 1]),[0 1],1e-6);
%!assert(srvf_to_plf([1;1],[0 1]),[0 sqrt(2);0 sqrt(2)],1e-6);
%!assert(srvf_to_plf([1 -1],[0 1/2 1]),[0 1/2 0],1e-6);
%!assert(srvf_to_plf([1 -1;0 1],[0 1/2 1]),[0 1/2 -.20711;0 0 sqrt(2)/2],1e-4);
