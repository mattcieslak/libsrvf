% Computes the square of the L^2 norm of an SRVF.
%
% Inputs
%  Q,T : the SRVF
%
% Outputs
%  r : the square of the L^2 norm of the SRVF
% --------------------------------------------------------------------------
function r = srvf_squared_l2norm( Q, T )
  r = sum( sum( Q.*Q, 1 ) .* diff( T ) );
end


%!assert(srvf_squared_l2norm([0],[0 1]),0);
%!assert(srvf_squared_l2norm([1],[0 1]),1);
%!assert(srvf_squared_l2norm([-1],[0 1]),1);
%!assert(srvf_squared_l2norm([1 -1],[0 1/2 1]),1);
%!
%!test
%! Q=[0 1 0 -1;\
%!    1 0 1  0];
%! T=[0 1/4 1/2 3/4 1];
%! assert(srvf_squared_l2norm(Q,T),1,1e-6);
