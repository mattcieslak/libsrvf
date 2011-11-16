% Finds a reparametrization gamma minimizing the preshape distance 
% between the SRV functions of p1 and p2(gamma).
%
% Parameters:
%   p1 is the reference curve
%   p2 is the target curve
%   p1 and p2 are assumed to have the same size.
%   Note:  p1 and p2 are the actual curves, not the q-functions.
%   ndp is the desired DP grid size.  Smaller values for faster running 
%   time, larger values for improved accuracy.
% 
% Return:
%   gamma will have the same length as p1 and p2
function gamma = dp( p1, p2, ndp )
  [M N] = size(p1);

  q1 = curve_to_q( p1 );
  q2 = curve_to_q( p2 );

  % dp() needs q2 oversampled, but not q1
  xv  = linspace(0,1,N);
  xvl = linspace(0,1,5*N);
  q1L = spline( xv, q1, xvl );
  q2L = spline( xv, q2, xvl );

  gamma = dp_mex( q1L, q2L, ndp );
  %gamma = spline( linspace(0,1,ndp), gamma, linspace(0,1,N) );
  %gamma = 1/N + (N-1)/N * gamma;
end
