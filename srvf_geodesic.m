% Computes several steps along the minimum-length great circle path between 
% the given SRVFs Q1,T1 and Q2,T2.  The SRVFs must have the same dimension 
% and the same L^2 norm, and must be defined on the same interval.
%
% Inputs:
%  Q1,T1  : the SRVF at the beginning of the geodesic
%  Q2,T2  : the SRVF at the end of the geodesic
%  Nsteps : the number of points (including the end points) along 
%           the geodesic to return
%
% Outputs
%  P : a 3-dimensional array containing the points on the geodesic.
%      P(:,:,k) is an SRVF representing the kth point.  All of these 
%      SRVFs will have the same change point parameters.
%  T : a 1-row matrix containing the change point parameters for the SRVFs in P.
function [P T] = srvf_geodesic( Q1, T1, Q2, T2, Nsteps )
  assert( rows(Q1) == rows(Q2) );
  assert( rows(T1) == 1 && columns(T1) == columns(Q1)+1 );
  assert( rows(T2) == 1 && columns(T2) == columns(Q2)+1 );
  assert( min(diff(T1)) > 0 );
  assert( min(diff(T2)) > 0 );
  assert( abs(T1(1)-T2(1)) < 1e-4 );
  assert( abs(T1(end)-T2(end)) < 1e-4 );

  T = unique( [T1 T2] );
  Q1r = srvf_refine( Q1, T1, T );
  Q2r = srvf_refine( Q2, T2, T );
  [dim nsegs] = size(Q1r);

  P = zeros( dim, nsegs, Nsteps );
  Ptv = linspace( 0, 1, Nsteps );

  ip = srvf_l2product( Q1r, T, Q2r, T );
  theta = acos( ip );

  for i=1:Nsteps
    P(:,:,i)=(sin(theta*(1-Ptv(i)))*Q1r + sin(theta*Ptv(i))*Q2r) ./ sin(theta);
  end
end


%!demo
%! load demos/horse-1.mat
%! load demos/horse-1a.mat
%! [F1 T1]=samps_to_plf(X1);
%! [F2 T2]=samps_to_plf(X2);
%! Q1=plf_to_srvf(F1,T1);
%! Q2=plf_to_srvf(F2,T2);
%! [G GT]=srvf_optimal_matching(Q1,T1,Q2,T2);
%! [F2r T2r]=plf_compose(F2,T2,G,GT);
%! Q2r=plf_to_srvf(F2r,T2r);
%! [P T] = srvf_geodesic(Q1,T1,Q2r,T2r,5);
%! plot_geodesic(P,T);
%! title("The curves corresponding to the geodesic");

%!#demo
%! load demos/rna1.mat
%! load demos/rna2.mat
%! [F1 T1]=samps_to_plf(X1);
%! [F2 T2]=samps_to_plf(X2);
%! Q1=plf_to_srvf(F1,T1);
%! Q2=plf_to_srvf(F2,T2);
%! [G GT]=srvf_optimal_matching(Q1,T1,Q2,T2);
%! [F2r T2r]=plf_compose(F2,T2,G,GT);
%! Q2r=plf_to_srvf(F2r,T2r);
%! [P T] = srvf_geodesic(Q1,T1,Q2r,T2r,5);
%! plot_geodesic(P,T);
%! title("The curves corresponding to the geodesic");
%! plot_geodesic(P,T,'q');
%! title("The SRVFs on the geodesic");
