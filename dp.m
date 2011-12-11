% Given two SRVFs q1 and q2, finds a reparametrization gamma that 
% minimizes the distance between q1 and q2*gamma (the action of 
% gamma on q2).
%
% Parameters:
%  Q1 - reference SRVF function values
%  T1 - reference SRVF changepoint parameter values
%  Q2 - target SRVF function values
%  T2 - target SRVF changepoint parameter values
%  ndp - DP grid size
% 
% Returns:
%  g - gamma function values
%  Tg - gamma changepoint parameter values
function [G T] = dp( Q1, T1, Q2, T2, ndp )
  [G T] = dp_mex( Q1, T1, Q2, T2, ndp );
end
