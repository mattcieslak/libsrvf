% Selects a single match from the set of Pareto-optimal partial matches.
% Inputs:
%  T1 - parameters for the first curve
%  T2 - parameters for the second curve
%  P - a Kx5 matrix representing the Pareto set, as returned by 
%      srvf_pmatch_find_matches()
%  thresh - (optional) specify the linear regression error threshold 
%      used to find the knee of the Pareto frontier curve.
%
% Outputs:
%  idx - the index of the selected match
% --------------------------------------------------------------------------
function idx = srvf_pmatch_select_single(T1,T2,P,thresh)
  if (nargin < 4) 
    thresh = 0.003; 
  end

  nmatches = size(P,1);
  tv = zeros(1,nmatches);
  yv = zeros(1,nmatches);

  for i=1:nmatches
    a = T1(P(i,1)); b = T1(P(i,2));
    c = T2(P(i,3)); d = T2(P(i,4));

    tv(i) = b-a + d-c;
    yv(i) = P(i,5);
  end

  nst = 0;
  tdy = 0;
  idx = 1;  % the return value

  for i=2:size(P,1)
    nst = nst + tv(i)*tv(i);
    tdy = tdy + tv(i)*yv(i);

    if (nst < 1e-6)
      continue;
    end

    m = tdy / nst;  % slope of current regression line
    err = 0;        % sum of squared residuals
    for j=1:i
      dy = yv(j) - m*tv(j);
      err = err + dy*dy;
    end

    if (err > thresh)
      break;
    else
      idx = i;
    end
  end
end
