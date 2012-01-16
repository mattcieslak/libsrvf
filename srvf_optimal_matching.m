% Given two SRVFs q1 and q2, finds a reparametrization gamma that 
% minimizes the distance between q1 and q2*gamma (the action of 
% gamma on q2).
%
% This function calls a C implementation of a matching algorithm due to 
% Dmitri Bertsekas.
%
% If tv1 is a vector, then the DP grid columns will be spaced according 
% to the values in tv1.  If tv1 is a scalar, then the DP grid will have 
% tv1 uniformly-spaced columns (i.e. linspace(0,1,tv1)).  tv2 is similar.
% If tv1 and tv2 are omitted, the default is tv1=T1 and tv2=T2.
%
% The optional parameter nseeds can be used to specify how many 
% different seedpoints should be tried.  This is only useful for closed
% curves.  If nseeds is omitted, or if nseeds <= 1, then no seedpoint 
% optimization is done.
%
% Parameters:
%  Q1     - reference SRVF function values
%  T1     - reference SRVF changepoint parameter values
%  Q2     - target SRVF function values
%  T2     - target SRVF changepoint parameter values
%  tv1    - DP grid column parameter values, or grid width
%           (default is tv1 = T1).
%  tv2    - DP grid row parameter values, or grid height
%           (default is tv2 = T2).
%  nseeds - number of seedpoints to try.  Set nseeds <= 1 for no seedpoint 
%           optimization (default is 1).
% 
% Returns:
%  G - gamma function values
%  T - gamma changepoint parameter values
%  I - best seed index
function [G T seed] = srvf_optimal_matching( Q1, T1, Q2, T2, tv1, tv2, nseeds )
  % Validation
  assert(nargin>=4);
  assert(rows(Q1)==rows(Q2));
  assert(length(T1)==length(Q1)+1);
  assert(length(T2)==length(Q2)+1);
  assert(rows(T1)==1);
  assert(rows(T2)==1);

  % Argument tv1
  if ( nargin < 5 ) 
    tv1=T1; 
  elseif ( isscalar(tv1) )
    tv1=linspace(0,1,tv1);
  end

  % Argument tv2
  if ( nargin < 6 ) 
    tv2=T2; 
  elseif ( isscalar(tv2) )
    tv2=linspace(0,1,tv2);
  end

  % Argument nseeds
  if ( nargin < 7 || nseeds < 2 )
    nseeds = 1;
  end
  assert( nseeds <= length(T2) );

  if ( nseeds > 1 )
    % Initialize return to the identity function
    G = T2;
    T = T2;
    seed = 1;

    % Do the elastic matching for each seedpoint, updating G, T, and seed 
    % whenever a new minimum distance is found.
    min_dist = 1e9;
    for seedidx = 1:floor(length(T2)/nseeds):length(T2)
      % Change seedpoint
      T2cur = [T2(seedidx:end) T2(1:(seedidx-1))];
      Q2cur = [Q2(:,seedidx:end) Q2(:,1:(seedidx-1))];

      % Do the elastic matching
      [Gcur Tcur cur_dist] = dp_mex( Q1, T1, Q2cur, T2cur, tv1, tv2 );

      if ( cur_dist < min_dist )
        G = Gcur;
        T = Tcur;
        seed = seedidx;
        min_dist = cur_dist;
      end
    end
  else
    [G T] = dp_mex( Q1, T1, Q2, T2, tv1, tv2 );
    seed = 1;
  end
end


%!demo
%! load demos/horse-1.mat
%! load demos/horse-2.mat
#%! load demos/1B0W.mat
#%! load demos/1BWW.mat
#%! uv = linspace(0,1,100);
#%! X1 = [uv; sin(pi*uv)];
#%! X2 = [uv; sin(2*pi*uv(1:50)) zeros(1,length(uv)/2)];
%!
%! X1 = [X1 X1(:,1)];
%! X2 = [X2 X2(:,1)];
%! T1 = linspace(0,1,length(X1));
%! T2 = linspace(0,1,length(X2));
%! printf('Arclengths: %0.3f, %0.3f\n', \
%!         plf_arclength(X1,T1), \
%!         plf_arclength(X2,T2));
%!
%! plot_registration( X1, T1, X2, T2, 'b-*', 'r-*' );
%! title('Curves before resampling');
%!
%! X1 = Resample_Uniform_Increase_4(X1,length(X1));
%! X2 = Resample_Uniform_Increase_4(X2,length(X2));
%!
%! plot_registration( X1, T1, X2, T2, 'b-*', 'r-*' );
%! title('Curves after resampling');
%!
%! Q1=plf_to_srvf(X1,T1);
%! Q2=plf_to_srvf(X2,T2);
%!
%! n1 = sqrt(srvf_squared_l2norm(Q1,T1));
%! n2 = sqrt(srvf_squared_l2norm(Q2,T2));
%! Q1=Q1 / n1;
%! X1=X1 / n1;
%! Q2=Q2 / n2;
%! X2=X2 / n2;
%! printf('L2 Norms: %0.3f, %0.3f\n',\
%!         srvf_squared_l2norm(Q1,T1), \
%!         srvf_squared_l2norm(Q2,T2));
%! R=srvf_optimal_rotation(Q1,T1,Q2,T2);
%! X2r=R*X2;
%! Q2r=R*Q2;
%! 
%! [G TG]=srvf_optimal_matching(Q1,T1,Q2r,T2);
%! [X2rr T2rr]=plf_compose(X2r,T2,G,TG);
%! Q2rr=plf_to_srvf(X2rr,T2rr);
%! Q2rr=Q2rr / sqrt(srvf_squared_l2norm(Q2rr,T2rr));
%!
%! plot_registration(X1,T1,X2r,T2);
%! xl=min([X1(1,:) X2r(1,:)]);
%! xu=max([X1(1,:) X2r(1,:)]);
%! yl=min([X1(2,:) X2r(2,:)]);
%! yu=max([X1(2,:) X2r(2,:)]);
%! axis([xl xu yl yu],"square");
%! title("Original matching");
%!
%! plot_registration(X1,T1,X2rr,T2rr);
%! xl=min([X1(1,:) X2rr(1,:)]);
%! xu=max([X1(1,:) X2rr(1,:)]);
%! yl=min([X1(2,:) X2rr(2,:)]);
%! yu=max([X1(2,:) X2rr(2,:)]);
%! axis([xl xu yl yu],"square");
%! title("Matching after DP");
%!
%! figure();
%! plot(TG,G);
%!
%! dist1 = srvf_squared_l2distance(Q1,T1,Q2,T2);
%! dist2 = srvf_squared_l2distance(Q1,T1,Q2r,T2);
%! dist3 = srvf_squared_l2distance(Q1,T1,Q2rr,T2rr);
%! printf("Distances:\n");
%! printf("\tOriginal: %0.3f\n",dist1);
%! printf("\tRotated: %0.3f\n",dist2);
%! printf("\tRotated and Matched: %0.3f\n",dist3);
