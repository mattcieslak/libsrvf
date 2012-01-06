% Plots the SRVF represented by Q,T.
%
% Inputs
%  Q,T : the SRVF
%  fmt : a format string recognized by the plot() function
% --------------------------------------------------------------------------
function plot_srvf( Q, T, fmt )
  if ( nargin < 3 )
    fmt = 'b';
  end

  [dim nsegs] = size(Q);
  
  if ( dim == 1 )
    hold on;
    for i=1:length(Q)
      plot([T(i) T(i+1)], [Q(i), Q(i)], fmt );
    end
  elseif ( dim == 2 )
    plot(Q(1,:),Q(2,:),fmt);
  elseif ( dim == 3 )
    plot3(Q(1,:),Q(2,:),Q(3,:),fmt);
  else
    error( "SRVFs must be dimension 3 or lower." );
  end
end
