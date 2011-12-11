function plot_geodesic( G, T, plot_type, colors )
  if ( nargin < 3 )
    plog_type = 'f';  % default: plot functions
  end
  if ( nargin < 4 )
    colors = { 'k', 'b', 'c', 'm', 'r' };
  end

  for i=1:size(G,1)
    cidx = mod( i-1, 5 ) + 1;
    if ( nargin >= 3 && plot_type == 'q' )
      plot_srvf( G(i,:), T, colors{cidx} );
    else
      GFi = srvf_to_plf( G(i,:), T );
      if ( i==1 || i==size(G,1) )
        plot( T, GFi, sprintf('%s;f%d;', colors{cidx}, i) );
      else
        plot( T, GFi, colors{cidx} );
      end
    end
  end
end
