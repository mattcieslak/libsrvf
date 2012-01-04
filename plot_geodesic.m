function plot_geodesic( G, T, plot_type, colors )
  if ( nargin < 3 )
    plot_type = 'f';  % default: plot functions
  end
  if ( nargin < 4 )
    colors = { 'k', 'b', 'c', 'm', 'r' }; % default colors
  end

  [dim nsegs nsteps] = size(G);

  for i=1:nsteps
    cidx = mod( i-1, 5 ) + 1;  % cycle through the colors
    if ( nargin >= 3 && plot_type == 'q' )
      plot_srvf( G(i,:), T, colors{cidx} );
    else
      GFi = srvf_to_plf( G(:,:,i), T );
      if ( i==1 || i==nsteps )
        keystr = sprintf('%s;f%d;', colors{cidx}, i);
      else
        keystr = colors{cidx};
      end
      if ( dim == 1 )
        plot( T, GFi, keystr );
      elseif ( dim == 2 )
        plot( GFi(1,:), GFi(2,:), keystr );
      elseif ( dim == 3 )
        plot3( GFi(1,:), GFi(2,:), GFi(3,:), keystr );
      else
        error( "Can't plot higher than 3-dimensional curves." );
      end
    end
  end
end
