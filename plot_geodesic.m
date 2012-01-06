% Plots the geodesic represented by G,T.  G is a 3-dimensional array 
% containing the SRVFs along the geodesic.  G(:,:,k) represents the kth
% SRVF, and all SRVFs have the same change point parameters.
%
% This function plots to the current plot window.  If you want to see 
% all of the functions at the same time, you'll need to set hold on before 
% calling this routine.
%
% This routine can either plot the SRVFs themselves, or it can plot 
% the corresponding curves.  The default is to plot curves; set the argument 
% plot_type to 'q' to draw SRVFs instead.
%
% The color of each curve (or SRVF) along the geodesic is determined by 
% the argument colors, which is a cell array containing at least one color 
% string recognized by the plot() function.  The default is 
% { 'k', 'b', 'c', 'm', 'r' }.  The first curve / SRVF gets colors{1}, 
% the second gets colors{2}, and so on.  Repeats cyclically if there are 
% more curves / SRVFs than colors.
%
% Inputs
%  G : the SRVFs.  G(:,:,k) represents the kth SRVF.
%  T : the change point parameters for the SRVFs in G
%  plot_type : 'f' for functions (the default), or 'q' for SRVFs
%  colors : a cell array containing color strings.  Default is 
%           colors = { 'k', 'b', 'c', 'm', 'r' }.
% --------------------------------------------------------------------------
function plot_geodesic( G, T, plot_type, colors )
  if ( nargin < 3 )
    plot_type = 'f';  % default: plot functions
  end
  if ( nargin < 4 )
    colors = { 'k', 'b', 'c', 'm', 'r' }; % default colors
  end

  [dim nsegs nsteps] = size(G);

  for i=1:nsteps
    cidx = mod( i-1, length(colors) ) + 1;  % cycle through the colors
    if ( i==1 || i==nsteps )
      keystr = sprintf('%s;f%d;', colors{cidx}, i);
    else
      keystr = colors{cidx};
    end

    if ( plot_type == 'q' )
      plot_srvf( G(:,:,i), T, keystr );
    else
      GFi = srvf_to_plf( G(:,:,i), T );
      plot_plf( GFi, T, keystr );
    end
  end
end
