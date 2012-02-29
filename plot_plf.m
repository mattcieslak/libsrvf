% Plots the PLF represented by F,T.
%
% Inputs
%  F,T : the PLF
%  fmt : a format string recognized by the plot() function
% --------------------------------------------------------------------------
function plot_plf( F, T, color )
  if ( nargin < 3 )
    color = 'b';
  end

  [dim nsegs] = size(F);
  
  if ( dim == 1 )
    plot(T,F,color);
  elseif ( dim == 2 )
    plot(F(1,:),F(2,:),color);
  elseif ( dim == 3 )
    plot3(F(1,:),F(2,:),F(3,:),color);
  else
    error( 'PLFs must be dimension 3 or lower.' );
  end
end
