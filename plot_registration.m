function plot_registration( X1, T1, X2, T2, color1, color2 )
  if ( nargin < 6 )
    color2 = 'r';
  end
  if ( nargin < 5 )
    color1 = 'b';
  end
  Tr = union( T1, T2 );

  dim = rows( X1 );
  nsamps = length( Tr );

  for i=1:dim
    X1r(i,:) = interp1( T1, X1(i,:), Tr );
    X2r(i,:) = interp1( T2, X2(i,:), Tr );
  end

  if ( dim == 2 )
    figure();
    hold on;

    plot( X1(1,:), X1(2,:), color1 );
    plot( X2(1,:), X2(2,:), color2 );

    for i=1:3:nsamps
      plot( [X1r(1,i) X2r(1,i)], [X1r(2,i) X2r(2,i)], 'k' );
    end
  elseif ( dim == 1 )
    plot( Tr, X1r, color1, Tr, X2r, color2 );
  else
    error "Unsupported dimension";
  end
end
