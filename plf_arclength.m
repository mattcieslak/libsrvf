function [l S] = plf_arclength( F, T, closed=0 )
  [dim npts] = size( F );

  if ( dim > 1 )
    dF = diff(F,1,2);
    if ( closed )
      dF = [dF F(:,1)-F(:,end)];
    end

    S = [0 cumsum( sqrt(sum( dF .* dF, 1 )), 2 )];
    l = S(end);
  else
    dF = diff(F);
    if ( closed ) 
      dF = [dF F(1)-F(end)]; 
    end

    S = [0 cumsum( abs( dF ) )];
    l = S(end);
  end
end
