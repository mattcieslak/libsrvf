% F1, T1:  the outer function
% gamma, Tgamma  the inner function (must be 1-D)
function [F T] = plf_compose( F1, T1, gamma, Tgamma )
  % T = change points of the composite function.
  % Start with T = change points of the inner function ...
  T = Tgamma;

  % ... and add in the points where the inner function hits a 
  % change point of the outer function.
  for i=1:(length(Tgamma)-1)
    if ( abs( gamma(i) - gamma(i+1) ) > 1e-5 ) 
      v = find_in_range( T1, [gamma(i) gamma(i+1)] );
      T = union( T, interp1( [gamma(i) gamma(i+1)], \
                    [Tgamma(i) Tgamma(i+1)], v, 'extrap' ) );
    end
  end

  % Evaluate the composite function at all of these points
  for i=1:rows(F1)
    F(i,:) = interp1(T1,F1(i,:),interp1( Tgamma,gamma,T,'extrap' ),'extrap');
  end
end

function w = find_in_range( v, R )
  R = sort( R );

  w = [];
  for i=1:length(v)
    if ( v(i) >= R(1) && v(i) <= R(2) )
      w = [w v(i)];
    end
  end
end
