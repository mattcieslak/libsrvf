function f = plf_to_samps( F, T, t )
  if nargin < 3
    t = linspace( 0, 1, 100 );
  elseif length(t) == 1
    t = linspace(0,1,t);
  end

  f = interp1( T, F, t );
end
