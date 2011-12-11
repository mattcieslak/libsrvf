function plot_srvf( Q, T, color )
  if ( nargin < 3 )
    color = 'b';
  end

  hold on;
  for i=1:length(Q)
    plot([T(i) T(i+1)], [Q(i), Q(i)], color );
  end
end
