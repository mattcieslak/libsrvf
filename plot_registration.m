function plot_registration( X1, X2 )
  figure();
  hold on;

  plot( X1(1,:), X1(2,:), 'b' );
  plot( X2(1,:), X2(2,:), 'r' );

  for i=1:length(X1)
    plot( [X1(1,i) X2(1,i)], [X1(2,i) X2(2,i)], 'k' );
  end
end
