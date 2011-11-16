% Given two SRVFs Q1 and Q2, returns the rotation R which best aligns Q2 to Q1.
function R = srvf_optimal_rotation( Q1, T1, Q2, T2, closed )
  dim = rows( Q1 );

  % Is the sort necessary?
  Tr = sort( union( T1, T2 ) );
  Q1r = srvf_refine( Q1, T1, Tr );
  Q2r = srvf_refine( Q2, T2, Tr );

  for i=1:dim
    Q2r(i,:) = Q2r(i,:) .* diff( Tr );
  end

  A = Q1r * Q2r'
  [U,S,V] = svd(A);
  if det(A) > 0
    S = eye(dim);
  else
    S = eye(dim);
    S(:,end) = -S(:,end);
  end
  R = U*S*V';
end
