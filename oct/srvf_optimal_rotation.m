% Given two SRVFs Q1 and Q2, returns the rotation R which best aligns Q2 to Q1.
function R = srvf_optimal_rotation( Q1, T1, Q2, T2 )
  dim = size(Q1,1);

  Tr = unique( [T1 T2] );
  Q1r = srvf_refine( Q1, T1, Tr );
  Q2r = srvf_refine( Q2, T2, Tr );

  for i=1:dim
    Q2r(i,:) = Q2r(i,:) .* diff( Tr );
  end

  A = Q1r * Q2r';
  [U,S,V] = svd(A);
  if det(A) > 0
    S = eye(dim);
  else
    S = eye(dim);
    S(:,end) = -S(:,end);
  end
  R = U*S*V';
end


%!demo
%! load demos/horse-1.mat
%! load demos/horse-2.mat
%! T1 = linspace(0,1,length(X1));
%! T2 = linspace(0,1,length(X2));
%! Q1=plf_to_srvf(X1,T1);
%! Q2=plf_to_srvf(X2,T2);
%! R=srvf_optimal_rotation(Q1,T1,Q2,T2);
%! X2r=R*X2;
%!
%! figure();
%! plot(X1(1,:),X1(2,:),'b',X2(1,:),X2(2,:),'r');
%! xl=min([X1(1,:) X2(1,:)]);
%! xu=max([X1(1,:) X2(1,:)]);
%! yl=min([X1(2,:) X2(2,:)]);
%! yu=max([X1(2,:) X2(2,:)]);
%! axis([xl xu yl yu],'square');
%! title('Curves before rotational alignment');
%!
%! figure();
%! plot(X1(1,:),X1(2,:),'b',X2r(1,:),X2r(2,:),'r');
%! xl=min([X1(1,:) X2r(1,:)]);
%! xu=max([X1(1,:) X2r(1,:)]);
%! yl=min([X1(2,:) X2r(2,:)]);
%! yu=max([X1(2,:) X2r(2,:)]);
%! axis([xl xu yl yu],'square');
%! title('Curves after rotational alignment');
