% Computes the composition of two piecewise-linear functions.
% The inner function must be 1-dimensional and non-decreasing.
% The range of the inner function must be contained in the domain 
% of the outer function.
%
% Inputs
% F1, T1:  the outer function
% F2, T2:  the inner function.  Must be a non-decreasing, 1-D function.
% Returns
% F, T:  the composite function
function [F T] = plf_compose( F1, T1, F2, T2 )
  assert(columns(F1)==columns(T1));  % PLF condition for F1
  assert(rows(T1)==1);
  assert(min(diff(T1))>=0);
  assert(columns(F2)==columns(T2));  % PLF condition for F2
  assert(rows(T2)==1);
  assert(min(diff(T2))>=0);
  assert(rows(F2)==1);       % F2 must be 1-D
  assert(min(diff(F2))>=0);  % F2 must be non-decreasing
  assert(F2(1)>T1(1)-1e-6 && F2(end)<T1(end)+1e-6);  % Functions composable?

  % To get the change points of the inner function, we take T2, along 
  % with the preimages of T1 under F2.
  T = union( T2, plf_preimages( F2, T2, T1 ) );

  % Evaluate the composite function at all of these points
  for i=1:rows(F1)
    F(i,:) = interp1(T1,F1(i,:),interp1( T2,F2,T,'extrap' ),'extrap');
  end
end


%!test
%! F1=linspace(0,1,5);
%! T1=linspace(0,1,5);
%! F2=linspace(0,1,5);
%! T2=linspace(0,1,5);
%! [F T]=plf_compose(F1,T1,F2,T2);
%! assert(F,F1,eps);
%! assert(T,T1,eps);
%!
%!test
%! F1=[0 3/4 1];
%! T1=[0 1/4 1];
%! F2=[0 1/4 1];
%! T2=[0 3/4 1];
%! [F T]=plf_compose(F1,T1,F2,T2);
%! assert(F,[0 3/4 1],1e-4);
%! assert(T,[0 3/4 1],1e-4);
%! [F T]=plf_compose(F2,T2,F1,T1);
%! assert(F,[0 1/4 1],1e-4);
%! assert(T,[0 1/4 1],1e-4);
%!#figure();
%!#plot(T1,F1,'b',T2,F2,'r',T,F,'k');
%!
%!test
%! F1=[0 1 1 0 0; \
%!     0 0 1 1 0];
%! T1=linspace(0,1,5);
%! F2=[0 1/3 1];
%! T2=[0 2/3 1];
%! [F T]=plf_compose(F1,T1,F2,T2);
%! Texp=[0 1/2 2/3 3/4 7/8 1];
%! Fexp=[0 1 1   1 0 0;
%!       0 0 1/3 1 1 0];
%! assert(T,Texp,1e-4);
%! assert(F,Fexp,1e-4);
%!#figure();
%!#plot(F(1,:),F(2,:),'b-*');
%!
%!test
%! F1=[0 3/8 5/8 1];
%! T1=[0 3/8 5/8 1];
%! F2=[0 0 1/4 1/2 3/4 3/4 1];
%! T2=linspace(0,1,7);
%! [F T]=plf_compose(F1,T1,F2,T2);

