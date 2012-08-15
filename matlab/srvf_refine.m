% Given a piecewise-constant SRVF Q with changepoint vector T, 
% and a new changepoint vector Tr which is a refinement of T, 
% compute the values qr corresponding to Tr.  Basically, this 
% just amounts to duplicating an element of Q whenever the corresponding 
% interval is split.
%
% Inputs
%  Q,T : the SRVF
%  Tr : a refinement of T
%
% Outputs
%  Qr : the SRVF values corresponding to the new changepoint vector
% --------------------------------------------------------------------------
function Qr = srvf_refine( Q, T, Tr )
  %Qr = [];

  %idx1 = 1;

  %for idx2 = 1:(length(Tr)-1)
  %  while ( Tr(idx2+1) > T(idx1+1) + 1e-6 )
  %    idx1 = idx1 + 1;
  %  end
  %  Qr = [Qr Q(:,idx1)];
  %end

  tv = (Tr(1:end-1)+Tr(2:end))/2;
  Qr = srvf_evaluate(Q,T,tv);
end


%!test
%! Q=[1];
%! T=[0 1];
%! Tr=[0 1/4 1/2 3/4 1];
%! Qrexp=[1 1 1 1];
%! Qr=srvf_refine( Q, T, Tr );
%! assert( Qr, Qrexp );
%!
%!test
%! Q=[1 -1 1 -1 1;\
%!    1 -1 1 -1 1];
%! T=[0 1/5 2/5 3/5 4/5 1];
%! Tr=unique([T 0.19 0.2001 0.65 0.999]);
%! Qrexp=[1 1 -1 -1 1 -1 -1 1 1;\
%!        1 1 -1 -1 1 -1 -1 1 1];
%! Qr=srvf_refine( Q, T, Tr );
%! assert( Qr, Qrexp );
