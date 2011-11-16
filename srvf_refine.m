% Given a piecewise-constant function q with changepoint vector T, 
% and a new changepoint vector Tr which is a refinement of T, 
% compute the values qr corresponding to Tr.  Basically, this 
% just amounts to duplicating an element of q whenever the corresponding 
% interval is split.
function Qr = srvf_refine( Q, T, Tr )
  Qr = [];

  idx1 = 1;

  for idx2 = 1:(length(Tr)-1)
    while ( Tr(idx2+1) > T(idx1+1) + 1e-4 )
      idx1 = idx1 + 1;
    end
    Qr = [Qr Q(:,idx1)];
  end
end
