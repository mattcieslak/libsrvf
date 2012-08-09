% Computes the action of a 1-D diffeomorphism on the given SRVF.
%
% Inputs:
%  Q  - SRVF function values
%  TQ - SRVF parameter values
%  G  - diffeomorphism function values.  Must be increasing, and must 
%       satisfy TQ(1)<=G(1) and TQ(end)>=G(end).
%  TG - diffeomorphism parameter values.
%
% Outputs:
%  Qr - the new SRVF function values
%  Tr - the new SRVF parameter values
% --------------------------------------------------------------------------
function [Qr Tr] = srvf_gamma_action(Q,TQ,G,TG)
  TGx = plf_preimages(G,TG,TQ);
  Tr = unique([TG TGx]);
  Gr = unique([TQ G]);
  DGr = (Gr(2:end)-Gr(1:(end-1))) ./ (Tr(2:end)-Tr(1:(end-1)));
  Qr = srvf_refine(Q,TQ,Gr);
  Qr = Qr .* repmat(sqrt(DGr),size(Qr,1),1);
end


%!test
%! Q=[1];
%! TQ=[0 1];
%! G=[0 1/3 1];
%! TG=[0 1/2 1];
%! [Qr Tr]=srvf_gamma_action(Q,TQ,G,TG);
%! Qrexp=[sqrt(2/3) sqrt(4/3)];
%! Trexp=[0 1/2 1];
%! assert(Qr,Qrexp,1e-6);
%! assert(Tr,Trexp,1e-6);
%!
%!test
%! Q=[1 -1];
%! TQ=[0 1/2 1];
%! G=[0 1/3 1];
%! TG=[0 2/3 1];
%! [Qr Tr]=srvf_gamma_action(Q,TQ,G,TG);
%! Qrexp=[sqrt(1/2) sqrt(2) -sqrt(2)];
%! Trexp=[0 2/3 3/4 1];
%! assert(Qr,Qrexp,1e-6);
%! assert(Tr,Trexp,1e-6);
