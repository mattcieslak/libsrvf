function Fc=plf_trans_to_origin(F,T);
  [dim npts]=size(F);
  ctr=eye(npts) - ones(npts) * (1/npts);
  Fc=F*ctr;
end
