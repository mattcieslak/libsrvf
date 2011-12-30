function q = srvf_to_samps( Q, T, TS )
  if ( length(TS) == 1 )
    tv = linspace(0,1,TS);
  else
    tv = TS;
  end

  q = [];
  Qidx = 1;
  nsamps = length(tv);

  for qidx = 1:nsamps
    while ( Qidx < length(T)-1 && tv(qidx) > T(Qidx+1) ) 
      Qidx = Qidx + 1; 
    end

    q = [q Q(:,Qidx)];
  end
end
