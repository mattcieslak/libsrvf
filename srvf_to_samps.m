function q = srvf_to_samps( Q, T, nsamps )
  q = [];
  tv = linspace( 0, 1, nsamps );

  Qidx = 1;

  for qidx = 1:nsamps
    while ( Qidx < length(T)-1 && tv(qidx) > T(Qidx+1) ) 
      Qidx = Qidx + 1; 
    end

    q = [q Q(Qidx)];
  end
end
