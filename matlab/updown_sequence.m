function s = updown_sequence( f )
  df = diff( f );

  s = [df(1)];

  for i=2:length(df)
    if ( sign(df(i)) == sign(s(end)) )
      s(end) = s(end) + df(i);
    else
      s = [s df(i)];
    end
  end
end
