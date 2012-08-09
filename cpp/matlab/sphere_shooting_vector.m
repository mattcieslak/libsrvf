% Returns the shooting vector from Q1 to Q2.
% --------------------------------------------------------------------------
function [V TV] = sphere_shooting_vector(Q1, T1, Q2, T2)
  Q1nrm = srvf_squared_l2norm(Q1, T1);
  [DQ TDQ] = srvf_linear_combination(Q1, T1, Q2, T2, -1, 1);
  proj_nrm = srvf_l2product(DQ, TDQ, Q1, T1) / Q1nrm;
  d = srvf_preshape_distance(Q1, T1, Q2, T2);
  [V TV] = srvf_linear_combination(DQ, TDQ, Q1, T1, 1, -proj_nrm);
  Vnrm = sqrt(srvf_squared_l2norm(V, TV));
  if (Vnrm > 1e-6)
    V = (d/Vnrm)*V;
  else
    V = zeros(size(Q1,1), length(TV)-1);
  end
end


%!test
%! Q1 = [1];
%! T1 = [0 1];
%! Q2 = [1 -1];
%! T2 = [0 0.5 1];
%! Vexp = (pi/2)*[1 -1];
%! TVexp = [0 0.5 1];
%! [V TV] = sphere_shooting_vector(Q1, T1, Q2, T2);
%! assert(V,Vexp);
%! assert(TV,TVexp);

