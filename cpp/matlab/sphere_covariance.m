% Computes the covariance of a collection of points on the unit sphere 
% in L^2.
%
% Inputs:
%  Qs : the collection of SRVFs.  Must all have the same domain, and 
%       must already be aligned to Qm.
%  Ts : the corresponding parameter values for Qs
%  Qm : the Karcher mean of the functions in Qs.  Must have same domain as 
%       the SRVFs in Qs.
%  Ts : the corresponding parameter values for Qm
% --------------------------------------------------------------------------
function [K, Mu] = sphere_covariance(Qs, Ts, Qm, Tm)
  uv = linspace(Tm(1), Tm(end), length(Tm));

  % TODO: to get a bunch of vectors of the same length, we're evaluating 
  % all of the SRVF's at uniformly-spaced parameters.  Almost certainly 
  % not the best way to do this.
  Dat = zeros(length(Qs), size(Qm,1)*length(uv));
  for i=1:length(Qs)
    [V TV] = sphere_shooting_vector(Qm, Tm, Qs{i}, Ts{i});
    Vsamps = srvf_evaluate(V, TV, uv);
    Dat(i,:) = Vsamps'(:);
  end

  K = cov(Dat);
  Mu = mean(Dat);
end


%!demo
%! debug_on_error(1);
%! load demos/rna1.mat
%! load demos/rna2.mat
%! load demos/rna3.mat
%! load demos/rna4.mat
%! load demos/rna5.mat
%! nfuncs=5;
%! colors = {'b', 'g', 'c', 'm', 'r'};
%! 
%! [Fs{1}, Ts{1}] = poly_to_plf(X1);
%! [Fs{2}, Ts{2}] = poly_to_plf(X2);
%! [Fs{3}, Ts{3}] = poly_to_plf(X3);
%! [Fs{4}, Ts{4}] = poly_to_plf(X4);
%! [Fs{5}, Ts{5}] = poly_to_plf(X5);
%! figure();
%! hold on;
%! for i=1:nfuncs
%!   Ts{i} = Ts{i} / Ts{i}(end);
%!   Qs{i} = plf_to_srvf(Fs{i}, Ts{i});
%!   nrm = srvf_squared_l2norm(Qs{i}, Ts{i});
%!   Qs{i} = Qs{i} / sqrt(nrm);
%!   Fs{i} = Fs{i} / nrm;
%!   plot3(Fs{i}(1,:), Fs{i}(2,:), Fs{i}(3,:), colors{i});
%! end
%! [Qm, Tm] = sphere_karcher_mean(Qs, Ts, 1e-3, 30, 0.3);
%! Fm = srvf_to_plf(Qm, Tm);
%! figure();
%! hold on;
%! plot3(Fm(1,:), Fm(2,:), Fm(3,:), 'k;Fm;');
%! for i=1:nfuncs
%!   R = srvf_optimal_rotation(Qm, Tm, Qs{i}, Ts{i});
%!   Qs{i} = R*Qs{i};
%!   [G T] = srvf_optimal_matching(Qm, Tm, Qs{i}, Ts{i});
%!   [Qs{i}, Ts{i}] = srvf_gamma_action(Qs{i}, Ts{i}, G, T);
%!   R = srvf_optimal_rotation(Qm, Tm, Qs{i}, Ts{i});
%!   Qs{i} = R*Qs{i};
%!   Fsr{i} = srvf_to_plf(Qs{i}, Ts{i});
%!   plot3(Fsr{i}(1,:), Fsr{i}(2,:), Fsr{i}(3,:), colors{i});
%! end
%! [K, Mu] = sphere_covariance(Qs, Ts, Qm, Tm)
