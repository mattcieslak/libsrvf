% Computes the Karcher mean of a collection of points on the unit sphere
% in L^2.
% --------------------------------------------------------------------------
function [Qm Tm] = sphere_karcher_mean(Qs, Ts, tol, max_iters, stp)
  if (nargin < 3) tol=1e-3; end
  if (nargin < 4) max_iters=100; end
  if (nargin < 5) stp=0.3; end

  Qm = Qs{1};
  Tm = Ts{1};

  for iter = 1:max_iters
    V=zeros(size(Qm,1),1);
    TV=[0 1];

    for i=1:length(Qs)
      % Rotate, reparametrize, and rotate Qi to align to Qm
      Qitmp = Qs{i};
      TQitmp = Ts{i};

      R = srvf_optimal_rotation(Qm, Tm, Qitmp, TQitmp);
      Qitmp = R*Qitmp;

      [G T] = srvf_optimal_matching(Qm, Tm, Qitmp, TQitmp);
      [Qitmp TQitmp] = srvf_gamma_action(Qitmp, TQitmp, G, T);

      R = srvf_optimal_rotation(Qm, Tm, Qitmp, TQitmp);
      Qitmp = R*Qitmp;

      % Add current shooting vector to V
      [Vcur TVcur] = sphere_shooting_vector(Qm, Tm, Qitmp, TQitmp);
      [V TV] = srvf_linear_combination(V, TV, Vcur, TVcur, 1, 1);
    end

    % Stop if norm of V < tol
    V = V / length(Qs);
    Vnrm = sqrt(srvf_squared_l2norm(V, TV))
    if (Vnrm < tol) break; end

    % Otherwise, update Qm
    w1 = cos(stp*Vnrm);
    w2 = sin(stp*Vnrm) / Vnrm;
    [Qm Tm] = srvf_linear_combination(Qm, Tm, V, TV, w1, w2);
  end
end

%!test
%! Qs{1}=[1];
%! Ts{1}=[0 1];
%! Qs{2}=[1 -1];
%! Ts{2}=[0 0.5 1];
%! Qmexp=[0];
%! Tmexp=[0 1];
%! [Qm Tm]=sphere_karcher_mean(Qs,Ts);

%!demo
%! colors = {'b', 'g', 'c', 'm', 'r'};
%! nfuncs=3;
%! figure();
%! hold on;
%! for i=1:nfuncs
%!   fs{i} = randfunc(12, 100, 0);
%!   [Fs{i} Ts{i}] = poly_to_plf(fs{i});
%!   Qs{i} = plf_to_srvf(Fs{i}, Ts{i});
%!   plot(Ts{i}, Fs{i}, colors{i});
%! end
%! [Qm Tm] = sphere_karcher_mean(Qs, Ts);
%! Fm = srvf_to_plf(Qm, Tm);
%! figure();
%! hold on;
%! plot(Tm, Fm, 'k;Fm;');
%! for i=1:nfuncs
%!   [G T] = srvf_optimal_matching(Qm, Tm, Qs{i}, Ts{i});
%!   [Qrs{i} Trs{i}] = srvf_gamma_action(Qs{i}, Ts{i}, G, T);
%!   Frs{i} = srvf_to_plf(Qrs{i}, Trs{i});
%!   plot(Trs{i}, Frs{i}, colors{i});
%! end
%!
%!demo
%! load demos/rna1.mat
%! load demos/rna2.mat
%! nfuncs=2;
%! colors = {'b', 'g', 'c', 'm', 'r'};
%! 
%! [Fs{1}, Ts{1}] = poly_to_plf(X1);
%! [Fs{2}, Ts{2}] = poly_to_plf(X2);
%! figure();
%! hold on;
%! for i=1:nfuncs
%!   Qs{i} = plf_to_srvf(Fs{i}, Ts{i});
%!   plot3(Fs{i}(1,:), Fs{i}(2,:), Fs{i}(3,:), colors{i});
%! end
%! [Qm, Tm] = sphere_karcher_mean(Qs, Ts, 1e-3, 10, 0.3);
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
