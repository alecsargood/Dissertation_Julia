function [reLambda] = DispersRel_dist(tau, ak, bk, gk, dk, chik, E)

dom = 10;
[D] = funcs_dist(tau, ak, bk, gk, dk, chik, dom, E);

values = roots(D,'complex');
%real_comps = values(:,1);
real_comps = real(values);

reLambda = max(real_comps);


end