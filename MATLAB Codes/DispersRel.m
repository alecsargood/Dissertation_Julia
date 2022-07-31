function [reLambda] = DispersRel(tau, ak, bk, gk, dk, chik)

dom = 10;
[D] = funcs(tau, ak, bk, gk, dk, chik, dom);

values = roots(D,'complex');
%real_comps = values(:,1);
real_comps = real(values);

reLambda = max(real_comps);


end