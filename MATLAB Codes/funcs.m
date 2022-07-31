function [D] = funcs(tau, ak, bk, gk, dk, chik, dom)
%Chebfun2 area to look for roots in
d = [-dom, dom, -0.1, dom];
d1 = [-dom, dom];
% Splitting characteristic eqn into Re and Im parts.

% Defining functions for ease of notation

E=@(x)exp(-x*tau);


D = chebfun(@(l)(l^2+(ak*l)+bk+((gk*l + dk)*E(l)+(chik*E(l)^2))), d1);
end