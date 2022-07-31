%% Select model

% Set desired model to 1. (Can only select one).

skewed = 0;

rho = 10;
a = 0.1;
b = 0.9;

u = a+b;
v = b/((a+b)^2);

epsi = sqrt(0.001);
L = sqrt(0.2); % domain size 
Du = epsi^2/L^2;
Dv = 1/(L^2);

scaling_factor = 0.99;
Psi_c = @(mu)Psi(mu,w,F);
phi = @(x)1/2*(1+erf(x/sqrt(2)));
F = @(x) phi(x) - (2*owens_T(x,rho));


tauvec = [0.1:0.2:1];
results_vec = zeros(length(tauvec),1);
for j = 1:length(tauvec)

tau = tauvec(j);

sig = @(tau)tau/3 * scaling_factor;
w = @(mu)mu/3 * scaling_factor;

rhoh = (1+rho^2)^(1/2);
phi = @(x)1/2*(1+erf(x/sqrt(2)));
F = @(x) phi(x) - (2*owens_T(x,rho));
skew_tau1 = @(mu)mu - (3*w(mu));
skew_tau2 = @(mu)mu + (3*w(mu));
Psi_c = @(mu)Psi(mu,wmax,F);
Kskew = @(s,mu)Psi_c(mu)/w(mu) * sqrt(2/pi) * exp((-1/2)*((s-mu)/w(mu)).^2) .* phi(rho*(s-mu)/w(mu));
func = @(mu) -tau + mu + ((w(mu)*Psi_c(mu))*(Kskew(skew_tau1(mu),mu)-Kskew(skew_tau2(mu),mu) + ((2*rho/(rhoh*sqrt(2*pi)))*(phi(rhoh*(skew_tau2(mu)-mu)/w(mu))-phi(rhoh*(skew_tau1(mu)-mu)/w(mu))))));


mu = fsolve(func,tau);
wmax = w(mu);


    kvec = [0:0.2:20];
    res = zeros(length(kvec),1);
    tau1 = mu+(3*wmax);
    tau2 = mu-(3*wmax);
    int_vec = linspace(tau1,tau2,1000);
    E = @(x)trapz(integrand(x,mu,wmax,int_vec,Psi_c,rho)); 

    for i = 1:length(kvec)
    k = kvec(i);
    [ak, bk, gk, dk, chik] = get_coeffs(k, Du, Dv, u, v, a, b, 1, 0, 0);
    res(i) = max(DispersRel_dist(tau, ak, bk, gk, dk, chik, E)); % storing lambda(k)
    end
results_vec(j) = max(res);


end

function y = integrand(x,mu,wmax,s,Psi_c,rho)

spacing = diff(s);
spacing = abs(spacing(1));
y = Psi_c(mu)/(wmax*sqrt(2*pi)) * (1+erf(rho*(s-mu)/wmax*sqrt(2))) .* exp(-0.5 * ((s-mu)/wmax).^2 - x*s);
y = y*spacing;
end

function y = Psi(mu,wmax,F)

tau1 = mu - (3*wmax);
tau2 = mu + (3*wmax);

 y = 1/(F((tau2-mu)/wmax)-F((tau1-mu)/wmax));

end

function T = owens_T(x,rho)

num = 10000;
s = linspace(0,rho,num);
ds = rho/(num-1);
integrand = exp(-1/2*(x^2)*(1+s.^2))./(1+s.^2);
integral = trapz(integrand);

T = ds * integral / (2*pi);


end