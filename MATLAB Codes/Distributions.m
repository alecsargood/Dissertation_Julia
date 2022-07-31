%% Run this section to solve eq12 with given tau and rho parameters.

clear all;
tau = 0.1;
rho = 10;
scaling_factor = 0.99; % vary this 0.99, 0.2, 0.01

sig = @(tau)tau/3 * scaling_factor;
w = @(mu)mu/3 * scaling_factor;

rhoh = (1+rho^2)^(1/2);
phi = @(x)1/2*(1+erf(x/sqrt(2)));
F = @(x) phi(x) - (2*owens_T(x,rho));
skew_tau1 = @(mu)mu - (3*w(mu));
skew_tau2 = @(mu)mu + (3*w(mu));
norm_tau1 = @(tau)tau - (3*sig(tau));
norm_tau2 = @(tau)tau + (3*sig(tau));
Psi_c = @(mu)Psi(mu,w,F);
Phi_c = @(tau)Phi(tau,sig,phi);
Kskew = @(s,mu)Psi_c(mu)/w(mu) * sqrt(2/pi) * exp((-1/2)*((s-mu)/w(mu)).^2) .* phi(rho*(s-mu)/w(mu));
Knorm = @(s,tau)(Phi_c(tau)/(sig(tau) * sqrt(2*pi))) * exp((-1/2)*((s-tau)/sig(tau)).^2);
func = @(mu) -tau + mu + ((w(mu)*Psi_c(mu))*(Kskew(skew_tau1(mu),mu)-Kskew(skew_tau2(mu),mu) + ((2*rho/(rhoh*sqrt(2*pi)))*(phi(rhoh*(skew_tau2(mu)-mu)/w(mu))-phi(rhoh*(skew_tau1(mu)-mu)/w(mu))))));


mu = fsolve(func,tau);
wmax = w(mu);
%% Plots skewed distribution (with given rho above)
st1 = skew_tau1(mu);
st2 = skew_tau2(mu);
skew_vec = linspace(st1,st2,200);

skew_vals = Kskew(skew_vec,mu);
plot(skew_vec,skew_vals,'LineWidth',5)
%% Plots symmetric distribution
sigmax = tau/3;
nt1 = norm_tau1(tau);
nt2 = norm_tau2(tau);
norm_vec = linspace(nt1,nt2,200);

dist_vals = Knorm(norm_vec,tau);

plot(norm_vec,dist_vals,'LineWidth',5)
xlabel('s')
ylabel('pdf')
legend('Symmetric Normal', 'Positively Skewed', 'Negatively Skewed')
hold on
%%

function y = Phi(tau,sig,phi)

tau1 = tau - (3*sig(tau));
tau2 = tau + (3*sig(tau));

 y = 1/(phi((tau2-tau)/sig(tau))-phi((tau1-tau)/sig(tau)));


end

function y = Psi(mu,w,F)

tau1 = mu - (3*w(mu));
tau2 = mu + (3*w(mu));

 y = 1/(F((tau2-mu)/w(mu))-F((tau1-mu)/w(mu)));

end


function T = owens_T(x,rho)

num = 10000;
s = linspace(0,rho,num);
ds = rho/(num-1);
integrand = exp(-1/2*(x^2)*(1+s.^2))./(1+s.^2);
integral = trapz(integrand);

T = ds * integral / (2*pi);


end


