%% Select model

% Set desired model to 1. (Can only select one).

LI = 1;
GM1 = 0;
GM2 = 0;

if LI == 1

u_star = @(a,b)a+b;
v_star = @(a,b)b/((a+b)^2);

avec = 1e-6:0.05:0.8;
bvec = 1e-6:0.05:1;

else

u_star = @(a,b)(a+1)/b;
v_star = @(a,b)((a+1)/b)^2;

avec = 1e-6:0.02:1;
bvec = 1e-6:0.02:4;

end

%% Bifurcation Diagram 

tau = 0;  % Set Delay
kmax = 10;

epsi = sqrt(0.001);
L = sqrt(4.5); 
Du = epsi^2/L^2;
Dv = 1/(L^2);

len_a = length(avec);
len_b = length(bvec);
mat = zeros(len_b,len_a);
% steady states

    parfor i = 1:len_a
        a = avec(i);
        for j = 1:len_b
            b = bvec(j);
              u = u_star(a,b);
              v = v_star(a,b);
            res = zeros(kmax+1,1);
for k = 0:kmax
    [ak, bk, gk, dk, chik] = get_coeffs(k, Du, Dv, u, v,a,b, LI, GM1, GM2);
    res(k+1) = max(DispersRel(tau, ak, bk, gk, dk, chik)); % storing lambda(k)
end
    val = max(res);
    mat(j,i) = val;
    zmat(j,i) = res(1);
        end
    end
    
imagesc(flip(mat))
if Schnak == 1
clim([-1, 1])
else
clim([0, 4])
end
hold on
contour(flip(mat),[0,0],'k','LineWidth',5)
contour(flip(zmat),[0,0],'k','LineWidth',5)

%% Dispersion relations
 
kmax = 5;
res = zeros(kmax+1,1);
a = 0.1;
b = 0.9;
tau = 0.1;
epsi = sqrt(0.001);
L = sqrt(0.2); % domain size 

% diffusive coefficients
Du = (epsi^2)/L^2;
Dv = 1/(L^2);
% steady states
u = u_star(a,b);
v = v_star(a,b);
kvec = 0:0.1:kmax;
for i = 1:length(kvec)
    k = kvec(i);
    [ak, bk, gk, dk, chik] = get_coeffs(k, Du, Dv, u, v, a, b, LI, GM1, GM2);
    res(i) = max(DispersRel(tau, ak, bk, gk, dk, chik)); % storing lambda(k)
end

    
