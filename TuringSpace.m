%% Schnakenberg model

avec = 1e-6:0.005:1.4;
bvec = 1e-6:0.005:2;
u_star = @(a,b) a+b;
v_star = @(a,b) b/((a+b)^2);
f_u = @(a,b) -1+(2*u_star(a,b)*v_star(a,b));
f_v = @(a,b) u_star(a,b)^2;
g_u = @(a,b) -2*u_star(a,b)*v_star(a,b);
g_v = @(a,b)-(u_star(a,b)^2);

%% GM model

avec = 1e-6:0.005:1;
bvec = 1e-6:0.005:4;
u_star = @(a,b) (a+1)/b;
v_star = @(a,b) ((a+1)/b)^2;
f_u = @(a,b) -b+(2*u_star(a,b)/v_star(a,b));
f_v = @(a,b) -(u_star(a,b)^2) / (v_star(a,b)^2);
g_u = @(a,b) 2*u_star(a,b);
g_v = @(a,b) -1;


%% Computing Turing space and plot

epsil= sqrt(0.001);
L = sqrt(4.5);
mat = zeros(length(bvec),length(avec));

for i = 1:length(avec)
    a = avec(i);
    for j = 1:length(bvec)
        b = bvec(j);
        fu = f_u(a,b); fv = f_v(a,b); gu = g_u(a,b); gv = g_v(a,b);

        if fu+gv<0 && (fu*gv - fv*gu)>0 && (fu/epsil^2 + gv)>0 && ...
            ((fu/(epsil^2))+gv)^2 -((4/epsil^2)*(fu*gv-fv*gu))>0

            mat(j,i) = 1;
        end
    end
end

imagesc(flip(mat));


