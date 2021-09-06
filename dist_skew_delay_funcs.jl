module dist_skew_delay_Funcs

using Parameters, Random, DifferentialEquations
using Distributions, LinearAlgebra


function run_deriv(p, jac_sparsity)

    h(p, t) = p.ic
    func = DDEFunction(deriv, jac_prototype=jac_sparsity)
    prob = DDEProblem(func, p.ic, h, p.tspan, p, constant_lags = p.lags)

    return prob

end

function deriv(du,U,h,p,t)

    u = U[1:p.m]; v = U[p.m+1:end]
    int = comp_simp(p, h, t)
    du .= [(p.d_u*p.A_u*u)+p.f(u,v,int); (p.d_v*p.A_v*v)+p.g(u,v)]

    return du
end

function comp_simp(p, hist, t) # Implementation of composite Simpson's rule to evaluate integral term in model


factor = (p.Fskew(p.int_b,p.tau,p.σ,p.λ)-p.Fskew(p.int_a,p.tau,p.σ,p.λ))
int = zeros(p.m)
 for i = 2:length(p.lags)-1
     ulag = hist(p, t-p.lags[i])[1:p.m]

     vlag = hist(p, t-p.lags[i])[p.m+1:end]
     pdf = p.kskew(p.lags[i], p.tau, p.σ, p.λ)
     if i % 2 == 0
     int += 2*pdf.* ulag.^2 .* vlag
     else
     int += 4*pdf.* ulag.^2 .* vlag
     end
 end
     ulag1 = hist(p, t-p.lags[1])[1:p.m]
     vlag1 = hist(p, t-p.lags[1])[p.m+1:end]

     pdf1 = p.kskew(p.lags[1], p.tau, p.σ, p.λ)
     ulag2 = hist(p, t-p.lags[end])[1:p.m]
     vlag2 = hist(p, t-p.lags[end])[p.m+1:end]

     pdf2 = p.kskew(p.lags[end], p.tau, p.σ, p.λ)

     int += (pdf1.* ulag1.^2 .* vlag1) + (pdf2.* ulag2.^2 .* vlag2)

     int = int * ((p.int_b-p.int_a)/(p.N)) / 3

     int = int ./ factor
 return int

end

end
