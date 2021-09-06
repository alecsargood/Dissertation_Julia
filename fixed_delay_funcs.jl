module fixed_delay_Funcs

using Parameters, Random, DifferentialEquations
using Distributions, LinearAlgebra

function run_deriv(p, jac_sparsity)

    h(p,t) = p.ic
    func = DDEFunction(deriv, jac_prototype=jac_sparsity)
    prob = DDEProblem(func, p.ic, h, p.tspan, p, constant_lags = p.tau)

    return prob

end

function deriv(du,U,hist,p,t)

    u = U[1:p.m]; v = U[p.m+1:end]

    ulag = hist(p,t-p.tau)[1:p.m]
    vlag = hist(p,t-p.tau)[p.m+1:end]

    if p.m == 1
        du .= [p.f(p.a,p.b,u,v,ulag,vlag); p.g(u,v)]
    else
        du .= [(p.d_u*p.A_u*u)+p.f(p.a,p.b,u,v,ulag,vlag); (p.d_v*p.A_v*v)+p.g(u,v)]
    end
    return du
end

end
