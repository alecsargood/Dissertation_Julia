module no_delay_Funcs

using Parameters, Random, DifferentialEquations
using Distributions, LinearAlgebra




function run_deriv(p, jac_sparsity)
    func = ODEFunction(deriv, jac_prototype=jac_sparsity)
    prob = ODEProblem(func, p.ic, p.tspan, p)

    return prob

end

function deriv(du,U,p,t)

    u = U[1:p.m]
    v = U[p.m+1:end]
    if p.m == 1
        du .= [p.f(p.a,p.b,u,v); p.g(u,v)]
    else
        du .= [(p.d_u*p.A_u*u)+p.f(p.a,p.b,u,v); (p.d_v*p.A_v*v)+p.g(p.a,p.v,u,v)]
    end

    return du
end




end
