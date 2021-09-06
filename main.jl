
using Random, Distributions
using LinearAlgebra, Plots
using SparseArrays, DifferentialEquations

include("func_mod.jl")
##
prob_type = 4 # 1 = No Delay, 2 = Fixed Delay, 3 = Dist. Delay, 4 = Skewed Dist. Delay


tlen = 100 # Time span of simulations
# Boundary conditions
bc_u = "Neumann"      # Or "Dirichlet"
bc_v = "Neumann"
m = 500 # Number of spatial discretisation points
L = 30 # Domain size
λ = 10 # Skew parameter : ρ in dissertation
# model parameters
a = 0.1
b = 0.9
tau = 1 # Time delay for problem types 2-4. If problem type 4: τ here corresponds to μ in dissertation,
        # namely the location (not mean of distribution)
num_sd = 3  # n of std deviations
σmax = tau/num_sd
σ = σmax * 0.99   # std for problem types 3-4. If problem type 4: σ here corresponds to ω in dissertation,
                  # namely the scale (not std of distribution)

Random.seed!(0)  # Setting a random seed

if prob_type == 1
    include("no_delay_funcs.jl")
    funcs = no_delay_Funcs
    f_kin = (a,b,u,v) -> a .-u .+ (u.^2 .* v)  # Schnakenberg Kinetics of f
    global p = Funcs.param_init(f = f_kin, tspan=[0, tlen], bc_u=bc_u, bc_v=bc_v, a=a, b=b, L=L)

elseif prob_type == 2

    include("fixed_delay_funcs.jl")
    funcs = fixed_delay_Funcs
    f_kin = (a,b,u,v,ulag,vlag) -> a .-u .- 2(u.^2 .* v) .+ 3(ulag.^2 .* vlag) # Schnakenberg kinetics of f
    f_gm1 = (a,b,u,v,ulag,vlag) ->a .- (b.*u) .+ ((ulag.^2)./vlag) # GM 1st variant for f
    f_gm2 = (a,b,u,v,ulag,vlag) ->a .- (b.*u) .+ ((ulag.^2)./v) # GM 2nd variant for f
    g_gm = (a,b,u,v,ulag,vlag) -> ulag.^2 .- v # GM kinetics for g
    global p = Funcs.param_init(f = f_kin, tspan=[0, tlen], bc_u=bc_u, bc_v=bc_v,tau=tau, a=a, b=b, L=L)

elseif prob_type == 3

    include("dist_delay_funcs.jl")
    funcs = dist_delay_Funcs
    global p = Funcs.param_init(tspan=[0, tlen], bc_u=bc_u, bc_v=bc_v, tau=tau,
                         σ=σ, num_sd=num_sd, a=a, b=b)

elseif prob_type == 4
    include("dist_skew_delay_funcs.jl")
    funcs = dist_skew_delay_Funcs
    global p = Funcs.param_init(tspan=[0, tlen], bc_u=bc_u, bc_v=bc_v, tau=tau,
                         σ=σ, num_sd=num_sd, λ = λ, a=a, b=b)

end


In = Diagonal(ones(Int(p.m)))
jac_sparsity = sparse([p.A_u+In In; In p.A_u+In])  # Setting up a sparsity matrix of Jacobian for numerical solver

prob = funcs.run_deriv(p, jac_sparsity)  # Defining Julia DE problem to solve

# Run solver
stepsize = tlen/10000   # Discretisation of solution points to save
sol = solve(prob, alg_hints=[:stiff], reltol=1e-6, abstol=1e-6, saveat = 0:stepsize:tlen, dtmax = 0.1)

# Format solution for plotting
u,v = Funcs.format_sol(sol,p.m)
Funcs.turing_plot(u,v,p,prob_type,sol)  # Plot
