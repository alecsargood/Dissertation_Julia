module Funcs

using LinearAlgebra, Random, DifferentialEquations
using Plots, SparseArrays, Distributions
using LaTeXStrings, Polynomials, FFTW
using SpecialPolynomials, SpecialFunctions, Parameters

export create_lap
export turing_plot
export format_sol
export plot_dist
export find_ttp
export fourier_coeff
export owenT


##
# Setting up structure for parameters
@with_kw struct param_init
    a = 0.1; b = 0.9; ϵ = sqrt(0.001);
    γ = 0.05; m = Int(500)
    L = 30; Δx = L/(m-1); xvec = (0:Δx:L)/L
    num_sd = 3

    d_u = ϵ^2 / γ
    d_v = 1 / γ

    bc_u = "Neumann"
    bc_v = "Neumann"
    # Create Laplacian dependent on BCs
    A_u = Funcs.create_lap(Δx, m, bc_u)
    A_v = Funcs.create_lap(Δx, m, bc_v)

    # Random perturbation for IC
    μ = 0
    σ_ic = 0.01
    d = Normal(μ,σ_ic)
    us = a + b
    vs = b/((a+b)^2)
    ic = [us .* (1 .+ rand(d,m)); vs.* (1 .+ rand(d,m))]

## Initial conditions from 2006 Gaffney Monk
    #Aa = -170 - (2/3)
    #Ab = -130.8 - (0.4/9)
    #Ba = 412 + (4/9)
    #Bb = 337.0 + (0.2/3)
    #Ca = -312 - (8/9)
    #Cb = -281.6
    #Da = 71 + (1/9)
    #Db = 75.3 + (0.7/9)
    #Ea = 0.00125
    #Eb = 0.006

    #u0_func = (xvec) -> 1 .- (Eb .* xvec.^7 .* (1 .- xvec.^2) .* (Ab.*xvec.^3 .+ Bb.*xvec.^2 .+ Cb.*xvec .+ Db))
    #v0_func = (xvec) -> 0.9 .+ (Ea .* xvec.^5 .* (1 .- xvec.^2) .* (Aa.*xvec.^3 .+ Ba.*xvec.^2 .+ Ca.*xvec .+ Da))
    #y0_func = (xvec) -> [u0_func(xvec);v0_func(xvec)]

    # ic = y0_func(1 .- xvec)
## Default parameter vals
    tau = 1e-6
    λ = -2
    δ = (1+λ^2)^(1/2)
    σ = 0.1
    tspan = [0,100]
    # Functions for distribution
    ϕ = (x) ->  (1/(sqrt(2*pi))).*exp.(-(x.^2)./2)
    Φ = (x) -> (1/2)*(1 .+ erf.(x./(sqrt(2))))
    # Integration limits
    int_a = tau - (num_sd*σ)
    int_b = tau + (num_sd*σ)
## Derivation of mean and std of skew truncated Gaussian
    u = (int_a - tau) / σ
    v = (int_b - tau) / σ
    k = (s,τ,σ) -> (1/σ).*ϕ.((s.-τ)./σ)./(Φ.((int_b-τ)./σ).-Φ.((int_a-τ)./σ))
    Fskew = (x,τ,σ,λ) -> Φ((x-τ)/σ) - (2*owensT((x-τ)/σ, λ))
    kskew = (s,τ,σ,λ) -> 2/σ .* ϕ.((s-τ)./σ) .* Φ.(λ .* (s-τ) ./ σ)
    skew_fact = Fskew(v,0,1,λ) - Fskew(u,0,1,λ)
    r1 = (kskew(u,0,1,λ)-kskew(v,0,1,λ))/skew_fact + ((2/sqrt(2*pi))*(λ/δ)*(Φ(v*δ)-Φ(u*δ))/skew_fact)
    skew_μ = tau + (σ*r1)
    m1 = (ϕ(v)-ϕ(u))/(Φ(v)-Φ(u))
    r2 = ((u*kskew(u,0,1,λ))-(v*kskew(v,0,1,λ)))/skew_fact + ((2/sqrt(2*pi))*(λ/δ)*(Φ(v*δ)-Φ(u*δ))*m1/skew_fact)
    E2 = tau^2 +(2*σ*tau*r1) + (σ^2*(1+r2))
    skew_var = E2 - (skew_μ^2)
    skew_σ = sqrt(skew_var)
## Qaudrature params
    N = 50
    Δt = (int_b-int_a)/(N)
    lags = int_a:Δt:int_b
## Default kinetics for Schnakenberg
    f = (u,v,int) -> a.-u.-((2*u.^2).*v).+3int
    g = (u,v) -> b.-((u.^2).*v)

end


## Function to create Laplacian based on Finite Difference scheme. bc = "Neumann" or "Dirichlet"
# Arguments: Δx - spatial discresisation
#            m - number of spatial points
#            bc - Boundary conditions type

function create_lap(Δx,m,bc)
    if m != 1
    A = Tridiagonal(ones(m-1),-2*ones(m),ones(m-1))
if bc == "Neumann"
    A[1,2]=2
    A[end,end-1]=2
elseif bc == "Dirichlet"
    A[1,1] = 0
    A[1,2] = 0
    A[end,end] = 0
    A[end,end-1] = 0
end

    A = A/((Δx)^2)
    A = sparse(A)
    else
    A = I
    end
end

## Function to format numerical solution for plotting
#  Arguments: sol - numerical solution structure outputted from solver
#             m - number of spatial points
function format_sol(sol,m)

    u = zeros(m,length(sol))
    v = zeros(m,length(sol))
    for i in 1:length(sol)
        u[:,i] = sol.u[i][1:m]
        v[:,i] = sol.u[i][m+1:end]
    end

    return u,v
end

## Function to plot a heatmap of activator, u (or v for phaseplane with no diffusion)
#  Arguments: u - activator soln outputted from format_sol
#             v - inhibitor soln outputted from format_sol
#             p - parameter structure
#             prob_type - Type of problem 1,2,3,4
#             sol - Soln structure outputted from solver

function turing_plot(u,v,p,prob_type,sol)
    tlen = p.tspan[2]
    if p.m != 1
        xlen = size(u,2)
        ylen = size(u,1)

        if prob_type == 1
            title = "No Delay"
            tau = 0
            σ = 0
        elseif prob_type == 2
            title = "Fixed Delay"
            tau= p.tau
            σ = 0
        elseif prob_type == 3
            title = "Dist Delay"
            tau = p.tau
            σ = round(p.σ,digits=3)
    elseif prob_type == 4
            title = "Skewed Dist Delay"
            tau = p.tau
            σ = round(p.σ,digits=3)
            λ = p.λ
        end

        Plots.display(heatmap(u, xticks=([0,xlen/4,xlen/2,xlen*3/4,xlen], [0, tlen/4, tlen/2, tlen*3/4, tlen]),
                  yticks=([0,ylen/4,ylen/2,ylen*3/4,ylen], [0, 1/4, 1/2, 3/4, 1]),
                  xlabel ="Time", ylabel="X", title = "BC = $(p.bc_u), type = $(title), τ=$(tau), σ=$(σ)"))
    else
        u = reshape(u,length(u))
        v = reshape(v,length(v))
        phaseplane = (plot(u,v,xlabel="u",ylabel="v",linewidth=5))
        timeseries = (plot(sol.t, u, xlabel="Time", ylabel="u", linewidth=5))
        return phaseplane, timeseries
    end
end

## Function to plot distribution (for prob_type = 3,4)
#  Arguments: p - parameter structure
#             prob_type - type of problem 3,4
# Currently set up to plot skewed truncated Gaussian
function plot_dist(p, prob_type)

    N = 1000
    τ = p.tau
    σ = p.σ
    Δx = (p.int_b-p.int_a)/(N-1)
    xvec = p.int_a:Δx:p.int_b


    if prob_type == 3
    plot!(xvec, p.k.(xvec,τ,σ), xlabel="τ", ylabel="pdf", label = "σ=$(round(σ,digits=2))",linewidth=5)
    vline!([τ-σ,τ+σ], label = "1 σ")
    vline!([τ-(2*σ),τ+(2*σ)], label = "2 σ")
    vline!([τ-(3*σ),τ+(3*σ)], label = "3 σ")
    else
    factor = (p.Fskew(p.int_b,p.tau,p.σ,p.λ)-p.Fskew(p.int_a,p.tau,p.σ,p.λ))
    plot!(xvec, p.kskew.(xvec,τ,σ,p.λ) ./ factor, xlabel="τ", ylabel="pdf", label = "ω=$(round(σ,digits=2))", linewidth=5)
    #vline!([p.tau], label="μ",linewidth=5)
    #vline!([p.skew_μ], label = "τ", linewidth=5)
    end

end

## Function to find time-to-pattern
#  Arguments: u - activator soln outputted from format_sol
#             thresh - threshold value to use for pattern
#             stepsize - stepsize used for solution points saved
#             us - steady state of activator

function find_ttp(u, thresh, stepsize, us)

len = size(u,2)
val = 0
i = 0
while val < thresh

i = i+1
val = maximum(abs.(u[:,i] .- us))

end

t_taken = stepsize*i

return val, i, t_taken

end


## Function to compute Fourier coefficients based on FFT
#  Arguments: u - activator soln outputted from format_sol
#             m - number of spatial discretisation points

function fourier_coeff(u,m)

    Fy = fft(u)[1:Int(m/2)]
    ak =  2/m * real.(Fy)
    bk = -2/m * imag.(Fy)


    return ak,bk
end


## Function to evaluate Owen's T integral
# Composite Simpson's rule applied
#  Arguments: x,y - arguments of Owen's T function T(x,y)

function owensT(x,y)

    int = 0
    N = 100000
    h = y / N
    points = 0:h:y

    for i = 2:length(points)-1
        s = points[i]
        val = exp(((-(x^2))/2)*(1+(s^2))) / (1+s^2)
        if i % 2 == 0
            int += 2*val
        else
            int += 4*val
        end
    end

    int += (exp(((-(x^2))/2)*(1+(points[1]^2))) / (1+points[1]^2)) + (exp(((-(x^2))/2)*(1+(points[end]^2))) / (1+points[end]^2))

    int = (int * h / 3) /(2*pi)

    return int
end

## Function to evaluate mean of distribution   ∫s k(s) ds
#  Composite Simpson's rule applied
#  Arguments: pdf - function handle of pdf being used, which takes one arguments, such that pdf(s)
#             a,b - lower and upper integration limits

function compute_mean(pdf,a,b)

    int = 0
    N = 100000
    h = (b-a) / N
    points = a:h:b

    for i = 2:length(points)-1
        s = points[i]
        val = pdf(s).*s
        if i % 2 == 0
            int += 2*val
        else
            int += 4*val
        end
    end

    int += (pdf(points[1])*points[1]) + (pdf(points[end])*points[end])

    int = (int * h / 3)

    return int
end


end
