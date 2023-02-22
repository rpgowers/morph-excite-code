using OrdinaryDiffEq, PyPlot, Parameters
using BSON: @save

include("sim-functions.jl")

M = 50 # number of compartments per neuron

args_s = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 89.0, ρ = 0.8, M=M, λ=100.0)
args_h = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 131.79697, ρ = 1.4, M=M, λ=100.0)

v0 = -20.0 .*ones(M+1)
n0 = 0.0
x0 = vcat(n0, v0)
cutoff = 1000.0
T = 10000.0
tspan = (0.0,T+cutoff)

ΔVs = 6.65e-2
ΔVh = 2.665e-4

δθ = 0.01
θin = 0.5*δθ:δθ:1-0.5*δθ

@time outref,tref,iref = LC_extract(tspan, x0, args_s; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_s = outref[end-2]
T0_s = 0.5*(diff(tref)[end-1]+diff(tref)[end])
@time Δθs = PRC_sim(T0_s, xp_s, args_s, ΔVs, θin; Ps = 10, solver=Tsit5(), soma_idx=2)

@time outref,tref,iref = LC_extract(tspan, x0, args_h; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_h = outref[end-2]
T0_h = 0.5*(diff(tref)[end-1]+diff(tref)[end])
@time Δθh = PRC_sim(T0_h, xp_h, args_h, ΔVh, θin; Ps = 10, solver=Tsit5(), soma_idx=2)

fig = figure(figsize=(9.0, 3.0))
ax1 = fig.add_subplot(121)
plot(θin, Δθs)

ax2 = fig.add_subplot(122)
plot(θin, Δθh)

show()

@save "data/MLDS-exc-prcs.bson" θin Δθs Δθh 