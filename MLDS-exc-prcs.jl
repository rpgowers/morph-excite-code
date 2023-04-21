using OrdinaryDiffEq, PyPlot, Parameters
using BSON: @save

include("sim-functions.jl")

M = 50 # number of compartments per neuron

args_snic = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 89.0, ρ = 0.8, M=M, λ=100.0)
args_hom = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 131.79697, ρ = 1.4, M=M, λ=100.0)
args_hopf = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 169.13389325141907, ρ = 1.88, M=M, λ=100.0)

v0 = -20.0 .*ones(M+1)
n0 = 0.0
x0 = vcat(n0, v0)
cutoff = 1000.0
T = 10000.0
tspan = (0.0,T+cutoff)

ΔVsnic = 6.65e-2
ΔVhom = 2.665e-4
ΔVhopf = 5.639e-5

δθ = 0.01
θin = 0.5*δθ:δθ:1-0.5*δθ

@time outref,tref,iref = LC_extract(tspan, x0, args_snic; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_snic = outref[end-2]
T0_snic = 0.5*(diff(tref)[end-1]+diff(tref)[end])
@time Δθsnic = PRC_sim(T0_snic, xp_snic, args_snic, ΔVsnic, θin; Ps = 10, solver=Tsit5(), soma_idx=2)

@time outref,tref,iref = LC_extract(tspan, x0, args_hom; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_hom = outref[end-2]
T0_hom = 0.5*(diff(tref)[end-1]+diff(tref)[end])
@time Δθhom = PRC_sim(T0_hom, xp_hom, args_hom, ΔVhom, θin; Ps = 10, solver=Tsit5(), soma_idx=2)

@time outref,tref,iref = LC_extract(tspan, x0, args_hopf; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_hopf = outref[end-2]
T0_hopf = 0.5*(diff(tref)[end-1]+diff(tref)[end])
@time Δθhopf = PRC_sim(T0_hopf, xp_hopf, args_hopf, ΔVhopf, θin; Ps = 10, solver=Tsit5(), soma_idx=2)

fig = figure(figsize=(13.5, 3.0))
ax1 = fig.add_subplot(131)
plot(θin, Δθsnic)

ax2 = fig.add_subplot(132)
plot(θin, Δθhom)

ax3 = fig.add_subplot(133)
plot(θin, Δθhopf)

show()

@save "data/MLDS-exc-prcs.bson" θin Δθsnic Δθhom Δθhopf