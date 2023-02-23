using OrdinaryDiffEq, Parameters, Dates, PyPlot, BSON
using BSON: @save, @load

include("sim-functions.jl")

date = string(Dates.today())

rfind = 1e-3
tmax = 10000.0
cutoff = 2000.0
tspan = (0, tmax+cutoff)

M = 50
L = 1000.0
λ = 100.0

@load "data/MLDS-bif-bt.bson" τδ gbt Ibt gc ρc Ic τbtc
τsel = 0.5:0.5:12.5
idx = findfirst.(isequal.(τsel), (τδ,))
niter = 10

gon = zeros(length(τsel))
ρon = zeros(length(τsel))
ron = zeros(length(τsel))
Ion = zeros(length(τsel))

M = 50
	
v0 = -60.0.*ones(M+1)
n0 = 0.0
x0 = vcat(n0, v0)

args_test = MLMDS_Param(gL=2.0, M=M, λ=λ, L=L)

for j in eachindex(τsel)
  println("τδ = $(τsel[j]) ms")
  glow = gbt[idx[j]]-1.0
  gup = gbt[idx[j]]
  ρlow = discrete_rho(glow, args_test)

  args = MLMDS_Param(gL=2.0, ρ=ρlow, τδ = τsel[j], M=M, λ=λ, L=L)
  local vsn, Isn = sn(args)
  local Isn_max = maximum(Isn)
  Ibound = [Isn_max-0.5, Isn_max+0.5]

  @time Ilow, rlow = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e7, saveat=[])
  if Ilow > Isn_max
    println("good lower bound")
  end

  ρup = discrete_rho(gup, args_test)
  args = MLMDS_Param(gL=2.0, ρ=ρup, τδ = τsel[j], M=M, λ=λ, L=L)
  local vsn, Isn = sn(args)
  local Isn_max = maximum(Isn)
  Ibound = [Isn_max-0.5, Isn_max+0.5]

  @time Iup, rup = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e7, saveat=[])
  # println("Iup = $(Iup) pA")
  if Iup < Isn_max
    println("good upper bound")
  end
  println([glow, gup])

  @time ρon[j], Ion[j], ron[j] = bsnl_extract(:ρ, ρup, ρlow, niter, args; maxiters=1e7)
  gon[j] = discrete_gin(MLMDS_Param(M=M, λ=λ, ρ=ρon[j]))

  
end

plot(τsel, gon)
show()

@save "data/MLDS-onset-bsnl.bson" gon ρon ron Ion M λ L τsel date tspan cutoff