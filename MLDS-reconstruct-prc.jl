using OrdinaryDiffEq, Parameters, PyPlot, Dates, BSON, JSON
using BSON: @save, @load

include("sim-functions.jl")

date = string(Dates.today())

@load "data/MLDS-reconstruct-onset.bson" gon ρon ron Ion M λ L τδ date tspan cutoff

δθ = 0.05

θin = 0.5*δθ:δθ:1-0.5*δθ
Δθsim = zeros(length(gon), length(θin))
ΔV = [5e-3, 5e-3, 1e-4, 1e-2]

for i in eachindex(gon)
  println(gon[i])
  v0 = -60.0.*ones(M[i]+1)
  n0 = 0.0
  x0 = vcat(n0, v0)
  args = MLMDS_Param(gL=2.0, ρ=ρon[i], τδ = τδ, M=M[i], λ=λ, L=L[i], Iext=Ion[i])

  @time outref,tref,iref = LC_extract(tspan, x0, args; soma_idx=2, verbose=true, vth = -8.0, mincross=10, solver=Tsit5(), maxiters=1e8, saveat=1e-1)
  xp = outref[end-2]
  T0 = 0.5*(diff(tref)[end-1]+diff(tref)[end])
  @time Δθsim[i,:] = PRC_sim(T0, xp, args, ΔV[i], θin; Ps = 10, solver=Tsit5(), soma_idx=2, maxiters=1e8, saveat=1e-1)
end

@save "data/MLDS-reconstruct-prcs.bson" gon ρon τδ ΔV Δθsim θin M λ L date

fig = figure(figsize=(9,6))
fig.subplots_adjust(hspace = 0.3, wspace=0.3)

for i=1:length(gon)
  fig.add_subplot(230+i)
  plot(θin, Δθsim[i,:], ".")
end
show()