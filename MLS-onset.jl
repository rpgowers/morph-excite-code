using OrdinaryDiffEq, Parameters, Dates, PyPlot, BSON
using BSON: @save, @load

include("sim-functions.jl")

date = string(Dates.today())

rfind = 1e-3
tmax = 10000.0
cutoff = 2000.0
tspan = (0, tmax+cutoff)

v0 = -60.0
n0 = 0.0
x0 = vcat(n0, v0)

@load "data/MLS-bif-general.bson" gin Isn Ih Sh ωh gbt Ibt gc Ic Cbtc

gon = 3.0:0.1:6.0
Ion = zeros(length(gon))
ron = zeros(length(gon))
idx = findfirst.(isequal.(gon), (gin,))

for i in eachindex(gon)

  args = MLS_Param(gL=gon[i])
  if gon[i] < gbt[1]
    Ibif = maximum(Isn[idx[i],:])
    Ibound = [Ibif-1.0, Ibif+1.0]
    Ion[i], ron[i] = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e6, saveat=[])
  else
    Ibif = minimum(Ih[idx[i],:])
    Ibound = [Ibif-5.0, Ibif+1.0]
    @time Ion[i], ron[i] = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, ISI_min = 9,  Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=0.8, solver=Tsit5(), maxiters=1e6, saveat=[])
  end

end

println(Ion)
println(ron)

@save "data/MLS-onset.bson" gon ron Ion tspan cutoff