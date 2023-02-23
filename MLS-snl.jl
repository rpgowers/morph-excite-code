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
niter = 10
gup = gbt[1]
glow = gbt[1]-1.0

args_low = MLS_Param(gL=glow)
vsn, Isn = sn(args_low)
Isn_max = maximum(Isn)
Ibound = [Isn_max-0.5, Isn_max+0.5]
@time Ilow, rlow = find_onset(tspan, x0, args_low, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e6, saveat=[])

if Ilow > Isn_max
  println("good lower bound")
end

args_up = MLS_Param(gL=gup)
vsn, Isn = sn(args_up)
Isn_max = maximum(Isn)
Ibound = [Isn_max-0.5, Isn_max+0.5]
@time Iup, rup = find_onset(tspan, x0, args_up, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e6, saveat=[])

if Iup < Isn_max
  println("good upper bound")
end

println([glow, gup])

@time gon, Ion, ron = bsnl_extract(:gL, gup, glow, niter, args_up; maxiters=1e7)

println(gon)

@save "data/MLS-onset-bsnl.bson" gon ron Ion tspan cutoff date