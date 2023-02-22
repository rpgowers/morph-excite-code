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

@load "data/MLDS-bif-general.bson" gin Isn Ih Sh gbt Ibt gc Ic τδ

gon = [3.0, 4.0, 4.7, 4.8, 5.0, 5.4, 5.5, 5.9]
Ion = zeros(length(τδ), length(gon))
ron = zeros(length(τδ), length(gon))
idx = findfirst.(isequal.(gon), (gin,))

ghmax = zeros(length(τδ))
for i in eachindex(τδ)
	hmax = findfirst(isnan.(Ih[i,:,2]))
	ghmax[i] = gin[hmax]
end

args_test = MLMDS_Param(gL=2.0, M=M, λ=λ, L=L)
ρon = discrete_rho.(gon, Ref(args_test))

v0 = -60.0.*ones(M+1)
n0 = 0.0
x0 = vcat(n0, v0)

for j in eachindex(τδ)
  println("τδ = $(τδ[j]) ms")
  for i in eachindex(gon)
    println("gin = $(gon[i]) nS")
    args = MLMDS_Param(gL=2.0, ρ=ρon[i], τδ = τδ[j], M=M, λ=λ, L=L)
    if gon[i] < gbt[j]
      # println(Isn[idx[i],:])
      Ibif = maximum(Isn[idx[i],:])
      Ibound = [Ibif-0.5, Ibif+0.5]
      

      @time Ion[j,i], ron[j,i] = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e6, saveat=[])
    elseif gon[i] < ghmax[j]
      Ibif = minimum(Ih[j,idx[i],:])
      Ibound = [Ibif-5.0, Ibif+1.0]
      @time Ion[j,i], ron[j,i] = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, ISI_min = 9, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=0.9, solver=Tsit5(), maxiters=1e6, saveat=[])
    else
    end

  end
end

println(Ion)
println(ron)

@save "data/MLDS-onset.bson" gon ρon ron Ion M λ L τδ date tspan cutoff