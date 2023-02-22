using OrdinaryDiffEq, Parameters, PyPlot, Dates
using BSON: @save, @load

date = string(Dates.today())

include("sim-functions.jl")

@load "data/MLDS-onset.bson" gon ρon ron Ion M λ L τδ tspan

v0 = -60.0.*ones(M+1)
n0 = 0.0
x0 = vcat(n0, v0)

δθ = 0.01
θin = 0.5*δθ:δθ:1-0.5*δθ

Δθsim = zeros(length(τδ), length(gon), length(θin))

ΔV = [5e-2 2e-2 3e-7 1e-7 1e-7 2e-3 1e-3 8e-2;
	    5e-2 2e-2 5e-3 5e-4 5e-7 1e-3 5e-5 2e-1;
	    5e-2 5e-2 3e-2 2e-2 1e-2 1e-4 2e-5 1e-1;
	    5e-2 5e-2 5e-2 5e-2 5e-2 5e-2 2e-2 2e-1]

for j in eachindex(τδ)
  println("τδ = $(τδ[j]) ms")
	for i in eachindex(gon)
		println("gin = $(gon[i]) nS")
		args = MLMDS_Param(gL=2.0, ρ=ρon[i], τδ = τδ[j], M=M, λ=λ, L=L, Iext=Ion[j,i])

		@time outref,tref,iref = LC_extract(tspan, x0, args; soma_idx=2, verbose=true, vth = -8.0, mincross=10, solver=Tsit5())
		xp = outref[end-2]
		T0 = 0.5*(diff(tref)[end-1]+diff(tref)[end])
		@time Δθsim[j,i,:] = PRC_sim(T0, xp, args, ΔV[j,i], θin; Ps = 10, solver=Tsit5(), soma_idx=2)
	end
end

fig = figure(figsize=(16,16))
fig.subplots_adjust(hspace = 0.3, wspace=0.3)

for i in eachindex(gon)
	fig.add_subplot(330+i)
	for j in eachindex(τδ)
		plot(θin, Δθsim[j,i,:], "-")
	end
end
show()

@save "data/MLDS-prcs.bson" gon ρon τδ ΔV Δθsim θin M λ L date
