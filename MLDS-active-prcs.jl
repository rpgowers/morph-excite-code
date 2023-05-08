using OrdinaryDiffEq, Parameters, PyPlot, Dates
using BSON: @save, @load

date = string(Dates.today())

include("sim-functions.jl")

@load "data/MLDS-active-onset.bson" gon ρon ron Ion M λ L τδ tspan ϵ

x0 = [zeros(M+1) -60.0*ones(M+1)]

δθ = 0.01
θin = 0.5*δθ:δθ:1-0.5*δθ

# Δθsim = zeros(length(ϵ), length(gon), length(θin))
@load "data/MLDS-active-prcs.bson" Δθsim


ΔV = [5e-2 5e-2 5e-2 5e-2 5e-2 5e-3 5e-3 2e-2;
			5e-2 5e-2 5e-2 5e-2 5e-2 5e-3 5e-3 1e-5;
	    5e-2 5e-2 3e-2 2e-2 1e-2 1e-4 2e-5 1e-1;
	    5e-2 5e-2 5e-2 5e-2 5e-2 5e-2 5e-2 5e-2]

for i=2:2 # in eachindex(ϵ)
  println("ϵ = $(ϵ[i])")
	for j=1:7 # in eachindex(gon)
		println("gin = $(gon[j]) nS")
		args = MLADS_Param(gL=2.0, ρ=ρon[j], τδ = τδ, M=M, λ=λ, L=L, Iext=Ion[i,j], ϵ = ϵ[i])

		@time outref,tref,iref = LC_extract(tspan, x0, args; soma_idx=M+2, verbose=true, vth = -8.0, mincross=10, solver=Tsit5())
		xp = outref[end-2]
		T0 = 0.5*(diff(tref)[end-1]+diff(tref)[end])

		# prob = ODEProblem(neuron_sim!, xp, 2*T0, args)
    # @time sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8, dense=false)
    # v = sol[1,end,:]
		# plot(sol.t, v)
		# show()

		@time Δθsim[i,j,:] = PRC_sim(T0, xp, args, ΔV[i,j], θin; Ps = 10, solver=Tsit5(), soma_idx=M+2)
	end
end

fig = figure(figsize=(16,16))
fig.subplots_adjust(hspace = 0.3, wspace=0.3)

for j in eachindex(gon)
	fig.add_subplot(330+j)
	for i in eachindex(ϵ)
		plot(θin, Δθsim[i,j,:], "-")
	end
end
show()

@save "data/MLDS-active-prcs.bson" gon ρon τδ ΔV Δθsim θin M λ L date ϵ
