using OrdinaryDiffEq, Parameters, PyPlot
using BSON: @save, @load

include("sim-functions.jl")

@load "data/MLS-onset.bson" gon Ion

	tspan = (0.0, 11000.0)
	v0 = -60.0
	n0 = 0.0
	x0 = vcat(n0, v0)

	δθ = 0.01
	θin = 0.5*δθ:δθ:1-0.5*δθ
	gsel = [3.0, 4.0, 4.5, 4.7, 4.8, 5.0, 5.4, 5.5, 5.9]
	idx = findfirst.(isequal.(gsel), (gon,))
	Δθsim = zeros(length(gsel), length(θin))

	ΔV = [1.0e-2, 1.0e-2, 2.0e-7, 5.0e-7, 5.0e-7, 1.0e-6, 1.0e-4, 1.0e-5, 2.0e-3]

for i in eachindex(gsel)
	println(gsel[i])
	args = MLS_Param(gL = gsel[i], Iext = Ion[idx[i]])

	@time outref,tref,iref = LC_extract(tspan, x0, args; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
	xp = outref[end-2]
	T0 = 0.5*(diff(tref)[end-1]+diff(tref)[end])
	@time Δθsim[i,:] = PRC_sim(T0, xp, args, ΔV[i], θin; Ps = 10, solver=Tsit5(), soma_idx=2)
end

fig = figure(figsize=(9,9))
fig.subplots_adjust(hspace = 0.3, wspace=0.3)

for i in eachindex(gsel)
	fig.add_subplot(330+i)
	title("\$g_{\\mathrm{in}}\$ = $(gsel[i])")
	plot(θin, Δθsim[i,:], ".")
end

@save "data/MLS-prcs.bson" gsel ΔV Δθsim θin

show()