using OrdinaryDiffEq, Dates, PyPlot, BSON, JSON
using BSON: @save, @load

include("sim-functions.jl")

date = string(Dates.today())

rfind = 1e-3
tmax = 10000.0
cutoff = 2000.0
tspan = (0, tmax+cutoff)
τδ = 2.5

num = "00892"
cell = "NMO_$(num)"
depth_sel = ["depth8", "depth13", "depth14", "full"]

open("data/$(cell)_cylinder_parameters.json","r") do f
  global param_dict
  param_dict = JSON.parse(f)
end

l = param_dict["leff"]
d = param_dict["deff"]
ρ∞ = param_dict["rhoinf"]
depths = param_dict["depths"]
idx = findfirst.(isequal.(depth_sel), (depths,))
lsel = l[idx]

Ion = zeros(length(depth_sel))
ron = zeros(length(depth_sel))
ρon = zeros(length(depth_sel))

gσ = 2.0
gon = gσ.*(ρ∞[idx].+1)
Mmax = 50
Lmax = 1000

λ = Lmax/lsel[end]
L = lsel.*λ
M = Int.(round.(L.*Mmax./Lmax))

for i=1:1 # in eachindex(depth_sel)
	v0 = -60.0.*ones(M[i]+1)
	n0 = 0.0
	x0 = vcat(n0, v0)

	println("gin = $(gon[i]) nS")

	vbt, Ibt, ρbt = bt(MLFDS_Param(gL=gσ, λ=λ, L=L[i]))
	gbt = gσ*(1+ρbt[1]*tanh(L[i]/λ))
	println(gbt)

	args_test = MLMDS_Param(gL=gσ, M=M[i], λ=λ, L=L[i], τδ=τδ)
	ρon[i] = discrete_rho(gon[i], args_test)

	args = MLMDS_Param(gL=gσ, ρ=ρon[i], τδ = τδ, M=M[i], λ=λ, L=L[i])

	if gon[i] < gbt

		vsn, Isn = sn(args)
		Ibound = [maximum(Isn)-0.5, maximum(Isn)+0.5]

		@time Ion[i], ron[i] = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=2, ϵr = 0.1, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e8, saveat=[])
		println(Ion[i])
		println(ron[i])

	else
		_, idx_l = findmin(abs.(lsel[i].-l_bif))
		_, idx_g = findmin(abs.(gon[i].-gin_bif))
		println(l_bif[idx_l])
		println(gin_bif[idx_g])
		Ibound = [minimum(Ih_bif[idx_l, idx_g, :])-2.0, minimum(Ih_bif[idx_l, idx_g, :])+1.0]

		@time Ion[i], ron[i] = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, ISI_min = 9, vth=-8.0, Nmax=1, ϵr = 0.1, ϵisi = 10, ϵsp=0.9, solver=Tsit5(), maxiters=1e8, saveat=[])
		println(Ion[i])
		println(ron[i])

	end

	probf = ODEProblem(neuron_sim!, x0, tspan, MLMDS_Param(gL=2.0, ρ=ρon[i], τδ = τδ, M=M[i], λ=λ, L=L[i], Iext = Ion[i]))
	@time solf = solve(probf, Tsit5(), reltol=1e-8, abstol=1e-8, maxiters=1e8, save_idxs=2)
	vf = solf[1,:]
	t = solf.t
	plot(t, vf)
	show()

end