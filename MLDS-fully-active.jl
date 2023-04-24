using OrdinaryDiffEq, PyPlot, Parameters
using BSON: @save

include("sim-functions.jl")

test = ARGS[1]

L = 1000.0
M = 50

if test == "rest"
  ϵ = [0.1, 0.2, 0.5, 1.0]
  v_act = zeros(length(ϵ), M+1)

  Iext = 20.0
  ρ = 1.0
  args_pas = MLMDS_Param(M=M, Iext = Iext, ρ = ρ)
  y = 0:L/M:L

  tspan = (0.0, 1000.0)

  x1 = vcat(0.0, -60.0*ones(M+1))
  prob1 = ODEProblem(neuron_sim!, x1, tspan, args_pas)
  @time sol1 = solve(prob1, Tsit5(), reltol=1e-8, abstol=1e-8, dense=false)
  v_pas = sol1[2:end,end]

  x2 = vcat(zeros(M+1), -60.0*ones(M+1))
  for i in eachindex(ϵ)
    args_act = MLADS_Param(M=M, Iext = Iext, ρ = ρ, ϵ = ϵ[i])
    prob2 = ODEProblem(neuron_sim!, x2, tspan, args_act)
    @time sol2 = solve(prob2, Tsit5(), reltol=1e-8, abstol=1e-8, dense=false)
    v_act[i,:] = sol2[M+2:end,end]
  end

  # println(sum(sol1[2:end,end].-sol2[M+2:end,end]))

  plot(y, v_pas, alpha=0.75, label="Passive, ϵ = 0")
  for i in eachindex(ϵ)
    plot(y, v_act[i,:], "-", alpha=0.75, label="ϵ = $(ϵ[i])")
  end
  legend(frameon=false)
  xlabel("Position (μm)", fontsize=14)
  ylabel("Steady state voltage (mV)")

  # plot(sol1.t, sol1[2,:])
  # plot(sol2.t, sol2[M+2,:], "--")
  # savefig("figures/active-dendrite-rest.svg", bbox_inches="tight")

  show()

elseif test == "onset"
  ϵ = [0.1, 0.2, 0.5, 1.0]
  Ion_act = similar(ϵ)
  ron_act = similar(ϵ)
  ρ = 1.0

  rfind = 1e-3
  tmax = 10000.0
  cutoff = 2000.0
  tspan = (0, tmax+cutoff)

  x1 = vcat(0.0, -60.0*ones(M+1))
  Ibound = [100.0, 110.0]
  args_pas = MLMDS_Param(M=M, ρ = ρ)
  @time Ion_pas, ron_pas = find_onset(tspan, x1, args_pas, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e6, saveat=[])
  println(Ion_pas)
  println(ron_pas)

  x2 = vcat(zeros(M+1), -60.0*ones(M+1))
  # args_act = MLADS_Param(M=M, ρ = ρ, ϵ = ϵ[1], Iext=100.0)
  # @time r, t, vf = rate_measure(tspan, x2, args_act, cutoff; soma_idx=M+2, vth=-8.0, ISImin = 4, ISIfactor=1.1, ϵisi = 10, ϵsp = 1.0, verbose=false, solver=Tsit5(), maxiters=1e6, saveat=[], save_vt = true)
  
  # println(r)
  # plot(t,vf)
  # show()

  for i=1:1 # in eachindex(ϵ)
    args_act = MLADS_Param(M=M, ρ = ρ, ϵ = ϵ[i])
    @time Ion_act[i], ron_act[i] = find_onset(tspan, x2, args_act, cutoff, rfind, Ibound; soma_idx=M+2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e6, saveat=[])   
  end
  println(Ion_act)
  println(ron_act)
end