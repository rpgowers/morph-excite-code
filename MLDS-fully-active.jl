using OrdinaryDiffEq, PyPlot, Parameters, BSON, Dates
using BSON: @save

date = string(Dates.today())

include("sim-functions.jl")

test = ARGS[1]

L = 1000.0
λ = 100.0
M = 50
τδ = 10.0

if test == "benchmark"
  tmax = 2000.0
  cutoff = 0.0
  tspan = (0, tmax+cutoff)

  x2 = [zeros(M+1) -60.0*ones(M+1)]

  args_test = MLADS_Param(M=M, ρ = 5.0, ϵ = 1.0, Iext = 357.0)
  probtest = ODEProblem(neuron_sim!, x2, tspan, args_test)
  @time soltest = solve(probtest, Tsit5(), reltol=1e-8, abstol=1e-8, dense=false, save_idxs = [M+2])

  # @time r, t, vf = rate_measure(tspan, x2, args_test, cutoff; soma_idx=[1,2], vth=-8.0, ISImin = 4, ISIfactor=1.1, ϵisi = 10, ϵsp = 1.0, verbose=false, solver=Tsit5(), maxiters=1e6, saveat=[], save_vt = true)
  
  plot(soltest.t, soltest[1,:])
  show()

elseif test == "rest"
  ϵ = [0.1, 0.2, 0.5, 1.0]
  v_act = zeros(length(ϵ), M+1)

  Iext = 20.0
  ρ = 1.0
  args_pas = MLMDS_Param(M = M, λ = λ, Iext = Iext, ρ = ρ)
  y = 0:L/M:L

  tspan = (0.0, 1000.0)

  x1 = vcat(0.0, -60.0*ones(M+1))
  prob1 = ODEProblem(neuron_sim!, x1, tspan, args_pas)
  @time sol1 = solve(prob1, Tsit5(), reltol=1e-8, abstol=1e-8, dense=false)
  v_pas = sol1[2:end,end]

  x2 = [zeros(M+1) -60.0*ones(M+1)]
  for i in eachindex(ϵ)
    args_act = MLADS_Param(M = M, λ = λ, Iext = Iext, ρ = ρ, ϵ = ϵ[i])
    prob2 = ODEProblem(neuron_sim!, x2, tspan, args_act)
    @time sol2 = solve(prob2, Tsit5(), reltol=1e-8, abstol=1e-8, dense=false)
    v_act[i,:] = sol2[:,end,end]
  end

  plot(y, v_pas, alpha=0.75, label="Passive, ϵ = 0")
  for i in eachindex(ϵ)
    plot(y, v_act[i,:], "-", alpha=0.75, label="ϵ = $(ϵ[i])")
  end
  legend(frameon=false)
  xlabel("Position (μm)", fontsize=14)
  ylabel("Steady state voltage (mV)")

  show()

elseif test == "onset"
  
  # load the onset currents for the passive dendrite model
  data = BSON.load("data/MLDS-onset.bson")
  ρon = data[:ρon]
  gon = data[:gon]
  Ion_pas = data[:Ion]
  
  data_bif = BSON.load("data/MLDS-bif-general.bson")
  gbt = data_bif[:gbt]
  gc = data_bif[:gc]

  ϵ = [0.1, 0.2, 0.5, 1.0]
  # Ion_act = zeros(length(ϵ), length(ρon))
  # ron_act = zeros(length(ϵ), length(ρon))

  rfind = 1e-3
  tmax = 10000.0
  cutoff = 2000.0
  tspan = (0, tmax+cutoff)

  # x1 = vcat(0.0, -60.0*ones(M+1))
  # Ibound = [100.0, 110.0]
  # args_pas = MLMDS_Param(M=M, ρ = ρ)
  # @time Ion_pas, ron_pas = find_onset(tspan, x1, args_pas, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e6, saveat=[])
  # println(Ion_pas)
  # println(ron_pas)

  prev_data = BSON.load("data/MLDS-active-onset.bson")
  ron_act = prev_data[:ron]
  Ion_act = prev_data[:Ion]

  x2 = [zeros(M+1) -60.0*ones(M+1)]
  for i=2:2 # in eachindex(ϵ) 
    
    for j in eachindex(ρon)
      # Ibound = [Ion_pas[3,j]-10.0, Ion_pas[3,j]+10.0]
      Ibound = [Ion_pas[3,j]-30.0, Ion_pas[3,j]+10.0]
      args_act = MLADS_Param(M = M, λ = λ, τδ = τδ, ρ = ρon[j], ϵ = ϵ[i])
      if gon[j] < gc # gbt[3]
        @time Ion_act[i,j], ron_act[i,j] = find_onset(tspan, x2, args_act, cutoff, rfind, Ibound; soma_idx=[M+2], vth=-8.0, Nmax=25, ϵisi = 10, ϵsp=1.0, solver=Tsit5())   
      else
        @time Ion_act[i,j], ron_act[i,j] = find_onset(tspan, x2, args_act, cutoff, rfind, Ibound; soma_idx=[M+2], ISI_min = 9, vth=-8.0, Nmax=25, ϵisi = 10, ϵsp=0.9, solver=Tsit5())
      end

      args_temp = MLADS_Param(M = M, λ = λ, τδ = τδ, ρ = ρon[j], ϵ = ϵ[i], Iext=Ion_act[i,j])

      # @time rtemp, t, vtemp = rate_measure(tspan, x2, args_temp, cutoff; soma_idx=[M+2], vth=-8.0, ISImin = 9, ϵisi = 10, ϵsp = 1.0, verbose=true, solver=Tsit5(), save_vt = true)
      # println(rtemp)
      # plot(t,vtemp)
      # show()
      
    end
  end

  println(Ion_act)
  println(ron_act)

  bson("data/MLDS-active-onset.bson", Dict(:gon => gon, :ρon => ρon, :ron => ron_act, :Ion => Ion_act, :M=>M, :λ=>λ, 
  :L=>L, :τδ=>τδ, :date =>date, :tspan => tspan, :cutoff => cutoff, :ϵ => ϵ) )
end