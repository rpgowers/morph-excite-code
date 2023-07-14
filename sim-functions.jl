include("activation.jl")

using Statistics, DiffEqCallbacks, Setfield, Parameters

function neuron_sim!(du, u, p, t)
	du .= Fall_dyn(u,p)
  nothing
end

function network_sim(du, u, p, t)
  @unpack M, dims = p[1]
  N = length(p)
  D = M+dims # number of dimensions in each cell
  for j=1:N
    du[D*(j-1)+1:D*j] .= Fall_dyn(u[D*(j-1)+1:D*j],p[j])
  end
  nothing
end

function rate_measure(tspan, x0, args, cutoff; soma_idx=2, vth=-5.0, ISImin = 4, ISIfactor=1.1, ϵisi = 10, ϵsp = 1.0, verbose=false, solver=Tsit5(), maxiters=1e6, saveat=[], save_vt = true)
  T = tspan[2]
  probf = ODEProblem(neuron_sim!, x0, tspan, args)
  solf = solve(probf, solver, reltol=1e-8, abstol=1e-8, maxiters=maxiters, saveat=saveat, save_idxs=soma_idx, dense=false)
  # solf = solve(probf, reltol=1e-8, abstol=1e-8, maxiters=maxiters, saveat=saveat, save_idxs=soma_idx)
  vf = solf[1,:]
  t = solf.t
  cutidx = findfirst(t .> cutoff)

  Vred = vf[cutidx+1:end]
  tred = t[cutidx+1:end]

  idx_max = tp_locations(Vred, 1)
  v_max = Vred[idx_max]
  idx_sp = idx_max[v_max .> vth]
  v_sp = v_max[v_max .> vth]

  idx_min = tp_locations(Vred, -1)
  v_min = Vred[idx_min]

  tp_bound = minimum([length(v_max), length(v_min)])

  t_sp = tred[idx_sp]
  ISI = diff(t_sp)

  if length(ISI) < ISImin
    r = 0.0
    verbose == true ? println("Too few ISIs") : 0
  else
    v_osc = v_max[end-ISImin:end].-v_min[end-ISImin:end]
    Q = sum(diff(v_osc)/v_osc[1] .< 0.0)/length(diff(v_osc)) # proportion of oscillations that are decreasing
    verbose == true ? println(Q) : 0
      
    if T-t_sp[end] > ISIfactor*ISI[end]
      r = 0.0
      verbose == true ? println("Firing eventually ceases") : 0
    elseif maximum(abs.(ISI[end].-ISI[ISImin-2:end-1])) > ϵisi
      r = 0
      verbose == true ? println("ISI length not consistent") : 0
    elseif Q > ϵsp # minimum(v_osc[end]./v_osc) < 1-ϵsp
      r = 0
      verbose == true ? println("Spike height decreases") : 0
    else
      r = (length(ISI)+1)/(T-cutoff)
    end
  end

  if save_vt == false
    vf = nothing
    t = nothing
  end
  return r, t, vf
end

function find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, ISI_min = 4, vth=-5.0, Nmax=50, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=1e6, saveat=[])
  rout = 0.0
  Iout = 0.0

  args_l = @set args.Iext = Ibound[1]
  r_l, _, _ = rate_measure(tspan, x0, args_l, cutoff; soma_idx=soma_idx, vth=vth, ISImin = ISI_min, ISIfactor=1.1, ϵisi = ϵisi, ϵsp = ϵsp, verbose=false, solver=Tsit5(), maxiters=maxiters, saveat=saveat, save_vt = false)
  if r_l > rfind
    println("Invalid lower value")
    return Iout, rout
  end
  args_u = @set args.Iext = Ibound[2]
  r_u, _, _ = rate_measure(tspan, x0, args_u, cutoff; soma_idx=soma_idx, vth=vth, ISImin = ISI_min, ISIfactor=1.1, ϵisi = ϵisi, ϵsp = ϵsp, verbose=false, solver=Tsit5(), maxiters=maxiters, saveat=saveat, save_vt = false)
  
  if r_u < rfind
    println("Initial invalid upper value")
    return Iout, rout
  end

  Itri = [Ibound[1], 0.5*(Ibound[1]+Ibound[2]), Ibound[2]]
  for i=1:Nmax
    args_m = @set args.Iext = Itri[2]
    r_m, _, _ = rate_measure(tspan, x0, args_m, cutoff; soma_idx=soma_idx, vth=vth, ISImin = ISI_min, ISIfactor=1.1, ϵisi = ϵisi, ϵsp = ϵsp, verbose=false, solver=Tsit5(), maxiters=maxiters, saveat=saveat, save_vt = false)
    if r_m > rfind # firing rate too high
      rout = r_m
      Iout = Itri[2]
      Itri[3] = Itri[2] # midpoint becomes the new upper bound
      Itri[2] = 0.5*(Itri[1]+Itri[2]) # average of lower and midpoint
    else # firing rate too low
      Itri[1] = Itri[2] # midpoint becomes the new lower bound
      Itri[2] = 0.5*(Itri[2]+Itri[3])
    end
    if abs(rout-rfind)/rfind < ϵr
      break
    end
  end
  
  return Iout, rout
end

function bsnl_extract(μ, μup, μlow, niter, args_in, x0; maxiters = 1e6, rfind = 1e-3, cutoff = 2000.0, tmax = 10000.0)
  μon = 0.0
  ron = 0.0
  Ion = 0.0
  args = deepcopy(args_in)
  tspan = (0.0, tmax+cutoff)
  for i=1:niter
    μmid = 0.5*(μlow+μup)
    setfield!(args, μ, μmid)
    local vsn, Isn = sn(args)
    local Isn_max = maximum(Isn)
    Ibound = [Isn_max-0.5, Isn_max+0.5]
    Iout, rout = find_onset(tspan, x0, args, cutoff, rfind, Ibound; soma_idx=2, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=1.0, solver=Tsit5(), maxiters=maxiters, saveat=[])

    if Iout < Isn_max # midpoint becomes new upper bound
      μon = μmid
      ron = rout
      Ion = Iout
      μup = μmid
      μmid = 0.5*(μlow+μmid)
    else # midpoint becomes new lower bound
      μlow = μmid
      μmid = 0.5*(μup+μmid)
    end
  end
  return μon, Ion, ron
end

function tp_locations(V, sign; shift=0)
  # sign = +1 for a maximum and -1 for a minimum
  dV = diff(V) # finds consecutive differences between voltage time-course
  d2V = diff(dV)
  dVπ = [dV[i]*dV[i+1] for i=1:length(dV)-1] # this finds the product of consective dV/dt values
  idx = findall(dVπ .< 0)
  idx_tp = []
  for i in idx
    if sign*d2V[i] < 0
      push!(idx_tp, i+shift)
    end
  end
  return idx_tp
end

function LC_extract(tspan, x0, args; soma_idx=2, verbose=true, vth = 0.0, mincross=5, solver=Tsit5(), dtmax=1e-1, maxiters=1e6, saveat = [])
  probf = ODEProblem(neuron_sim!, x0, tspan, args)
  solf = solve(probf, solver, reltol=1e-8, abstol=1e-8, dtmax=dtmax, maxiters=maxiters, saveat=saveat)
  vf = solf[soma_idx,:]
  t = solf.t

  idx_tp = tp_locations(vf, 1, shift=1)

  outref = []
  iref = []
  tref = []

  for i in idx_tp
    if vf[i] > vth
      push!(outref, solf[:,i])
      push!(tref, solf.t[i])
      push!(iref,i)
    end
  end

  if length(iref) < mincross # not firing maxima, some subthreshold oscillations
    if verbose == true
      println("Insufficient firing maxima, only $(length(iref)) found")
    end
    return 0, 0, 0
  end

  return outref, tref, iref
end

function PRC_sim(T0, xp, args, ΔV, θin; soma_idx=2, Ps = 10, solver=Tsit5(), dtmax=0.1, maxiters=1e6, tp_type = -1, saveat=[])
  tin = T0.*θin
  Δθ = zeros(length(θin))
  probf = ODEProblem(neuron_sim!, xp, Ps*T0, args) # this is the origial trace
  solf = solve(probf, solver, reltol=1e-8, abstol=1e-8, maxiters=maxiters, dtmax=dtmax, saveat=saveat)
  to = solf.t
  Vo = solf[soma_idx,:]
  cuto = findfirst(to .> T0)
  idx_mino = tp_locations(Vo[cuto:end], tp_type)
  t_mino = to[cuto.+idx_mino]

  for i in eachindex(θin)
    affect!(integrator) = integrator.u[soma_idx] += ΔV
    cb2 = PresetTimeCallback(tin[i],affect!)
    sols = solve(probf, solver, reltol=1e-8, abstol=1e-8, callback=cb2, maxiters=maxiters, dtmax=dtmax, saveat=saveat) # this is the shifted trace
    Vs = sols[soma_idx,:]
    ts = sols.t
    cuts = findfirst(ts .> T0)
    idx_mins = tp_locations(Vs[cuts:end], tp_type)
    t_mins = ts[cuts.+idx_mins]

    if length(idx_mins) < Ps-2
      println("Not enough minima for θin = $(θin[i])")
    end
        

    ΔT = t_mino[Ps-2]-t_mins[Ps-2]
    ΔT < -0.5*T0 ? Δθ[i] = 1+ΔT/T0 : Δθ[i] = ΔT/T0

  end
  return Δθ
end

# Create functions such that, when a spike occurs in neuron i, the potential is increased by ΔV in neuron j
function spike(u,t,integrator,idx,vth)
  integrator.u[idx] > vth && integrator.uprev[idx] < vth
end

function synapse!(integrator, idx, ΔV)
  integrator.u[idx] += ΔV
end

function all_to_all(N, ΔV, Vsp;D=2)
  synapses = []
  spikes = []
  connections = []
  for i=1:N
    syn!(integrator) = synapse!(integrator, D*(i-1)+2, ΔV)
    sp(u,t,integrator) = spike(u,t,integrator,D*(i-1)+2,Vsp)
    push!(synapses, syn!)
    push!(spikes, sp)
  end
  for i=1:N
    for j=1:N
      if i != j
        connect = DiscreteCallback(spikes[i], synapses[j])
        push!(connections, connect)
      end
    end
  end
  return connections
end

function find_spikes(u,t,xp)
  t_sp = []
  for i=2:length(u)
    if u[i] > xp-1e-1 && u[i-1] < xp-1e-1
      push!(t_sp, t[i])
    end
  end
  return Float64.(t_sp)
end

function array_forward(spike_set)
  N = length(spike_set)
  L = minimum([length(spike_set[i]) for i=1:N])
  spike_array = zeros(N,L)
  for i=1:N
    for j=1:L
      spike_array[i,j] = spike_set[i][j]
    end
  end
  return spike_array
end

function array_backward(spike_set)
  N = length(spike_set)
  L = minimum([length(spike_set[i]) for i=1:N])
  spike_array = zeros(N,L)
  for i=1:N
    for j=1:1:L
      spike_array[i,L-j+1] = spike_set[i][end-j+1]
    end
  end
  return spike_array
end

function phase_differences(spike_set; forward=true)
  if forward == true
    spike_array = array_forward(spike_set)
  else
    spike_array = array_backward(spike_set)
  end
  N, L = size(spike_array)
  Δϕ = zeros(N,L)
  Ti = diff(mean(spike_array, dims=1), dims=2)
  for j=1:L-1
    E = sort(spike_array[:,j])
    Δϕ[:,j] .= vcat(diff(E)./Ti[j], 1-sum(diff(E)./Ti[j]))
  end
  E = sort(spike_array[:,L])
  Δϕ[:,L] .= vcat(diff(E)./Ti[end], 1-sum(diff(E)./Ti[end]))
  return Δϕ, L, Ti
end

function sync_measure(Δϕ, L)
  N = length(Δϕ[:,1])
  ϵ2 = zeros(L)
  for j=1:L
    ϵ2[j] = sum((Δϕ[i,j]-1/N)^2 for i=1:N)
  end
  R = sqrt.(N.*ϵ2./(N-1))
  return R, ϵ2
end