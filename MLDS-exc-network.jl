using OrdinaryDiffEq, PyPlot, Parameters
using BSON: @save

include("sim-functions.jl")

N = 5 # number of neurons
M = 50 # number of compartments per neuron

args_s = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 89.0, ρ = 0.8, M=M, λ=100.0)
args_h = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 131.79697, ρ = 1.4, M=M, λ=100.0)

v0 = -20.0 .*ones(M+1)
n0 = 0.0
x0 = vcat(n0, v0)
cutoff = 1000.0
T = 10000.0
tspan = (0.0,T+cutoff)

ΔVs = 6.65e-2
ΔVh = 2.665e-4

Np = 25

@time outref,tref,iref = LC_extract(tspan, x0, args_s; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_s1 = outref[end-2]
T0_s = 0.5*(diff(tref)[end-1]+diff(tref)[end])

prob = ODEProblem(network_sim, xp_s1, 2*T0_s, [args_s])
@time sol = solve(prob, Tsit5())

xp_s2 = sol(0.16*T0_s)
xp_s3 = sol(0.32*T0_s)
xp_s4 = sol(0.48*T0_s)
xp_s5 = sol(0.64*T0_s)

x0_s = vcat(xp_s1, xp_s2, xp_s3, xp_s4, xp_s5)
synset = all_to_all(N, ΔVs, xp_s1[2]-1e-1; D=M+2)
cb = CallbackSet(synset...)

probf = ODEProblem(network_sim, x0_s, Np*T0_s, [args_s for i=1:N])
@time solf = solve(probf, Tsit5(), reltol=1e-10, abstol=1e-10, callback=cb, save_idxs=[2, M+4, 2*(M+2)+2, 3*(M+2)+2, 4*(M+2)+2])
t = solf.t

spike_set_s = [find_spikes(solf[i,:], solf.t, xp_s1[2]) for i=1:N]
Δϕs, Ls, Ti = phase_differences(spike_set_s; forward=true)
Ts = Ti[end]

@time outref,tref,iref = LC_extract(tspan, x0, args_h; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_h1 = outref[end-2]
T0_h = 0.5*(diff(tref)[end-1]+diff(tref)[end])

prob = ODEProblem(network_sim, xp_h1, 2*T0_h, [args_h])
@time sol = solve(prob, Tsit5())

xp_h2 = sol(0.16*T0_h)
xp_h3 = sol(0.32*T0_h)
xp_h4 = sol(0.48*T0_h)
xp_h5 = sol(0.64*T0_h)

x0_h = vcat(xp_h1, xp_h2, xp_h3, xp_h4, xp_h5)
synset = all_to_all(N, ΔVh, xp_h1[2]-1e-1; D=M+2)
cb = CallbackSet(synset...)

probf = ODEProblem(network_sim, x0_h, Np*T0_h, [args_h for i=1:N])
@time solf = solve(probf, Tsit5(), reltol=1e-10, abstol=1e-10, callback=cb, save_idxs=[2, M+4, 2*(M+2)+2, 3*(M+2)+2, 4*(M+2)+2])
t = solf.t

spike_set_h = [find_spikes(solf[i,:], solf.t,xp_h1[2]) for i=1:N]
Δϕh, Lh, Ti = phase_differences(spike_set_h; forward=true)
Th = Ti[end]

fig = figure(figsize=(12.0,18.0))
fig.subplots_adjust(hspace = 0.3, wspace=0.3)

ax1 = fig.add_subplot(221)

eventplot(spike_set_s./Ts, linelengths=0.8, colors=["C0", "C1", "C2", "C3", "C4"], lineoffsets=1:1:N)
xlabel("Time (\$T_n\$)", fontsize=16)
ylabel("Neuron Number", fontsize=16)
ax1.minorticks_on()
ax1.set_xticks(0:2:Np+1)
ax1.set_xticks(1:2:Np+1, minor=true)

ax1.tick_params(axis="y", which="minor", left=false)

ax1.set_yticks(1:1:5)
xticks(fontsize=14)
yticks(fontsize=14)
grid(axis="x", which="major")
grid(axis="x", which="minor", linestyle="--")

ax3 = fig.add_subplot(223)

eventplot(spike_set_h./Th, linelengths=0.8, colors=["C0", "C1", "C2", "C3", "C4"], lineoffsets=1:1:N)
xlabel("Time (\$T_n\$)", fontsize=16)
ylabel("Neuron Number", fontsize=16)
ax3.minorticks_on()
ax3.set_xticks(0:2:Np+1)
ax3.set_xticks(1:2:Np+1, minor=true)

ax3.tick_params(axis="y", which="minor", left=false)

ax3.set_yticks(1:1:5)
xticks(fontsize=14)
yticks(fontsize=14)
grid(axis="x", which="major")
grid(axis="x", which="minor", linestyle="--")


ax2 = fig.add_subplot(222)
plot(1:1:Ls, Δϕs[1,:], label="\$ \\theta_1-\\theta_2 \$")
plot(1:1:Ls, Δϕs[2,:], label="\$ \\theta_2-\\theta_3 \$")
plot(1:1:Ls, Δϕs[3,:], label="\$ \\theta_3-\\theta_4 \$")
plot(1:1:Ls, Δϕs[4,:], label="\$ \\theta_4-\\theta_5 \$")
plot(1:1:Ls, Δϕs[5,:], label="\$ \\theta_5-\\theta_1 \$")

axhline(y=1/5, linestyle="--", alpha=0.5, color="k")
xlabel("Spike number", fontsize=16)
ylabel("Phase difference", fontsize=16)
xticks(fontsize=14)
yticks(fontsize=14)
legend(frameon=false, fontsize=12, ncol=2)


ax4 = fig.add_subplot(224)
plot(1:1:Lh, Δϕh[1,:], label="\$ \\theta_1-\\theta_2 \$")
plot(1:1:Lh, Δϕh[2,:], label="\$ \\theta_2-\\theta_3 \$")
plot(1:1:Lh, Δϕh[3,:], label="\$ \\theta_3-\\theta_4 \$")
plot(1:1:Lh, Δϕh[4,:], label="\$ \\theta_4-\\theta_5 \$")
plot(1:1:Lh, Δϕh[5,:], label="\$ \\theta_5-\\theta_1 \$")

axhline(y=1/5, linestyle="--", alpha=0.5, color="k")
xlabel("Spike number", fontsize=16)
ylabel("Phase difference", fontsize=16)
xticks(fontsize=14)
yticks(fontsize=14)
legend(frameon=false, fontsize=12, ncol=2)

@save "data/MLDS-exc-network.bson" spike_set_s spike_set_h Ts Th N Np Δϕs Δϕh Lh Ls

show()