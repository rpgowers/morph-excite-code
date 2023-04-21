using OrdinaryDiffEq, PyPlot, Parameters
using BSON: @save

include("sim-functions.jl")

N = 5 # number of neurons
M = 50 # number of compartments per neuron

args_snic = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 89.0, ρ = 0.8, M=M, λ=100.0)
args_hom = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 131.79697, ρ = 1.4, M=M, λ=100.0)
args_hopf = MLMDS_Param(gL = 2.0, τδ = 2.5, Iext= 169.13389325141907, ρ = 1.88, M=M, λ=100.0)

v0 = -20.0 .*ones(M+1)
n0 = 0.0
x0 = vcat(n0, v0)
cutoff = 2000.0
T = 10000.0
tspan = (0.0,T+cutoff)

# rfind = 1e-3
# Ibound = [165.0, 175.0]
# Iout, rout = find_onset(tspan, x0, args_hopf, cutoff, rfind, Ibound; soma_idx=2, ISI_min = 9, vth=-8.0, Nmax=25, ϵr = 0.01, ϵisi = 10, ϵsp=0.9)

# println(Iout)
# println(rout)

# r, t, vf = rate_measure(tspan, x0, args_hopf, cutoff; soma_idx=2, vth=-8.0, ISImin = 4, ISIfactor=1.1, ϵisi = 10, ϵsp = 0.9, verbose=false, solver=Tsit5(), maxiters=1e6, saveat=[], save_vt = true)

# println(r)
# plot(t, vf)
# show()

ΔVsnic = 6.65e-2
ΔVhom = 2.665e-4
ΔVhopf = 5.639e-5

Np = 25

@time outref,tref,iref = LC_extract(tspan, x0, args_snic; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_snic1 = outref[end-2]
T0_snic = 0.5*(diff(tref)[end-1]+diff(tref)[end])

prob = ODEProblem(network_sim, xp_snic1, 2*T0_snic, [args_snic])
@time sol = solve(prob, Tsit5())

xp_snic2 = sol(0.16*T0_snic)
xp_snic3 = sol(0.32*T0_snic)
xp_snic4 = sol(0.48*T0_snic)
xp_snic5 = sol(0.64*T0_snic)

x0_snic = vcat(xp_snic1, xp_snic2, xp_snic3, xp_snic4, xp_snic5)
synset = all_to_all(N, ΔVsnic, xp_snic1[2]-1e-1; D=M+2)
cb = CallbackSet(synset...)

probf = ODEProblem(network_sim, x0_snic, Np*T0_snic, [args_snic for i=1:N])
@time solf = solve(probf, Tsit5(), reltol=1e-10, abstol=1e-10, callback=cb, save_idxs=[2, M+4, 2*(M+2)+2, 3*(M+2)+2, 4*(M+2)+2])
t = solf.t

spike_set_snic = [find_spikes(solf[i,:], solf.t, xp_snic1[2]) for i=1:N]
Δϕsnic, Lsnic, Ti = phase_differences(spike_set_snic; forward=true)
Tsnic = Ti[end]

@time outref,tref,iref = LC_extract(tspan, x0, args_hom; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_hom1 = outref[end-2]
T0_hom = 0.5*(diff(tref)[end-1]+diff(tref)[end])

prob = ODEProblem(network_sim, xp_hom1, 2*T0_hom, [args_hom])
@time sol = solve(prob, Tsit5())

xp_hom2 = sol(0.16*T0_hom)
xp_hom3 = sol(0.32*T0_hom)
xp_hom4 = sol(0.48*T0_hom)
xp_hom5 = sol(0.64*T0_hom)

x0_hom = vcat(xp_hom1, xp_hom2, xp_hom3, xp_hom4, xp_hom5)
synset = all_to_all(N, ΔVhom, xp_hom1[2]-1e-1; D=M+2)
cb = CallbackSet(synset...)

probf = ODEProblem(network_sim, x0_hom, Np*T0_hom, [args_hom for i=1:N])
@time solf = solve(probf, Tsit5(), reltol=1e-10, abstol=1e-10, callback=cb, save_idxs=[2, M+4, 2*(M+2)+2, 3*(M+2)+2, 4*(M+2)+2])
t = solf.t

spike_set_hom = [find_spikes(solf[i,:], solf.t, xp_hom1[2]) for i=1:N]
Δϕhom, Lhom, Ti = phase_differences(spike_set_hom; forward=true)
Thom = Ti[end]

@time outref,tref,iref = LC_extract(tspan, x0, args_hopf; soma_idx=2, verbose=true, vth = 0.0, mincross=10, solver=Tsit5())
xp_hopf1 = outref[end-2]
T0_hopf = 0.5*(diff(tref)[end-1]+diff(tref)[end])

prob = ODEProblem(network_sim, xp_hopf1, 2*T0_hopf, [args_hopf])
@time sol = solve(prob, Tsit5())

xp_hopf2 = sol(0.16*T0_hopf)
xp_hopf3 = sol(0.32*T0_hopf)
xp_hopf4 = sol(0.48*T0_hopf)
xp_hopf5 = sol(0.64*T0_hopf)

x0_hopf = vcat(xp_hopf1, xp_hopf2, xp_hopf3, xp_hopf4, xp_hopf5)
synset = all_to_all(N, ΔVhopf, xp_hopf1[2]-1e-1; D=M+2)
cb = CallbackSet(synset...)

probf = ODEProblem(network_sim, x0_hopf, Np*T0_hopf, [args_hopf for i=1:N])
@time solf = solve(probf, Tsit5(), reltol=1e-10, abstol=1e-10, callback=cb, save_idxs=[2, M+4, 2*(M+2)+2, 3*(M+2)+2, 4*(M+2)+2])
t = solf.t

spike_set_hopf = [find_spikes(solf[i,:], solf.t,xp_hopf1[2]) for i=1:N]
Δϕhopf, Lhopf, Ti = phase_differences(spike_set_hopf; forward=true)
Thopf = Ti[end]


fig = figure(figsize=(15.0,18.0))
fig.subplots_adjust(hspace = 0.3, wspace=0.3)

ax1 = fig.add_subplot(231)
eventplot(spike_set_snic./Tsnic, linelengths=0.8, colors=["C0", "C1", "C2", "C3", "C4"], lineoffsets=1:1:N)
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

ax3 = fig.add_subplot(232)
eventplot(spike_set_hom./Thom, linelengths=0.8, colors=["C0", "C1", "C2", "C3", "C4"], lineoffsets=1:1:N)
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

ax5 = fig.add_subplot(233)
eventplot(spike_set_hopf./Thopf, linelengths=0.8, colors=["C0", "C1", "C2", "C3", "C4"], lineoffsets=1:1:N)
xlabel("Time (\$T_n\$)", fontsize=16)
ylabel("Neuron Number", fontsize=16)
ax5.minorticks_on()
ax5.set_xticks(0:2:Np+1)
ax5.set_xticks(1:2:Np+1, minor=true)

ax5.tick_params(axis="y", which="minor", left=false)
ax5.set_yticks(1:1:5)
xticks(fontsize=14)
yticks(fontsize=14)
grid(axis="x", which="major")
grid(axis="x", which="minor", linestyle="--")

ax2 = fig.add_subplot(234)
plot(1:1:Lsnic, Δϕsnic[1,:], label="\$ \\theta_1-\\theta_2 \$")
plot(1:1:Lsnic, Δϕsnic[2,:], label="\$ \\theta_2-\\theta_3 \$")
plot(1:1:Lsnic, Δϕsnic[3,:], label="\$ \\theta_3-\\theta_4 \$")
plot(1:1:Lsnic, Δϕsnic[4,:], label="\$ \\theta_4-\\theta_5 \$")
plot(1:1:Lsnic, Δϕsnic[5,:], label="\$ \\theta_5-\\theta_1 \$")

axhline(y=1/5, linestyle="--", alpha=0.5, color="k")
xlabel("Spike number", fontsize=16)
ylabel("Phase difference", fontsize=16)
xticks(fontsize=14)
yticks(fontsize=14)
legend(frameon=false, fontsize=12, ncol=2)

ax4 = fig.add_subplot(235)
plot(1:1:Lhom, Δϕhom[1,:], label="\$ \\theta_1-\\theta_2 \$")
plot(1:1:Lhom, Δϕhom[2,:], label="\$ \\theta_2-\\theta_3 \$")
plot(1:1:Lhom, Δϕhom[3,:], label="\$ \\theta_3-\\theta_4 \$")
plot(1:1:Lhom, Δϕhom[4,:], label="\$ \\theta_4-\\theta_5 \$")
plot(1:1:Lhom, Δϕhom[5,:], label="\$ \\theta_5-\\theta_1 \$")

axhline(y=1/5, linestyle="--", alpha=0.5, color="k")
xlabel("Spike number", fontsize=16)
ylabel("Phase difference", fontsize=16)
xticks(fontsize=14)
yticks(fontsize=14)
legend(frameon=false, fontsize=12, ncol=2)

ax6 = fig.add_subplot(236)
plot(1:1:Lhopf, Δϕhopf[1,:], label="\$ \\theta_1-\\theta_2 \$")
plot(1:1:Lhopf, Δϕhopf[2,:], label="\$ \\theta_2-\\theta_3 \$")
plot(1:1:Lhopf, Δϕhopf[3,:], label="\$ \\theta_3-\\theta_4 \$")
plot(1:1:Lhopf, Δϕhopf[4,:], label="\$ \\theta_4-\\theta_5 \$")
plot(1:1:Lhopf, Δϕhopf[5,:], label="\$ \\theta_5-\\theta_1 \$")

axhline(y=1/5, linestyle="--", alpha=0.5, color="k")
xlabel("Spike number", fontsize=16)
ylabel("Phase difference", fontsize=16)
xticks(fontsize=14)
yticks(fontsize=14)
legend(frameon=false, fontsize=12, ncol=2)

@save "data/MLDS-exc-network.bson" spike_set_snic spike_set_hom spike_set_hopf Tsnic Thom Thopf N Np Δϕsnic Δϕhom Δϕhopf Lsnic Lhom Lhopf

show()