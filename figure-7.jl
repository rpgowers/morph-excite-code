using Parameters
using BSON: @load

include("sim-functions.jl")
include("fig-master.jl")
snic_col = colors[5]
hom_col = colors[3]

@load "data/MLDS-exc-network.bson" spike_set_snic spike_set_hom Tsnic Thom N Np Lhom Δϕhom Lsnic Δϕsnic Lhopf Δϕhopf
@load "data/MLDS-exc-prcs.bson" θin Δθsnic Δθhom Δθhopf

Rhom, ϵ2 = sync_measure(Δϕhom, Lhom)
Rsnic, ϵ2 = sync_measure(Δϕsnic, Lsnic)
Rhopf, ϵ2 = sync_measure(Δϕhopf, Lhopf)

fig, ax = subplots(2, 2, figsize=(9,6), gridspec_kw=Dict("height_ratios" => [1, 2]) )
fig.subplots_adjust(hspace = 0.4, wspace=0.4)

fig.text(0.075, 0.9, "A", fontsize=font_title)
fig.text(0.5, 0.9, "B", fontsize=font_title)
fig.text(0.075, 0.575, "C", fontsize=font_title)
fig.text(0.5, 0.575, "D", fontsize=font_title)

ax1 = ax[1,1] # subplot(221)
Tsnic_max = ceil(maximum(spike_set_snic./Tsnic)[end])
ax1.eventplot(spike_set_snic./Tsnic, linelengths=0.8, linewidths=2.5, colors=[snic_col for i=1:N], lineoffsets=1:1:N)
ax1.set_title("SNIC Network", fontsize=font_title)
ax1.spines["right"].set_visible(false)
ax1.spines["top"].set_visible(false)
ax1.set_xlabel("Time (# ISI)", fontsize=font_axis)
ax1.set_ylabel("Neuron", fontsize=font_axis)
ax1.set_yticks([])
ax1.set_xlim([Tsnic_max-6, Tsnic_max])
ax1.set_xticks(Tsnic_max-6:1:Tsnic_max)
ax1.set_xticklabels(0:1:6)
ax1.tick_params(labelsize=font_axis)

ax2 = ax[1,2] # subplot(222)
Thom_max = ceil(maximum(spike_set_hom./Thom)[end])
ax2.eventplot(spike_set_hom./Thom, linelengths=0.8, linewidths=2.5, colors=[hom_col for i=1:N], lineoffsets=1:1:N)
ax2.set_title("HOM Network", fontsize=font_title)
ax2.spines["right"].set_visible(false)
ax2.spines["top"].set_visible(false)
ax2.set_xlabel("Time (# ISI)", fontsize=font_axis)
ax2.set_ylabel("Neuron", fontsize=font_axis)
ax2.set_yticks([])
ax2.set_xlim([Thom_max-6, Thom_max])
ax2.set_xticks(Thom_max-6:1:Thom_max)
ax2.set_xticklabels(0:1:6)
ax2.tick_params(labelsize=font_axis)

ax3 = ax[2,1] # subplot(223)
ax3.plot(1:1:Lsnic, Rsnic, linewidth=linewidth, color=snic_col, label="SNIC")
ax3.plot(1:1:Lhom, Rhom, linewidth=linewidth, color=hom_col, label="HOM")
# ax3.plot(1:1:Lhopf, Rhopf, linewidth=linewidth, color="black", label="Hopf")

ax3.spines["right"].set_visible(false)
ax3.spines["top"].set_visible(false)
ax3.set_xlabel("Time (# ISI)", fontsize=font_axis)
ax3.set_ylabel("Network synchrony measure, \$ R \$", fontsize=font_axis)
ax3.legend(frameon=false, fontsize=font_axis)
ax3.tick_params(labelsize=font_axis)

ax4 = ax[2,2] #subplot(224)
ax4.plot(θin, Δθsnic, linewidth=linewidth, color=snic_col, label="SNIC")
ax4.plot(θin, Δθhom, linewidth=linewidth, color=hom_col, label="HOM")
ax4.spines["right"].set_visible(false)
ax4.spines["top"].set_visible(false)
ax4.legend(frameon=false, fontsize=font_axis)
ax4.set_xlabel("Phase", fontsize=font_axis)
ax4.set_ylabel("PRC", fontsize=font_axis)
ax4.tick_params(labelsize=font_axis)

savefig("figures/figure-7.svg", bbox_inches="tight")
show()