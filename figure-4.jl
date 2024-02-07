include("fig-master.jl")
using BSON: @load

@load "data/MLDS-bif-bt.bson" τδ gbt Ibt gc ρc Ic τbtc gbaut ghmax

data_S_snl = BSON.load("data/MLS-onset-bsnl.bson")
gsnl_S = data_S_snl[:gon]

data_DS_snl = BSON.load("data/MLDS-onset-bsnl.bson")
τsnl_DS = data_DS_snl[:τsel]
gsnl_DS = data_DS_snl[:gon]

τsnl = vcat(0.0, τsnl_DS)
gsnl = vcat(gsnl_S, gsnl_DS)
τsel = [0, 2.5, 5, 10, 20.0]
gsel = [4.7, 5.4, 5.9]

τmax = 20.0
idx_max = findfirst.(isequal.(τmax), (τδ,))
idx_snl = findfirst.(isequal.(τsnl), (τδ,))
_, idx_btc = findmin(abs.(τδ.-τbtc[1]))
ax = subplot(111)

for i in eachindex(τsel)
  plot([τsel[i] for j=1:length(gsel)], gsel, "o", color=colors[i], alpha=alpha)
end

plot(τδ[1:idx_btc-1], gbt[1:idx_btc-1], linewidth = linewidth, color="k", label="BT (higher)")
plot(τδ[idx_btc:idx_max[1]], gbt[idx_btc:idx_max[1]], "--", linewidth = linewidth, color="k", label="BT (lower)")

plot(τsnl, gsnl, linewidth=linewidth, color="darkgrey", alpha=0.75, label="SNL")
fill_between(τsnl, gsnl, gbt[idx_snl], color="darkgrey", alpha=0.75)
fill_between(τδ[1:idx_max[1]], gbaut[1:idx_max[1]], gbt[1:idx_max[1]], color="lightgrey", alpha=0.5)
axhline(y=gc, color="C6", linestyle="-", linewidth=linewidth, alpha=0.75, label="Cusp")
xlabel("\$\\tau_\\delta\$ (ms)", fontsize=font_axis)
ylabel("\$G_{\\mathrm{in}}\$ (nS)", fontsize=font_axis)
xticks(fontsize=font_axis)
yticks(fontsize=font_axis)
legend(frameon=false, fontsize=font_axis)
text(0.1, 5.2, "Subcr-Hopf", fontsize=font_axis)
text(3.0, 4.9, "HOM", fontsize=font_axis)
text(7.0, 4.4, "SNIC", fontsize=font_axis)
text(15.0, 5.7, "Supercr-Hopf", fontsize=font_axis)

ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)
ylim([4.1, 6.0])

savefig("figures/figure-4.svg", bbox_inches="tight")
show()
