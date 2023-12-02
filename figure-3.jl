include("fig-master.jl")

S_data = BSON.load("data/MLS-bif-general.bson")
gin_S = S_data[:gin]
Ibt_S = S_data[:Ibt]
gbt_S = S_data[:gbt]
Ih_S = S_data[:Ih]
Sh_S = S_data[:Sh]

DS_data = BSON.load("data/MLDS-bif-general.bson")
Isn = DS_data[:Isn]
gin_DS = DS_data[:gin]
Ic = DS_data[:Ic]
gc = DS_data[:gc]
Ibt_DS = DS_data[:Ibt]
gbt_DS = DS_data[:gbt]
Ih_DS = DS_data[:Ih]
Sh_DS = DS_data[:Sh]
τδ = DS_data[:τδ]

fig = figure(figsize=(9, 12))
fig.subplots_adjust(hspace=0.4, wspace=0.4)

ax = subplot2grid(shape=(3, 3), loc=(0, 0), colspan=3, rowspan=2)
fig.text(0.075, 0.9, "A", fontsize=font_title)

plot(Isn[:, 2], gin_DS, "--", color="k", alpha=alpha, label="SN (lower)")
plot(Isn[:, 1], gin_DS, color="k", alpha=alpha, label="SN (upper)")
plot(Ic, gc, "^", color="k", alpha=alpha, label="Cusp")

plot(Ibt_S, gbt_S, "o", color=colors[1], alpha=alpha, label="BT, S (\$\\tau_\\delta\$ = 0)")
plot(Ih_S[:, 1][findall(x -> x > 0, Sh_S[:, 1])], gin_S[findall(x -> x > 0, Sh_S[:, 1])], ":", color=colors[1], alpha=alpha)
plot(Ih_S[:, 1][findall(x -> x < 0, Sh_S[:, 1])], gin_S[findall(x -> x < 0, Sh_S[:, 1])], "-", color=colors[1], alpha=alpha)

plot(Ih_S[:, 2][findall(x -> x > 0, Sh_S[:, 2])], gin_S[findall(x -> x > 0, Sh_S[:, 2])], ":", color=colors[1], alpha=alpha)
plot(Ih_S[:, 2][findall(x -> x < 0, Sh_S[:, 2])], gin_S[findall(x -> x < 0, Sh_S[:, 2])], "-", color=colors[1], alpha=alpha)

idx_S = length(filter(!isnan, Ih_S[:, 2]))
plot([Ih_S[idx_S, 1], Ih_S[idx_S, 2]], [gin_S[idx_S], gin_S[idx_S]], "-", color=colors[1], alpha=alpha)

for j in eachindex(τδ)
  plot(Ibt_DS[j], gbt_DS[j], "o", color=colors[j+1], alpha=0.5, label="BT, \$\\tau_\\delta\$ = $(roundint(τδ[j])) ms")
  if j == length(τδ)
    plot(Ih_DS[j, :, 1][findall(x -> x > 0, Sh_DS[j, :, 1])], gin_DS[findall(x -> x > 0, Sh_DS[j, :, 1])], ":", color=colors[j+1], alpha=alpha, label="Subcr-Hopf")
    plot(Ih_DS[j, :, 1][findall(x -> x < 0, Sh_DS[j, :, 1])], gin_DS[findall(x -> x < 0, Sh_DS[j, :, 1])], "-", color=colors[j+1], alpha=alpha, label="Supercr-Hopf")
  else
    plot(Ih_DS[j, :, 1][findall(x -> x > 0, Sh_DS[j, :, 1])], gin_DS[findall(x -> x > 0, Sh_DS[j, :, 1])], ":", color=colors[j+1], alpha=alpha)
    plot(Ih_DS[j, :, 1][findall(x -> x < 0, Sh_DS[j, :, 1])], gin_DS[findall(x -> x < 0, Sh_DS[j, :, 1])], "-", color=colors[j+1], alpha=alpha)
  end

  plot(Ih_DS[j, :, 2][findall(x -> x > 0, Sh_DS[j, :, 2])], gin_DS[findall(x -> x > 0, Sh_DS[j, :, 2])], ":", alpha=alpha, color=colors[j+1])
  plot(Ih_DS[j, :, 2][findall(x -> x < 0, Sh_DS[j, :, 2])], gin_DS[findall(x -> x < 0, Sh_DS[j, :, 2])], "-", alpha=alpha, color=colors[j+1])

  idx = length(filter(!isnan, Ih_DS[j, :, 2]))
  plot([Ih_DS[j, idx, 1], Ih_DS[j, idx, 2]], [gin_DS[idx], gin_DS[idx]], "-", color=colors[j+1], alpha=alpha)

end
legend(frameon=false, fontsize=font_legend)
xlabel("\$I_{\\mathrm{ext}}\$ (pA)", fontsize=font_axis)
ylabel("\$G_{\\mathrm{in}}\$ (nS)", fontsize=font_axis)
title("Bifurcation diagram", fontsize=font_title)

ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)

data_in = BSON.load("data/MLDS-bif-bt.bson")
τδ_in = data_in[:τδ]
gbt_in = data_in[:gbt]
Ibt_in = data_in[:Ibt]
τbtc = data_in[:τbtc]
_, idx_btc = findmin(abs.(τδ_in .- τbtc[1]))
τmax = 20.0
idx_max = findfirst.(isequal.(τmax), (τδ_in,))

axin = subplot2grid(shape=(3, 3), loc=(2, 0))
fig.text(0.075, 0.33, "B", fontsize=font_title)

axin.plot(τbtc[1], gc, "^", color="k", alpha=alpha, label="Cusp")
axin.plot(τδ_in[1:idx_btc-1], gbt_in[1:idx_btc-1], linewidth=linewidth, color="k", label="BT (upper)")
axin.plot(τδ_in[idx_btc:idx_max[1]], gbt_in[idx_btc:idx_max[1]], "--", linewidth=linewidth, color="k", label="BT (lower)")

axin.plot(0.0, gbt_S, "o", color=colors[1], alpha=alpha)
for j in eachindex(τδ)
  axin.plot(τδ[j], gbt_DS[j], "o", color=colors[j+1], alpha=0.5)
end
axin.set_xlabel("\$\\tau_{\\delta}\$ (ms)", fontsize=font_legend)
axin.set_ylabel("\$G_{\\mathrm{in}}^{\\mathrm{BT}}\$ (nS)", fontsize=font_legend)

axin.spines["right"].set_visible(false)
axin.spines["top"].set_visible(false)

Sbt_data = BSON.load("data/MLS-bif-bt-region.bson")
Isn = Sbt_data[:Isn]
gin_S = Sbt_data[:gin]
Ibt_S = Sbt_data[:Ibt]
gbt_S = Sbt_data[:gbt]
Ih_S = Sbt_data[:Ih]
Sh_S = Sbt_data[:Sh]

DS_data = BSON.load("data/MLDS-bif-bt-region.bson")
Isn_DS = DS_data[:Isn]
gin_DS = DS_data[:gin]
Ibt_DS = DS_data[:Ibt]
gbt_DS = DS_data[:gbt]
Ih_DS = DS_data[:Ih]
Sh_DS = DS_data[:Sh]

axbt = subplot2grid(shape=(3, 3), loc=(2, 1))
fig.text(0.375, 0.33, "C", fontsize=font_title)

axbt.set_xlabel("\$I_{\\mathrm{ext}}-I_{\\mathrm{SN (upper)}}\$ (pA)", fontsize=font_legend)
axbt.set_ylabel("\$G_{\\mathrm{in}}\$ (nS)", fontsize=font_legend)
axbt.spines["right"].set_visible(false)
axbt.spines["top"].set_visible(false)

axbt.plot(0.0, gbt_S, "o", color=colors[1], alpha=alpha)
# axbt.plot(Ih_S[:, 1] .- Isn[:, 1], gin_S, color=colors[1])
plot((Ih_S[:, 1].-Isn[:, 1])[findall(x -> x > 0, Sh_S[:, 1])], gin_S[findall(x -> x > 0, Sh_S[:, 1])], ":", color=colors[1], alpha=alpha)
plot((Ih_S[:, 1].-Isn[:, 1])[findall(x -> x < 0, Sh_S[:, 1])], gin_S[findall(x -> x < 0, Sh_S[:, 1])], "-", color=colors[1], alpha=alpha)

for j = 1:length(τδ)-1
  axbt.plot(0.0, gbt_DS[j], "o", color=colors[j+1], alpha=alpha)
  axbt.plot((Ih_DS[j, :, 1].-Isn_DS[j, :, 1])[findall(x -> x > 0, Sh_DS[j, :, 1])], gin_DS[j, :][findall(x -> x > 0, Sh_DS[j, :, 1])], ":", color=colors[j+1], alpha=alpha)
  axbt.plot((Ih_DS[j, :, 1].-Isn_DS[j, :, 1])[findall(x -> x < 0, Sh_DS[j, :, 1])], gin_DS[j, :][findall(x -> x < 0, Sh_DS[j, :, 1])], "-", color=colors[j+1], alpha=alpha)
  # axbt.plot(Ih_DS[j, :, 1] .- Isn_DS[j, :, 1], gin_DS[j, :], color=colors[j+1], alpha=alpha)
end
axbt.axvline(x=0.0, color="k")

axbt2 = subplot2grid(shape=(3, 3), loc=(2, 2))
fig.text(0.65, 0.33, "D", fontsize=font_title)

axbt2.set_xlabel("\$I_{\\mathrm{ext}}-I_{\\mathrm{SN (lower)}}\$ (pA)", fontsize=font_legend)
axbt2.set_ylabel("\$G_{\\mathrm{in}}\$ (nS)", fontsize=font_legend)
axbt2.spines["right"].set_visible(false)
axbt2.spines["top"].set_visible(false)

axbt2.plot(0.0, gbt_DS[end], "o", color=colors[end], alpha=alpha)
axbt2.plot((Ih_DS[end, :, 1].-Isn_DS[end, :, 2])[findall(x -> x > 0, Sh_DS[end, :, 1])], gin_DS[end, :][findall(x -> x > 0, Sh_DS[end, :, 1])], ":", color=colors[end], alpha=alpha)
axbt2.plot((Ih_DS[end, :, 1].-Isn_DS[end, :, 2])[findall(x -> x < 0, Sh_DS[end, :, 1])], gin_DS[end, :][findall(x -> x < 0, Sh_DS[end, :, 1])], "-", color=colors[end], alpha=alpha)
axbt2.axvline(x=0.0, color="k")

savefig("figures/figure-3.svg", bbox_inches="tight")
show()