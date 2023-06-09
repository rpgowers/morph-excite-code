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

ax = subplot(111)

plot(Isn[:,2], gin_DS, "--", color="k", alpha=alpha, label="SN (lower)")
plot(Isn[:,1], gin_DS, color="k", alpha=alpha, label="SN (upper)")
plot(Ic, gc, "^", color="k", alpha=alpha, label="Cusp")

plot(Ibt_S, gbt_S, "o", color=colors[1], alpha=alpha, label="BT, S (\$\\tau_\\delta\$ = 0)")
plot(Ih_S[:,1][findall(x->x>0, Sh_S[:,1])], gin_S[findall(x->x>0, Sh_S[:,1])], ":", color=colors[1], alpha=alpha)
plot(Ih_S[:,1][findall(x->x<0, Sh_S[:,1])], gin_S[findall(x->x<0, Sh_S[:,1])], "-", color=colors[1], alpha=alpha)

plot(Ih_S[:,2][findall(x->x>0, Sh_S[:,2])], gin_S[findall(x->x>0, Sh_S[:,2])], ":", color=colors[1], alpha=alpha)
plot(Ih_S[:,2][findall(x->x<0, Sh_S[:,2])], gin_S[findall(x->x<0, Sh_S[:,2])], "-", color=colors[1], alpha=alpha)

idx_S = length( filter(!isnan, Ih_S[:,2]) )
plot([Ih_S[idx_S, 1], Ih_S[idx_S, 2]], [gin_S[idx_S],  gin_S[idx_S]], "-", color=colors[1], alpha=alpha)

for j in eachindex(τδ)
	plot(Ibt_DS[j], gbt_DS[j], "o", color=colors[j+1], alpha=0.5, label="BT, \$\\tau_\\delta\$ = $(roundint(τδ[j])) ms")
	if j==length(τδ)
		plot(Ih_DS[j,:,1][findall(x->x>0, Sh_DS[j,:,1])], gin_DS[findall(x->x>0, Sh_DS[j,:,1])], ":", color=colors[j+1], alpha=alpha, label="Subcr-Hopf")
		plot(Ih_DS[j,:,1][findall(x->x<0, Sh_DS[j,:,1])], gin_DS[findall(x->x<0, Sh_DS[j,:,1])], "-", color=colors[j+1], alpha=alpha, label="Supercr-Hopf")
	else
		plot(Ih_DS[j,:,1][findall(x->x>0, Sh_DS[j,:,1])], gin_DS[findall(x->x>0, Sh_DS[j,:,1])], ":", color=colors[j+1], alpha=alpha)
		plot(Ih_DS[j,:,1][findall(x->x<0, Sh_DS[j,:,1])], gin_DS[findall(x->x<0, Sh_DS[j,:,1])], "-", color=colors[j+1], alpha=alpha)
	end

	plot(Ih_DS[j,:,2][findall(x->x>0, Sh_DS[j,:,2])], gin_DS[findall(x->x>0, Sh_DS[j,:,2])], ":", alpha=alpha, color=colors[j+1])
	plot(Ih_DS[j,:,2][findall(x->x<0, Sh_DS[j,:,2])], gin_DS[findall(x->x<0, Sh_DS[j,:,2])], "-", alpha=alpha, color=colors[j+1])

  idx = length( filter(!isnan, Ih_DS[j,:,2]) )
  plot([Ih_DS[j,idx,1], Ih_DS[j,idx,2]], [gin_DS[idx],  gin_DS[idx]], "-", color=colors[j+1], alpha=alpha)
		
end
legend(frameon=false, fontsize=font_legend)
xlabel("\$I_{\\mathrm{ext}}\$ (pA)", fontsize=font_axis)
ylabel("\$G_{\\mathrm{in}}\$ (nS)", fontsize=font_axis)
title("Bifurcation diagram", fontsize=font_title)

data_in = BSON.load("data/MLDS-bif-bt.bson")
τδ_in = data_in[:τδ]
gbt_in = data_in[:gbt]
Ibt_in = data_in[:Ibt] 
τbtc = data_in[:τbtc]
_, idx_btc = findmin(abs.(τδ_in.-τbtc[1]))
τmax = 20.0
idx_max = findfirst.(isequal.(τmax), (τδ_in,))

ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)

axin = ax.inset_axes([0.73, 0.16, 0.27, 0.27])

axin.plot(τbtc[1],	gc, "^", color="k", alpha=alpha, label="Cusp")
axin.plot(τδ_in[1:idx_btc-1], gbt_in[1:idx_btc-1], linewidth = linewidth, color="k", label="BT (upper)")
axin.plot(τδ_in[idx_btc:idx_max[1]], gbt_in[idx_btc:idx_max[1]], "--", linewidth = linewidth, color="k", label="BT (lower)")

axin.plot(0.0, gbt_S, "o", color=colors[1], alpha=alpha)
for j in eachindex(τδ)
	axin.plot(τδ[j],	gbt_DS[j], "o", color=colors[j+1], alpha=0.5)
end
axin.set_xlabel("\$\\tau_{\\delta}\$ (ms)", fontsize=font_legend)
axin.set_ylabel("\$G_{\\mathrm{in}}^{\\mathrm{BT}}\$ (nS)", fontsize=font_legend)

axin.spines["right"].set_visible(false)
axin.spines["top"].set_visible(false)

savefig("figures/figure-3.svg", bbox_inches="tight")
show()