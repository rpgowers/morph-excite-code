include("fig-master.jl")
dir2 = "../spatial-neuron-bifurcations/data"


# C1_data = BSON.load("$(dir2)/bifs/MLC1_bif_general.bson")
S_data = BSON.load("data/MLS_bif_general.bson")

gin_S = S_data[:gin]
Ibt_S = S_data[:Ibt]
gbt_S = S_data[:gbt]
Ih_S = S_data[:Ih]
Sh_S = S_data[:Sh]

DS_data = BSON.load("data/MLDS_bif_general.bson")
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

plot(Isn[:,1], gin_DS, color="k", alpha=alpha, label="SN")
plot(Isn[:,2], gin_DS, color="k", alpha=alpha)
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

ax.spines["right"].set_visible(false)
ax.spines["top"].set_visible(false)

savefig("figures/figure-3.svg", bbox_inches="tight")
show()