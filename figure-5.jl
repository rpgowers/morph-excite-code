include("fig-master.jl")

DS_data = BSON.load("data/MLDS-prcs.bson")
gon_DS = DS_data[:gon]
Δθ_DS = DS_data[:Δθsim]
θin_DS = DS_data[:θin]
τδ = DS_data[:τδ]
println(τδ)

# C1_data = BSON.load("$(dir2)/prcs/MLC1_prcs_compare.bson")
# gon_C1 = C1_data[:gsel]
# Δθ_C1 = C1_data[:Δθsim]
# θin_C1 = C1_data[:θin]

gsel = [3.0, 4.7, 5.4, 5.9]
idx_DS = findfirst.(isequal.(gsel), (gon_DS,))
# idx_C1 = findfirst.(isequal.(gsel), (gon_C1,))

fig = figure(figsize=(9,9))
fig.subplots_adjust(hspace = 0.3, wspace=0.2)
fig.text(0.075, 0.9, "A", fontsize=font_title)
fig.text(0.5, 0.9, "B", fontsize=font_title)
fig.text(0.075, 0.45, "C", fontsize=font_title)
fig.text(0.5, 0.45, "D", fontsize=font_title)

for i in eachindex(gsel)
	ax = fig.add_subplot(220+i)
	title("\$G_{\\mathrm{in}} = \$ $(gsel[i]) nS", fontsize=font_title)
	axhline(y=0.0, linestyle="--", color="gray", alpha=0.5, linewidth=linewidth)

	# if i == 1 # length(gsel)
	# 	plot(θin_C1, Δθ_C1[idx_C1[i],:]./maximum(abs.(Δθ_C1[idx_C1[i],:])), "-", label="0", color=colors[1], alpha=0.75, linewidth=linewidth)
	# else
	#   plot(θin_C1, Δθ_C1[idx_C1[i],:]./maximum(abs.(Δθ_C1[idx_C1[i],:])), "-", color=colors[1], alpha=0.75, linewidth=linewidth)
	# end

	for j in eachindex(τδ)
		if i == 1 # length(gsel)
		  plot(θin_DS, Δθ_DS[j,idx_DS[i],:]./maximum(abs.(Δθ_DS[j,idx_DS[i],:])), "-", label="$(roundint(τδ[j]))", color=colors[j+1], alpha=0.75, linewidth=linewidth)
	  else
	  	plot(θin_DS, Δθ_DS[j,idx_DS[i],:]./maximum(abs.(Δθ_DS[j,idx_DS[i],:])), "-", color=colors[j+1], alpha=0.75, linewidth=linewidth)
	  end
	end
	
	if i === 1 || i == 3
		ylabel("PRC", fontsize=font_title)
		# ax.set_yticks([0.0])
	else
		# ax.set_yticks([0.0])
	end
	if i == 3 || i == 4
	  xlabel("Phase", fontsize=font_title)
	end
	if i==1
		text(0.4, 0.5, "\$\\tau_\\delta\$ (ms)", fontsize=font_axis)
		legend(frameon=false, ncol = 1, fontsize = font_axis, loc="lower center") # bbox_to_anchor = [0.65, 0.85]
	end
	xticks(fontsize=font_axis)
	yticks(fontsize=font_axis)
	ax.spines["right"].set_visible(false)
  ax.spines["top"].set_visible(false)
end

savefig("figures/figure-5.svg", bbox_inches="tight")
show()