include("fig-master.jl")

data_act = BSON.load("data/MLDS-active-prcs.bson")
θin_act = data_act[:θin]
Δθ_act = data_act[:Δθsim]
gon_act = data_act[:gon]

data_pas = BSON.load("data/MLDS-prcs.bson")
θin_pas = data_pas[:θin]
Δθ_pas = data_pas[:Δθsim]

gsel = [3.0, 5.4, 5.5, 5.9]
idx = findfirst.(isequal.(gsel), (gon_act,))

fig = figure(figsize=(9, 9))
fig.subplots_adjust(hspace=0.3, wspace=0.2)
fig.text(0.075, 0.9, "i", fontsize=font_title)
fig.text(0.5, 0.9, "ii", fontsize=font_title)
fig.text(0.075, 0.45, "iii", fontsize=font_title)
fig.text(0.5, 0.45, "iv", fontsize=font_title)

for i in eachindex(gsel)
  ax = fig.add_subplot(220 + i)
  title("\$G_{\\mathrm{in}} = \$ $(gsel[i]) nS", fontsize=font_title)
  axhline(y=0.0, linestyle="--", color="gray", alpha=0.5, linewidth=linewidth)

  plot(θin_pas, Δθ_pas[3, idx[i], :] ./ maximum(abs.(Δθ_pas[3, idx[i], :])), color=colors[1], alpha=0.75, linewidth=linewidth, label="0 (passive)")
  plot(θin_act, Δθ_act[1, idx[i], :] ./ maximum(abs.(Δθ_act[1, idx[i], :])), color=colors[3], alpha=0.75, linewidth=linewidth, label="0.1")
  plot(θin_act, Δθ_act[4, idx[i], :] ./ maximum(abs.(Δθ_act[4, idx[i], :])), color=colors[5], alpha=0.75, linewidth=linewidth, label="1.0")
  if i == 1
    text(0.425, 0.26, "\$\\epsilon\$", fontsize=font_title)
    legend(frameon=false, ncol=1, fontsize=font_axis, loc="lower center") # bbox_to_anchor = [0.65, 0.85]
  end

  if i === 1 || i == 3
    ylabel("PRC", fontsize=font_title)
  else
  end
  if i == 3 || i == 4
    xlabel("Phase", fontsize=font_title)
  end

  ax.spines["right"].set_visible(false)
  ax.spines["top"].set_visible(false)
  xticks(fontsize=font_axis)
  yticks(fontsize=font_axis)
end

savefig("figures/figure-s1.pdf", bbox_inches="tight")
show()