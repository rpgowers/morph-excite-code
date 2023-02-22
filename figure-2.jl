include("fig-master.jl")

data = BSON.load("data/DS-impedance.bson")
gin = data[:gin]
f = data[:f]
τδ = data[:τδ]
l = data[:l]
ZS = data[:ZS]
ZI = data[:ZI]
ZF = data[:ZF]

fig = figure(figsize=(1.25*fig_width, panel_height))
fig.subplots_adjust(wspace=0.4)

fig.text(0.07, 0.875, "A", fontsize=font_title)
ax1 = subplot(131)
for i in eachindex(gin)
  plot(f, abs.(ZS[:,i]), "--", color=colors[2*i-1], alpha=alpha, linewidth=linewidth)
  plot(f, abs.(ZI[:,i,2]), color=colors[2*i-1], alpha=alpha, linewidth=linewidth, label="$(Int(gin[i]))")
end
plot([],[], "--", color="k", alpha=0.5, linewidth=linewidth, label="S")
plot([],[], "-", color="k", alpha=0.5, linewidth=linewidth, label="DS")
ax1.spines["right"].set_visible(false)
ax1.spines["top"].set_visible(false)
xlim([0,250.0])
xlabel("Frequency (Hz)", fontsize=font_axis)
ylabel("\$|Z_{\\mathrm{in}}|\$ (GΩ)", fontsize=font_axis)
legend(fontsize=font_legend, frameon=false, ncol=2)
text(70.0, 0.26, "\$G_{\\mathrm{in}}\$ (nS)", fontsize=font_axis)
text(160.0, 0.26, "Model", fontsize=font_axis)

fig.text(0.35, 0.875, "B", fontsize=font_title)
ax2 = subplot(132)
for i in eachindex(τδ)
  plot(f, abs.(ZI[:,3,i]), color=colors[2*i-1], alpha=alpha, linewidth=linewidth, label="$(τδ[i]*1000)")
end
ax2.spines["right"].set_visible(false)
ax2.spines["top"].set_visible(false)
xlim([0,250.0])
xlabel("Frequency (Hz)", fontsize=font_axis)
legend(fontsize=font_legend, frameon=false)
text(170.0, 0.13, "\$\\tau_\\delta\$ (ms)", fontsize=font_axis)

fig.text(0.64, 0.875, "C", fontsize=font_title)
ax3 = subplot(133)
for i=1:length(τδ)
  plot(f, abs.(ZF[:,3,i]), color=colors[2*i-1], alpha=alpha, linewidth=linewidth, label="$(l[i])")
end
ax3.spines["right"].set_visible(false)
ax3.spines["top"].set_visible(false)
xlim([0,250.0])
xlabel("Frequency (Hz)", fontsize=font_axis)
legend(fontsize=font_legend, frameon=false)
text(210.0, 0.13, "\$\\ell\$", fontsize=font_axis)

savefig("figures/figure-2.svg", bbox_inches="tight")
show()