include("fig-master.jl")

S_data = BSON.load("data/MLS-bif-bt-region.bson")
Isn = S_data[:Isn]
gin_S = S_data[:gin]
Ibt_S = S_data[:Ibt]
gbt_S = S_data[:gbt]
Ih_S = S_data[:Ih]
Sh_S = S_data[:Sh]

DS_data = BSON.load("data/MLDS-bif-bt-region.bson")
Isn_DS = DS_data[:Isn]
gin_DS = DS_data[:gin]
Ic = DS_data[:Ic]
gc = DS_data[:gc]
Ibt_DS = DS_data[:Ibt]
gbt_DS = DS_data[:gbt]
Ih_DS = DS_data[:Ih]
Sh_DS = DS_data[:Sh]
τδ = DS_data[:τδ]

# Ih_S_shift = Ih_S[60:70, 1] .- Isn[60:70, 1]
# println(Isn[50:80, 1])
# println(Ih_S[50:80, 1])


# ax = subplot(111)
# plot(Isn[:, 1], gin_DS, color="k", alpha=alpha, label="SN (upper)")
plot(0.0, gbt_S, "o", color=colors[1], alpha=alpha, label="BT, S (\$\\tau_\\delta\$ = 0)")
plot(Ih_S[:, 1] .- Isn[:, 1], gin_S)
axvline(x=0.0)
# plot(Ih_S_shift[findall(x -> x > 0, Sh_S[:, 1])], gin_S[findall(x -> x > 0, Sh_S[:, 1])], ":", color=colors[1], alpha=alpha)
# plot(Ih_S_shift[findall(x -> x < 0, Sh_S[:, 1])], gin_S[findall(x -> x < 0, Sh_S[:, 1])], "-", color=colors[1], alpha=alpha)

for j in eachindex(τδ)
  plot(0.0, gbt_DS[j], "o", color=colors[j+1], alpha=0.5, label="BT, \$\\tau_\\delta\$ = $(roundint(τδ[j])) ms")
  plot(Ih_DS[j, :, 1] .- Isn_DS[j, :, 1], gin_DS[j, :])
  # if j == length(τδ)
  #   plot(Ih_DS[j, :, 1][findall(x -> x > 0, Sh_DS[j, :, 1])], gin_DS[findall(x -> x > 0, Sh_DS[j, :, 1])], ":", color=colors[j+1], alpha=alpha, label="Subcr-Hopf")
  #   plot(Ih_DS[j, :, 1][findall(x -> x < 0, Sh_DS[j, :, 1])], gin_DS[findall(x -> x < 0, Sh_DS[j, :, 1])], "-", color=colors[j+1], alpha=alpha, label="Supercr-Hopf")
  # else
  #   plot(Ih_DS[j, :, 1][findall(x -> x > 0, Sh_DS[j, :, 1])], gin_DS[findall(x -> x > 0, Sh_DS[j, :, 1])], ":", color=colors[j+1], alpha=alpha)
  #   plot(Ih_DS[j, :, 1][findall(x -> x < 0, Sh_DS[j, :, 1])], gin_DS[findall(x -> x < 0, Sh_DS[j, :, 1])], "-", color=colors[j+1], alpha=alpha)
  # end

end
# # ax.set_ylim([4.7, 5.8])
# # ax.set_xlim([140.0, 180.0])
# xlabel("\$I_{\\mathrm{ext}}\$ (pA)", fontsize=font_axis)
# ylabel("\$G_{\\mathrm{in}}\$ (nS)", fontsize=font_axis)

show()