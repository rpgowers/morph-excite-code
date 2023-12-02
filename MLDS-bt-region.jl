include("activation.jl")

using PyPlot
using BSON: @save

function hopf_finder(v0, ω0, args; ωtol=1e-3, vtol=1e-2)
  hstop = 0
  vhl, Ihl, ωhl = hopf(args; v0=v0[1], ω0=ω0[1], ωtol=ωtol)
  vhu, Ihu, ωhu = hopf(args; v0=v0[2], ω0=ω0[2], ωtol=ωtol)
  if abs(vhu - vhl) < vtol
    Ihl = NaN
    Ihu = NaN
    hstop = 1
  end
  return [Ihl, Ihu], [ωhl, ωhu], [vhl, vhu], hstop
end

function hopfstab_finder(vh, ωh, Ih, args)
  Sh = [0.0, 0.0]
  if Ih[1] > -Inf
    Sh[1] = hopf_stability(vh[1], ωh[1], args)
  end
  if Ih[2] > -Inf
    Sh[2] = hopf_stability(vh[2], ωh[2], args)
  end
  return Sh
end

gin = 2.0:0.05:9.0
gσ = 2.0
N = 25

τδ = [2.5, 5.0, 10.0, 20.0]

Isn = zeros(length(τδ), N, 2)
gin = zeros(length(τδ), N)
Ih = zeros(length(τδ), N, 2)
Sh = zeros(length(τδ), N, 2)
ωh = zeros(length(τδ), N, 2)
Ibt = zeros(length(τδ))
gbt = zeros(length(τδ))

vc, Ic, ρc = [cusp(MLDS_Param())[i][1] for i = 1:3]
gc = gσ * (1 + ρc)

hstop = 0
for j in eachindex(τδ)
  global hstop = 0
  vbt, Ibt[j], ρbt = [bt(MLDS_Param(τδ=τδ[j]))[i][1] for i = 1:3]
  gbt[j] = gσ * (1 + ρbt)
  gin[j, :] .= LinRange(gbt[j], gc, 25)
  ρ = gin[j, :] ./ gσ .- 1

  for i in eachindex(gin[j, :])
    args = MLDS_Param(gL=gσ, ρ=ρ[i], τδ=τδ[j])
    vout, Iout = sn(args)
    if length(Iout) > 1
      Isn[j, i, :] .= Iout[1:2]
    else
      Isn[j, i, :] .= NaN
    end
    if hstop == 0
      Ih[j, i, :], ωh[j, i, :], vh, hstop = hopf_finder([vbt, 10.0], [0.05, 0.05], args; ωtol=1e-3, vtol=1e-2)
      Sh[j, i, :] .= hopfstab_finder(vh, ωh[j, i, :], Ih[j, i, :], args)
    else
      Ih[j, i, :] .= NaN
    end
  end

end

# plot(Isn[:, 1], gin, color="C0", label="SN")
# plot(Isn[:, 2], gin, color="C0")
# plot(Ic, gc, "^", color="k", label="Cusp")

# for j in eachindex(τδ)
#   if j == 1
#     plot(Ih[j, :, 1][findall(x -> x > 0, Sh[j, :, 1])], gin[findall(x -> x > 0, Sh[j, :, 1])], ":", color="C$(j)", alpha=0.5, label="Sub-Hopf")
#     plot(Ih[j, :, 1][findall(x -> x < 0, Sh[j, :, 1])], gin[findall(x -> x < 0, Sh[j, :, 1])], "--", color="C$(j)", alpha=0.5, label="Sup-Hopf")
#   else
#     plot(Ih[j, :, 1][findall(x -> x > 0, Sh[j, :, 1])], gin[findall(x -> x > 0, Sh[j, :, 1])], ":", color="C$(j)", alpha=0.5)
#     plot(Ih[j, :, 1][findall(x -> x < 0, Sh[j, :, 1])], gin[findall(x -> x < 0, Sh[j, :, 1])], "--", color="C$(j)", alpha=0.5)
#   end

#   plot(Ih[j, :, 2][findall(x -> x > 0, Sh[j, :, 2])], gin[findall(x -> x > 0, Sh[j, :, 2])], ":", alpha=0.5, color="C$(j)")
#   plot(Ih[j, :, 2][findall(x -> x < 0, Sh[j, :, 2])], gin[findall(x -> x < 0, Sh[j, :, 2])], "--", alpha=0.5, color="C$(j)")
#   plot(Ibt[j], gbt[j], "o", color="C$(j)", alpha=0.5, label="BT, \$\\tau_\\delta\$ = $(τδ[j]) ms")

# end
# legend(frameon=false)
# show()

@save "data/MLDS-bif-bt-region.bson" gin τδ Isn Ih Sh ωh gbt Ibt gc ρc Ic
