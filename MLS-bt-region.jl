include("activation.jl")

using PyPlot
using BSON: @save

function hopf_finder(args)
  Ih = zeros(2)
  ωh = zeros(2)
  Sh = zeros(2)
  vout, Iout, ωout = hopf(args)
  if length(Iout) == 2
    Ih .= Iout[1:2]
    ωh .= ωout[1:2]
  else
    Ih .= NaN
    ωh .= NaN
  end
  if Ih[1] > -Inf
    Sh[1] = hopf_stability(vout[1], abs(ωout[1]), args)
  end
  if Ih[2] > -Inf
    Sh[2] = hopf_stability(vout[2], abs(ωout[2]), args)
  end

  return Ih, ωh, Sh
end

vc, Ic, gc = cusp(MLS_Param())
vbt, Ibt, gbt = bt(MLS_Param())

println(gbt)
println(gc)

gin = LinRange(gbt[1], gc[1], 25)
Isn = zeros(length(gin), 2)
Ih = zeros(length(gin), 2)
Sh = zeros(length(gin), 2)
ωh = zeros(length(gin), 2)

for i in eachindex(gin)
  args = MLS_Param(gL=gin[i])
  vout, Iout = sn(args)
  if length(Iout) > 1
    Isn[i, :] .= Iout[1:2]
  else
    Isn[i, :] .= NaN
  end
  Ih[i, :], ωh[i, :], Sh[i, :] = hopf_finder(args)
end

@save "data/MLS-bif-bt-region.bson" gin Isn Ih Sh ωh gbt Ibt gc Ic


