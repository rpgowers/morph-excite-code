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

gin = 2.0:0.05:9.4
Isn = zeros(length(gin), 2)
Ih = zeros(length(gin), 2)
Sh = zeros(length(gin), 2)
ωh = zeros(length(gin), 2)

vc, Ic, gc = cusp(MLS_Param())
vbt, Ibt, gbt = bt(MLS_Param())

vbtc, Ibtc, gbtc, Cbtc = btc(MLS_Param())
println(Cbtc)

for i in eachindex(gin)
  args = MLS_Param(gL=gin[i])
  vout, Iout = sn(args)
  if length(Iout) > 1
    Isn[i,:] .= Iout[1:2]
  else
    Isn[i,:] .= NaN
  end
  Ih[i,:], ωh[i,:], Sh[i,:] = hopf_finder(args)
end

plot(Isn[:,1], gin, color="C0", label="SN")
plot(Isn[:,2], gin, color="C0")
plot(Ih[:,1][findall(x->x>0, Sh[:,1])], gin[findall(x->x>0, Sh[:,1])], ":", color="C1", label="Sub-Hopf")
plot(Ih[:,1][findall(x->x<0, Sh[:,1])], gin[findall(x->x<0, Sh[:,1])], "--", color="C1", label="Sup-Hopf")

plot(Ih[:,2][findall(x->x>0, Sh[:,2])], gin[findall(x->x>0, Sh[:,2])], ":", color="C1")
plot(Ih[:,2][findall(x->x<0, Sh[:,2])], gin[findall(x->x<0, Sh[:,2])], "--", color="C1")
plot(Ic[1], gc[1], "^", color="C2", label="Cusp")
plot(Ibt, gbt, "o", color="C3", label="BT")

legend(frameon=false)
grid(which="both")
xlabel("\$I_{\\mathrm{ext}}\$ (pA)", fontsize=12)
ylabel("\$g_{\\mathrm{in}}\$ (nS)", fontsize=12)
title("MLS model bifurcation diagram", fontsize=14)

show()

@save "data/MLS-bif-general.bson" gin Isn Ih Sh ωh gbt Ibt gc Ic Cbtc
