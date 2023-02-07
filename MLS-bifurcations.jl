include("activation.jl")

using PyPlot, BSON

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
  # Ih[i,:], ωh[i,:], Sh[i,:] = hopf_finder(args)
end

plot(Isn[:,1], gin, color="C0", label="SN")
plot(Isn[:,2], gin, color="C0")
# plot(Ih[:,1][findall(x->x>0, Sh[:,1])], gin[findall(x->x>0, Sh[:,1])], ":", color="C1", label="Sub-Hopf")
# plot(Ih[:,1][findall(x->x<0, Sh[:,1])], gin[findall(x->x<0, Sh[:,1])], "--", color="C1", label="Sup-Hopf")

# plot(Ih[:,2][findall(x->x>0, Sh[:,2])], gin[findall(x->x>0, Sh[:,2])], ":", color="C1")
# plot(Ih[:,2][findall(x->x<0, Sh[:,2])], gin[findall(x->x<0, Sh[:,2])], "--", color="C1")
plot(Ic[1], gc[1], "^", color="C2", label="Cusp")
plot(Ibt, gbt, "o", color="C3", label="BT")

legend(frameon=false)
grid(which="both")
xlabel("\$I_{\\mathrm{ext}}\$ (pA)", fontsize=12)
ylabel("\$g_{\\mathrm{in}}\$ (nS)", fontsize=12)
title("MLS model bifurcation diagram", fontsize=14)

show()