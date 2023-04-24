include("activation.jl")
using BSON: @save

function hopf_finder(v0, ω0, args; ωtol=1e-3, vtol=1e-2)
  hstop = 0
  vhl, Ihl, ωhl = hopf(args; v0=v0[1], ω0=ω0[1], ωtol = ωtol)
  vhu, Ihu, ωhu = hopf(args; v0=v0[2], ω0=ω0[2], ωtol = ωtol)
  if abs(vhu-vhl) < vtol
    Ihl = NaN
    Ihu = NaN
    hstop = 1
  end
  return [Ihl, Ihu], [ωhl, ωhu], [vhl, vhu], hstop
end

function hopfstab_finder(vh,ωh,Ih,args)
  Sh = [0.0, 0.0]
  if Ih[1] > -Inf
    Sh[1] = hopf_stability(vh[1],ωh[1],args)
  end
  if Ih[2] > -Inf
    Sh[2] = hopf_stability(vh[2],ωh[2],args)
  end
  return Sh
end

τδ = 0.0:0.1:40.0
gσ = 2.0
  
vbtc, Ibtc, ρbtc, τbtc = btc(MLDS_Param())
vc, Ic, ρc = [cusp(MLDS_Param())[i][1] for i=1:3]
gc = gσ*(1+ρc)

gbt = zeros(length(τδ))
Ibt = zeros(length(τδ))
ghmax = zeros(length(τδ))
Ihmax = zeros(length(τδ))
gbaut = zeros(length(τδ))
Ibaut = zeros(length(τδ))

for i in eachindex(τδ)
  vbt, Ibt[i], ρbt = [bt(MLDS_Param(τδ=τδ[i]))[k][1] for k=1:3]
  gbt[i] = gσ*(1+ρbt)
  gset = gbt[i]:0.01:9.5
  global hstop = 0
  Ilast = [0.0, 0.0]
  Slast = 0.0
  if τδ[i] == 0

    for j in eachindex(gset)
      args_S = MLS_Param(gL=gset[j])
      vout, Iout, ωout = hopf(args_S)
      
      if length(Iout) == 2
        Sout = hopf_stability(vout[1], abs(ωout[1]), args_S)
        if Sout*Slast < 0.0
          gbaut[i] = 0.5*(gset[j-1]+gset[j])
          Ibaut[i] = 0.5*(Ilast[1]+Iout[1])
        end
        Ilast .= Iout
        Slast = Sout
      else
        ghmax[i] = gset[j-1]
        Ihmax[i] = 0.5*(Ilast[1]+Ilast[2])
        if gbaut[i] == 0.0
          gbaut[i] = NaN
          Ibaut[i] = NaN
        end
        break
      end

    end

  else
    
    
    for j in eachindex(gset)
      ρset = gset[j]/gσ -1
      args = MLDS_Param(gL=gσ, ρ=ρset, τδ=τδ[i])
      Iout, ωout, vout, hstop = hopf_finder([vbt, 10.0], [0.05, 0.05], args; ωtol=1e-3, vtol=1e-2)
      Sout = hopfstab_finder(vout, ωout, Iout, args)
      if Sout[1]*Slast < 0.0
        gbaut[i] = 0.5*(gset[j-1]+gset[j])
        Ibaut[i] = 0.5*(Ilast[1]+Iout[1])
      end
      if hstop == 1
        ghmax[i] = gset[j-1]
        Ihmax[i] = 0.5*(Ilast[1]+Ilast[2])
        if gbaut[i] == 0.0
          gbaut[i] = NaN
          Ibaut[i] = NaN
        end
        break
      else
        Ilast .= Iout
        Slast = Sout[1]
      end

    end

  end

end


@save "data/MLDS-bif-bt.bson" τδ gbt Ibt gc ρc Ic τbtc ghmax Ihmax gbaut Ibaut