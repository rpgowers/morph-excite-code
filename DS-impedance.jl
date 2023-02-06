using BSON

function S_impedance(G_σ, C_σ, ω)
  return 1/(G_σ+1im*ω*C_σ)
end

function DS_inf_impedance_soma(G_σ, C_σ, ρ, τ_δ, ω) 
  # input impedance at soma due to somatic current injection, semi-infinite dendrite
  return 1/(G_σ+1im*ω*C_σ+ρ*G_σ*sqrt(1+1im*ω*τ_δ))
end

function DS_fin_impedance_soma(G_σ, C_σ, ρ, l, τ_δ, ω)
  # input impedance at soma due to somatic current injection, finite dendrite
  γ = sqrt(1+1im*ω*τ_δ)
  return 1/(G_σ+1im*ω*C_σ+ρ*G_σ*γ*tanh(l*γ))
end

gσ = 2.0
Cσ = 20.0/1000.0 # to give the time constant in seconds
τδ = [2.5, 10.0, 40.0]./1000 # in seconds

gL = [4.0,6.0,8.0] # mS/cm^2
l = [0.2, 1.0, 5.0]

ρ = (gL.-gσ)./gσ

p = 0:0.1:5
ω = 10.0 .^p # rad/s
f = ω/(2*pi)

ZS = zeros(Complex, length(ω),length(gL))
ZI = zeros(Complex, length(ω),length(ρ), length(τδ))
ZF = zeros(Complex, length(ω), length(ρ), length(l))
for i in eachindex(ρ)
  ZS[:,i] = S_impedance.(gL[i], Cσ, ω)
  for j in eachindex(τδ)
  	ZI[:,i,j] = DS_inf_impedance_soma.(gσ, Cσ, ρ[i], τδ[j], ω)
  end
  for j in eachindex(l)
  	ρF = (gL[i]-gσ)/(gσ*tanh(l[j]))
  	ZF[:,i,j] = DS_fin_impedance_soma.(gσ, Cσ, ρF, l[j], τδ[2], ω)
  end
end

bson("data/DS-impedance.bson", Dict(:gin => gL, :f => f, :τδ => τδ, :l => l, :ZS => ZS, :ZI => ZI, :ZF => ZF))