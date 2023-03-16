import numpy as np

def l_equiv(Z0, ZL):
    return 0.5*np.arccosh(2*(Z0/ZL)**2-1)

def d_equiv(Z0, Ri, gL, l):
    return (2/(np.pi*Z0)*np.sqrt(Ri/gL)/np.tanh(l))**(2/3)

def DS_fin_impedance(Gσ, Ginf, τσ, τδ, l, ω):
    Yδ = sum( Ginf[i]*np.sqrt(1+ω*τδ*1j)*np.tanh(l[i]*np.sqrt(1+ω*τδ*1j)) for i in range(0,len(Ginf)) )
    return 1/(Gσ*(1+ω*τσ*1j)+Yδ)