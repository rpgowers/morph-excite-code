import brian2 as b2
from pymodules.ML_params import *
import pymodules.cable_theory as ct
import copy
import numpy as np

def gin(neuron, Iin, v0, tmax):
    neuron.v = v0
    neuron.I[0] = Iin
    mon = b2.StateMonitor(neuron, 'v', record=0)
    t = mon.t
    b2.run(tmax, report='text')
    return Iin/(mon.v[0][-1]-v0)

def Zin(neuron, q, f, v0, tmax): # returns the magnitude and phase of the impedance at a single frequency
    neuron.v = v0
    neuron.qin[0] = q # amplitude of sinusoid
    neuron.freq[0] = f # frequency of sinusoid
    mon = b2.StateMonitor(neuron, 'v', record=0)
    Tp = 1/f
    if Tp > tmax:
        tmax = Tp*2.5 # get at least two periods per simulation

    t = mon.t
    b2.run(tmax)
    dt = t[1]-t[0]
    Np = int(1//(f*dt)) # number of time steps in a period
    vsec = mon.v[0][-Np:]
    vbase = b2.mean(vsec) # need the baseline voltage 
    vmax = max(vsec)
    idxV = b2.where(vsec==vmax)
    Iin = b2.sin(2*b2.pi*f*t[-Np:])*q*neuron.area[0]
    imax = max(Iin)
    idxI = b2.where(Iin==imax)
    mag = abs((vmax-vbase)/imax)
    phase = (idxI[0][0]-idxV[0][0])*2*b2.pi/Np
    if phase > b2.pi: # avoid high positive phases
        phase = phase-2*b2.pi # make the phase negative
        
    return mag, phase

def Z0L(neuron, Iin, v0, tmax): # this gives Z00 and Z0L for a dendritic stem
    neuron.v = v0
    neuron.I[0] = Iin
    mon = b2.StateMonitor(neuron, 'v', record=range(0,len(neuron)))
    t = mon.t
    b2.run(tmax, report='text')
    Zt = (mon.v[:,-1]-v0)/Iin
    Z0 = Zt[0]
    ZL = min(Zt)
    return Z0, ZL

def stem_reduce(name, Nstem, Iin, v0, tmax, form="stem"):
    l = b2.zeros(Nstem)
    d = b2.zeros(Nstem)*b2.um
    for i in range(1,Nstem+1):
        print(i)
        morph = b2.spatialneuron.morphology.Morphology.from_swc_file(f"morphs/{name}_{form}{i}.swc", spherical_soma=True)
        stem = b2.SpatialNeuron(morphology=morph, model=eqs, Cm=Cm, Ri=Ri, method='exponential_euler')
        stem[0::].gL = gL0
        stem.EL = EL0
        Z0, ZL = Z0L(stem, Iin, v0, tmax)
        l[i-1] = ct.l_equiv(Z0, ZL)
        d[i-1] = ct.d_equiv(Z0, Ri, gL0, l[i-1])
        
    return l, d

def soma_equiv(name):
    morph = b2.spatialneuron.morphology.Morphology.from_swc_file(f"morphs/{name}_somas.swc", spherical_soma=True)
    soma = b2.SpatialNeuron(morphology=morph, model=eqs, Cm=Cm, Ri=Ri,method='exponential_euler')
    return sum(soma.area)

def f_measure(t, v, cutoff, vmh=-5.0, verbose=False, ϵisi = 10):
    T = t[-1]-t[cutoff]
    N = len(t)
    vout = v[cutoff:]-vmh
    p = b2.zeros(len(vout)-1)
    for i in range(0,len(vout)-1):
        p[i] = vout[i]*vout[i+1]
        
    idx = b2.where(p<0)[0]
    Δm = b2.diff(idx)
    ISI = [Δm[2*i]+Δm[2*i+1] for i in range(len(Δm)//2)]
    
    if verbose == True:
        print(ISI)
    
    if len(Δm) < 2:
        f = 0.0*b2.Hz
    elif abs(ISI[-1]-ISI[-2]) > ϵisi:
        f = 0.0*b2.Hz # ISI changes too much
    elif N-cutoff-idx[-1] > 2*(Δm[-1]+Δm[-2]):
        f = 0.0*b2.Hz
    else:
        f = sum(p <= 0)/(2*T)
    return f


def findrate(neuron, Ibound, cutoff, tmax, rfind, v0 = -60.0*b2.mV, vmh = -8.0, Nmax = 10, ϵr = 0.05, ϵisi = 10):
    rout = 0
    Iout = 0
    
    mon = b2.StateMonitor(neuron, 'v', record=0)
    neuron.v = v0
    neuron.n = 0.0
    neuron.I[0] = Ibound[0]
    b2.run(tmax)
    t = mon.t
    rl = f_measure(t, mon.v[0]/b2.mV, cutoff, vmh=vmh)
    if rl > rfind:
        print("Invalid lower bound")
    
    mon = b2.StateMonitor(neuron, 'v', record=0)
    neuron.v = v0
    neuron.n = 0.0
    neuron.I[0] = Ibound[1]
    b2.run(tmax)
    t = mon.t
    ru = f_measure(t, mon.v[0]/b2.mV, cutoff, vmh = vmh, ϵisi = ϵisi)
    if ru < rfind:
        print("Invalid upper bound")
        
    Itri = [Ibound[0], 0.5*(Ibound[0]+Ibound[1]), Ibound[1]]
    for i in range(Nmax):
        mon = b2.StateMonitor(neuron, 'v', record=0)
        neuron.v = v0
        neuron.n = 0.0
        neuron.I[0] = Itri[1]
        b2.run(tmax)
        t = mon.t
        rm = f_measure(t, mon.v[0]/b2.mV, cutoff, vmh = vmh, ϵisi = ϵisi)
        if rm > rfind:
            rout = rm
            Iout = Itri[1]
            Itri[2] = Itri[1]
            Itri[1] = 0.5*(Itri[0]+Itri[1])
        else:
            Itri[0] = Itri[1]
            Itri[1] = 0.5*(Itri[1]+Itri[2])
            
        if abs(rout-rfind)/rfind < ϵr:
            break
        
    return Iout, rout

def LC_setup(v, n, vth, mlc): # this finds the references voltages and active variable for a LC
    vσ = v[0]
    dv = b2.diff(vσ)
    dvπ = b2.array([dv[i]*dv[i+1] for i in range(0,len(dv)-1)])
    idx = b2.where(dvπ < 0)[0]
    vref = []
    nref = []
    iref = []
    for i in idx:
        if vσ[i] > vth:
            vref.append(v[:,i])
            nref.append(n[i])
            iref.append(i)
            
    Nlc = b2.diff(iref)
    Np = b2.mean(Nlc[-mlc:])
    vstart = b2.mean(vref[-mlc:],axis=0)
    nstart = b2.mean(nref[-mlc:])
    return vstart, nstart, Np


def PRC_extract(neuron, Iin, in_idx, Δv, P, Np, dt, vstart, nstart, vth):
    Tp = Np*dt # time in a period/ISI
    Npint = int(b2.ceil(Np))
    
    neuron.v = vstart*b2.mV
    neuron.n = nstart
    neuron.I[0] = Iin
    mon = b2.StateMonitor(neuron, ('v','n'), record=range(0,len(neuron)))
    t = mon.t
    b2.run(P*Tp)
    
    vbase = mon.v
    nbase = mon.n
    
    Δθ = b2.zeros(len(in_idx))
    for k in range(0,len(in_idx)):
        vnew = copy.copy(vbase[:,in_idx[k]])
        vnew[0]+=Δv
        neuron.v = vnew
        neuron.n = nbase[:, in_idx[k]]
        neuron.I[0] = Iin
        mon = b2.StateMonitor(neuron, ('v','n'), record=0)
        t = mon.t
        b2.run(P*Tp-in_idx[k]*dt)
        vshift = mon.v[0]/b2.mV
        
        idx_new = b2.argmin(vshift[-Npint:])
        idx_old = b2.argmin(vbase[0,-Npint:])
        Δθ[k] = (idx_old-idx_new)/Np
        if abs(Δθ[k]) > 0.5: # this keeps the PRC values in between -0.5 and 0.5
            Δθ[k] = Δθ[k]-b2.sign(Δθ[k])
            
    return Δθ