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