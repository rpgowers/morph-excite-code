import brian2 as b2

gL0 = 5e-5*b2.siemens/b2.cm**2
Cm = 0.5*b2.uF/b2.cm**2
Ri = 150*b2.ohm*b2.cm
gf0 = 2*gL0
gs0 = 4*gL0
EL0 = -60*b2.mV
Es = -80*b2.mV
Dn = 8.7*b2.mV
An = 12.0*b2.mV
Dm = 9.0*b2.mV
Am = -1.2*b2.mV
Ef = 120.0*b2.mV
phi = 1/(15*b2.ms)

eqs = '''
Im = gL * (EL-v) + gf*minf*(Ef-v)+gs*n*(Es-v) + Is: amp/meter**2
dn/dt = (ninf-n)*psi : 1
psi = phi*cosh( (v-An)/(4*Dn) ) : 1/second
minf = 1/(1+exp(-(v-Am)/Dm)) : 1
ninf = 1/(1+exp(-(v-An)/Dn)) : 1
Is = qin*sin(2*pi*freq*t) : amp/meter**2
qin : amp/meter**2
freq : 1/second
gf : siemens/meter**2
gs : siemens/meter**2
gL : siemens/meter**2
EL : volt
I : amp (point current)
'''