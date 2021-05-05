#!/usr/bin/python3
from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
import control

#System Parameters
N = 20                    # Num Coils
mu_0 = (4*np.pi*10**(-7)) # Absolute permeability constant 
R = 10                    # solenoid resistance [Ohms]
L = 0.05                  # solenoid length of coil [m]
m = 0.05                  # permanent magnet mass [kg]
g = 9.8                   # acceleration due to gravity [m^2]          
A = .01                   # Cross-sectional Area of Electromagnet Plunger [m^2]     
P2 = .005                 # Permanent Magnet Flux [Wb]

# Constant Defs
x_0 = -0.03               # distance of levitation [m]
ground = -0.2             # distance to ground [m] 

# finding v_0      
K_nov =  ((P2*N*A) / (4*np.pi*L*R)) * (1/(x_0**2))
Kv = (P2*N*A) / (4*np.pi*L*R) * (1/(x_0**2))
Kx_nov = ((P2*N*A) / (4*np.pi*L*R)) * (2/(x_0**3))        
v_0 = (m*g) / (x_0*Kx_nov - Kv + K_nov)   

# plug v_0 back in 
K = (P2*N*A) / (4*np.pi*L*R) * (v_0/(x_0**2))
Kx = (P2*N*A) / (4*np.pi*L*R) * (2*v_0/(x_0**3))

#Force Calculations
i = v_0/R                 # Current through electromagnet [A]
P1 = (mu_0*((N*A*i)/L))   # Electromagnet Flux [Wb]

x = x_0                   # Permanent Magnet Location [m]
v = v_0                   # Voltage through electromagnet [V]

F_m = 2*(((P2*N*A*v_0)/(4*np.pi*L*R*(x_0**2)))-((2*v_0*P2*N*A)/((x_0**3)*4*np.pi*L*R))*(x-x_0)+((P2*N*A)/((x_0**2)*4*np.pi*L*R))*(v-v_0))
F_ext = F_m-(m*g)


print(f"Gravity Force: {m*g}")
print(f"Magnetic Force: {F_m}")
print(f"Sum of Forces: {F_ext}")
print(f"v_0: {v_0}")


# Root Locus
L = control.TransferFunction((Kv), (m, 0, Kx))
rlist, klist = control.root_locus(L) 
plt.show()

# Uncontrolled System Simulation
step_vs = [-3, -1, -0.1, 0, 0.1, 1, 3]
plt.figure()
for v in step_vs:
    G = signal.TransferFunction([v*Kv], [m, 0, Kx])
    t, out = signal.step(G)
    out[out > -x_0] = -x_0
    out[out < ground] = ground
    plt.plot(t, out, label=f'v_rel={v}')

plt.axhline(y=ground, color='k', linestyle='--', label='ground')
plt.axhline(y=-x_0, color='k', linestyle='--', label='electromagnet')
plt.xlabel("Time [s]")
plt.ylabel("Distance [m]")
plt.title('Magnet Distance from Stable Position with Variable Starting Positions')
plt.legend()
plt.show()
