#!/usr/bin/python3
from matplotlib import pyplot as plt
import numpy as np
from scipy import signal
import control

#System Parameters
mu_0 = (4*np.pi*10**(-7)) # Absolute permeability constant 
m = 0.003                  # Permanent magnet mass [kg]
g = 9.8                   # Acceleration due to gravity [m^2]          
R = 2.41
C = .5*(8.70819*(10**-6))  # lumped parameter of P2NA/4piLR
# Constant Defs
x_0 = -0.03               # Distance of levitation [m]
ground = -0.2             # Distance to ground [m] 

# finding v_0      
K_nov =  C * (1/(x_0**2))
Kv = C * (1/(x_0**2))
Kx_nov = C * (2/(x_0**3))        
v_0 = (m*g) / (x_0*Kx_nov - Kv + K_nov)   

# plug v_0 back in 
K = C * (v_0/(x_0**2))
Kx = C * (2*v_0/(x_0**3))

#Force Calculations
i = v_0/R                 # Current through electromagnet [A]
# P1 = (mu_0*((N*A*i)/L))   # Electromagnet Flux [Wb]

x = x_0                   # Permanent Magnet Location [m]
v = v_0                   # Voltage through electromagnet [V]

#F_m = 2*(((P2*N*A*v_0)/(4*np.pi*L*R*(x_0**2)))-((2*v_0*P2*N*A)/((x_0**3)*4*np.pi*L*R))*(x-x_0)+((P2*N*A)/((x_0**2)*4*np.pi*L*R))*(v-v_0))
F_m = 2*((C*v_0)/(x_0**2))-((2*v_0*C/((x_0**3))*(x-x_0)))+((C)/((x_0**2))*(v-v_0))
F_ext = F_m-(m*g)


print(f"Gravity Force: {m*g}")
print(f"Magnetic Force: {F_m}")
print(f"Sum of Forces: {F_ext}")
print(f"v_0: {v_0}")


# Root Locus
L = control.TransferFunction((Kv), (m, 0, Kx))
# print(L)
# rlist, klist = control.root_locus(L) 
# plt.show()

# Uncontrolled System Simulation
step_vs = [-3, -1, -0.3, 0, 0.3, 1, 3]
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
# plt.show()


pole = np.sqrt(abs(Kx/m))    # OL poles @ +/- this value
print(pole)


# PID Controller
Kp = 1
Ki = 19
Kd = .01
Ck = 1
Vin = 3
v_bar = Vin-v_0

C = signal.TransferFunction([Ck*Kd, Ck*Kp, Ck*Ki], [1, 0])

Gcl = signal.TransferFunction([Ck*v_bar*Kv*Kd, Ck*v_bar*Kv*Kp, Ck*v_bar*Kv*Ki], 
                              [m, (v_bar*Kv*Kd), (Kx+v_bar*Kv*Kp), v_bar*Kv*Ki])

# L = control.TransferFunction([Ck*v_bar*Kv*Kd, Ck*v_bar*Kv*Kp, Ck*v_bar*Kv*Ki], 
#                               [m, (v_bar*Kv*Kd), (Kx+v_bar*Kv*Kp), v_bar*Kv*Ki])
# print(L)
# rlist, klist = control.root_locus(L)
# plt.show()

z1 = (-Kp + np.sqrt(Kp**2 - 4*Kd*Ki))/(2*Kd)
z2 = (-Kp - np.sqrt(Kp**2 - 4*Kd*Ki))/(2*Kd)

print(z1, z2)

#PD Controller
Kp = 2.5
Kd = .1
Ck = 1
Vin = 3
v_bar = Vin-v_0
Vin = v_0-X

C = signal.TransferFunction([Ck*Kd, Ck*Kp], [1])

Gcl = signal.TransferFunction(  [Ck*Kv*v_bar*Kd, Ck*Kv*v_bar*Kp],
                                [m, v_bar*Kv*Kd, (Kx + v_bar*Kv*Kp)])

L = control.TransferFunction(   [Ck*Kv*v_bar*Kd, Ck*Kv*v_bar*Kp],
                                [m, v_bar*Kv*Kd, (Kx + v_bar*Kv*Kp)])
print(L)
rlist, klist = control.root_locus(L)
plt.show()