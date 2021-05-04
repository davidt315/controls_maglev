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
i = 5/R                   # current through electromagnet [A]
A = .01                   # Cross-sectional Area of Electromagnet Plunger [m^2]
P1 = mu_0*(N*A*i)/L       # Electromagnet Flux [Wb]
P2 = .005                 # Permanent Magnet Flux [Wb]
Xout = 0                  # Output distance between magnet and x_0 [m]
Vin = 0                   # Difference in Voltage from v_0 [V]
Vh = 0                    # Measured Voltage from Hall-Effect Sensor [V]


# Constant Defs
x_0 = 0.03                # distance of levitation [m]

# finding v_0      
K_nov =  ((P2*N*A) / (4*np.pi*L*R)) * (1/(x_0**2))
Kv = (P2*N*A) / (4*np.pi*L*R) * (1/(x_0**2))
Kx_nov = ((P2*N*A) / (4*np.pi*L*R)) * (-2/(x_0**3))
v_0 = (-m*g) / (x_0*Kx_nov + Kv - K_nov)          # voltage at distance of levitation [V]

# plug v_0 back in 
K = (P2*N*A) / (4*np.pi*L*R) * (v_0/(x_0**2))
Kx = (P2*N*A) / (4*np.pi*L*R) * (-2*v_0/(x_0**3))

#Force Calculations
i = v_0/R
P1 = mu_0*(N*A*i)/L  
F_m = (1/(mu_0*4*np.pi))*((P1*P2)/(x_0**2))
F_ext = F_m-(m*g)
print(f"Gravity Force: {m*g}")
print(f"Magnetic Force: {F_m}")
print(f"External Force: {F_ext}")
print(f"v_0: {v_0}")


# whatever
L = control.TransferFunction((Kv), (1, 0, -Kx))
print(L)
rlist, klist = control.root_locus(L) 
plt.show()

# Uncontrolled TF
G = signal.TransferFunction([Kv], [m, 0, -Kx])
print(G)
t, x = signal.step(G)

plt.figure()
plt.plot(t, x, "b", label='Transfer Function Solution')
plt.xlabel("Time [s]")
plt.ylabel("Distance [m]")
plt.title('title')
plt.legend()
plt.show()
