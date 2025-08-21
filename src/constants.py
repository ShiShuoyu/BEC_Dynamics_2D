import numpy as np
import json

with open('constants.json', 'w') as f:
    const = json.load(f)

# physical constants #
hbar = const["hbar"] # reduced Planc constant # J*s = Kg*μm^2*ms^-1
m0 = const["m0"] # 1a.u. # Kg
ab = const["ab"] # Bohr radius # μm
mu_B = const["mu_B"] # Bohr magneton # A*μm^2
grav = const["grav"] # gravity # μm*ms^-2
Ar = const["Ar"] # relative atomic mass # dimensionless
a = const["a"] # scattering length # a_Bohr

m = Ar * m0  # atmoic mass # Kg

# dimensionless units
x0 = 1 # length # μm
t0 = Ar*m0*x0**2 / hbar # time # ms
e0 = hbar/t0 # energy # pJ = Kg*μm^2*ms^-2

# example: # 
# Grav = m*grav/e0 # dimensionless gravity #