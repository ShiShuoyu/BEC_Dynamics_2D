# Physical Constant #
# length: μm; time: ms; mass: Kg #
hbar = 1.05457e-25 # reduced Planc constant # J*s = Kg*μm^2*ms^-1
m0 = 1.6605e-27 # 1a.u. # Kg
ab = 5.2917721e-5 # Bohr radius # μm
mu_B = 9.2732e-12 # Bohr magneton # A*μm^2
grav = 9.8 # gravity # μm*ms^-2
m = 87*m0 # mass of Rb87 # Kg

# characteristic quantities #
x0 = 1 # length # μm
t0 = m*x0^2 / hbar # time # ms
e0 = hbar/t0 # energy # pJ = Kg*μm^2*ms^-2

# dimensionless quantities #
# x -> x/x0, t -> t/t0, E -> E/e0 #
# example: # 
Grav = m*grav/e0 # dimensionless gravity