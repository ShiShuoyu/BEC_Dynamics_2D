import numpy as np
import cupy as cp

from constants import hbar, m, x0, e0, t0, a, ab

def time_step(dt:np.float32, duration:np.float32, sampling_interval:np.int32
              ) -> tuple[np.float32, np.float32, np.int32, np.int32]:
    '''
    functionality:
        build up the time structure for the simulation
    input:
        dt: time step in ms
        duration: total duration of the simulation in ms
        sampling_interval: after how many time steps to sample the wave function
    output:
        dt: dimensionless time step
        duration: dimensionless total duration
        n_steps: number of time steps in the simulation
        n_samples: number of samples
    '''
    return (dt / t0, duration / t0, int(np.round(duration / dt)), int(np.round(duration / dt / sampling_interval)))

def time_evolution(psi:cp.ndarray, U:cp.ndarray, V_sqrt:cp.ndarray, T:cp.ndarray, 
                   dt:np.float32, Num:np.int32, omega_z:np.float32) -> np.ndarray:
    '''
    functionality:
        evolve the wave function psi under the potential V, kinetic operator T, and interaction erengy g|ψ|^2
    input:
        psi: wave function, shape (Nx, Ny)
        U: potential energy, shape (Ny, Nx) # dimensionless
        V_sqrt: potential energy operator over dt/2, shape (Nx, Ny)
        T: kinetic energy operator over dt, shape (Nx, Ny)
        dt: time step # ms
        Num: number of atoms
        omega_z: trapping frequency along z direction
    output:
        psi: evolved wave function, shape (Nx, Ny)
    '''
    # interacting strength #
    g = (4*np.pi*Num * hbar**2 * a*ab / m / e0 / x0**3
         ) * np.sqrt(m*omega_z/2/np.pi/hbar * x0
                     ) # correction for squeezing a 3D wavepacket (along z direction) into a x0 μm layer

    if g == 0:
        psi = psi * V_sqrt
        psi = cp.ifft2(cp.fft2(psi) * T)
        psi = psi * V_sqrt
    else:
        V_sqrt_g = cp.exp(-1j * (dt/2) * (U + g*cp.abs(psi)**2))
        psi = psi * V_sqrt_g
        psi = cp.ifft2(cp.fft2(psi) * T)
        psi = psi * V_sqrt_g
    return psi