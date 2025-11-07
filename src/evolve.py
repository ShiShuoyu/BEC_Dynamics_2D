import numpy as np
import cupy as cp

from constants import hbar, m, x0, e0, t0, a, ab

def time_step(dt:cp.float32, duration:cp.float32, sampling_interval:cp.int32
              ) -> tuple[cp.float32, cp.float32, cp.int32, cp.int32]:
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
    return (dt / t0, duration / t0, int(cp.round(duration / dt)), int(cp.round(duration / dt / sampling_interval)+1))

def time_evolution(psi:cp.ndarray, U:cp.ndarray, V_sqrt:cp.ndarray, T:cp.ndarray, 
                   dt:cp.float32, g:cp.float32, imaginary_time:np.bool) -> cp.ndarray:
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
        imaginary_time: if True, perform imaginary time evolution
    output:
        psi: evolved wave function, shape (Nx, Ny)
    '''
    loss = 1 - 0.4j if imaginary_time else 1

    if g == 0:
        psi = psi * V_sqrt
        psi = cp.fft.ifft2(cp.fft.fft2(psi) * T)
        psi = psi * V_sqrt
    else:
        V_sqrt_g = cp.exp(-1j * loss * (dt/2) * (U + g*cp.abs(psi)**2))
        psi = psi * V_sqrt_g
        psi = cp.fft.ifft2(cp.fft.fft2(psi) * T)
        V_sqrt_g = cp.exp(-1j * loss * (dt/2) * (U + g*cp.abs(psi)**2))
        psi = psi * V_sqrt_g
    return psi