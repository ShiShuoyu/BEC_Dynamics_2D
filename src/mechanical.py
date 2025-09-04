import numpy as np
import cupy as cp

from constants import hbar, m, x0, t0, e0
from wf import Norm

def COM(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray
        ) -> tuple[cp.float32, cp.float32]:
    '''
    functionality:
        find the center of mass of the wavepacket
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
    output:
        (Cx, Cy): center of mass of the wavepacket, shape (2,) # μm
    '''
    return (
        cp.sum(X * cp.abs(psi)**2) / cp.sum(cp.abs(psi)**2),
        cp.sum(Y * cp.abs(psi)**2) / cp.sum(cp.abs(psi)**2)
    )

def Iz(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, dx:cp.float32, dy:cp.float32, 
       Num:cp.int32) -> cp.float32:
    '''
    functionality:
        get the moment of inertia around z axis
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        dx: grid spacing in x direction # μm
        dy: grid spacing in y direction # μm
        Num: number of atoms
    output:
        Iz: moment of inertia around z axis (kg*μm^2/ms)
    '''
    return Num * cp.sum(cp.abs(psi)**2 * (X**2 + Y**2)) * (dx*dy) * (m*x0**2/t0)

def Iz_c(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, dx:cp.float32, dy:cp.float32, 
         Num:cp.int32) -> cp.float32:
    '''
    functionality:
        get the moment of inertia around COM along z axis
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        dx: grid spacing in x direction # μm
        dy: grid spacing in y direction # μm
        Num: number of atoms
    output:
        Iz_c: moment of inertia around COM along z axis (kg*μm^2/ms)
    '''
    (Cx, Cy) = COM(psi, X, Y)
    return Num * cp.sum(cp.abs(psi)**2 * ((X-Cx)**2 + (Y-Cy)**2)) * (dx*dy) * (m*x0**2/t0)

def wo_COM(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, Kx:cp.ndarray, Ky:cp.ndarray, 
           dx:cp.float32, dy:cp.float32) -> cp.ndarray:
    '''
    functionality:
        get the wavepacket without COM motion by cleaning the plane wave component
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        Kx: kx coordinates meshgrid, shape (Ny, Nx) # μm^-1
        Ky: ky coordinates meshgrid, shape (Ny, Nx) # μm^-1
    output:
        psi1: wavefunction without COM motion, shape (Ny, Nx)
    '''
    for _ in range(8): # iterate a few times to converge
        Fx = cp.conj(psi) * cp.fft.ifft2(Kx * cp.fft.fft2(psi,axes=(1,)),axes=(1,)) * (hbar/m)
        Fy = -cp.conj(psi) * cp.fft.ifft2(Ky * cp.fft.fft2(psi,axes=(0,)),axes=(0,)) * (hbar/m)
        Px = cp.sum(Fx).real * (dx*dy) * m # total momentum in x direction (kg*μm/ms)
        Py = cp.sum(Fy).real * (dx*dy) * m # total momentum in y direction (kg*μm/ms)
        psi = psi * cp.exp(-1j * (Px*X + Py*Y) / hbar)
    return psi

def flow_field(psi:cp.ndarray, psi1:cp.ndarray, Kx:cp.ndarray, Ky:cp.ndarray, 
               ) -> tuple[cp.ndarray, cp.ndarray, cp.ndarray, cp.ndarray]:
    '''
    functionality:
        get the flow field (velocity field times density field) of the wavepacket
    input:
        psi: wavefunction, shape (Ny, Nx)
        psi1: wavefunction without COM motion, shape (Ny, Nx)
        Kx: kx coordinates meshgrid, shape (Ny, Nx) # μm^-1
        Ky: ky coordinates meshgrid, shape (Ny, Nx) # μm^-1
    output:
        Fx: x component of the flow field, shape (Ny, Nx) # (μm*ms)^(-1)
        Fy: y component of the flow field, shape (Ny, Nx) # (μm*ms)^(-1)
        Fx1: x component of the flow field without COM motion, shape (Ny, Nx) # (μm*ms)^(-1)
        Fy1: y component of the flow field without COM motion, shape (Ny, Nx) # (μm*ms)^(-1)
    '''
    Fx = cp.conj(psi) * cp.fft.ifft2(Kx * cp.fft.fft2(psi,axes=(1,)),axes=(1,)) * (hbar/m)
    Fy = -cp.conj(psi) * cp.fft.ifft2(Ky * cp.fft.fft2(psi,axes=(0,)),axes=(0,)) * (hbar/m)
    Fx1 = cp.conj(psi1) * cp.fft.ifft2(Kx * cp.fft.fft2(psi1,axes=(1,)),axes=(1,)) * (hbar/m)
    Fy1 = -cp.conj(psi1) * cp.fft.ifft2(Ky * cp.fft.fft2(psi1,axes=(0,)),axes=(0,)) * (hbar/m)
    
    return (cp.real(Fx), cp.real(Fy), cp.real(Fx1), cp.real(Fy1))

def rotate(Fx:cp.ndarray, Fy:cp.ndarray, Fx1:cp.ndarray, Fy1:cp.ndarray, 
           dx:cp.float32, dy:cp.float32, X:cp.ndarray, Y:cp.ndarray, Num:cp.int32, 
           Iz_tot:cp.float32, Iz_sr:cp.float32
           ) -> tuple[cp.float32, cp.float32, cp.float32, cp.float32]:
    '''
    functionality:
        get the angular momentum and angular velocity, for self-rotation and total-rotation
    input:
        psi: wavefunction, shape (Ny, Nx)
        Fx: x component of the flow field, shape (Ny, Nx) # (μm*ms)^(-1)
        Fy: y component of the flow field, shape (Ny, Nx) # (μm*ms)^(-1)
        Fx1: x component of the flow field without COM motion, shape (Ny, Nx) # (μm*ms)^(-1)
        Fy1: y component of the flow field without COM motion, shape (Ny, Nx) # (μm*ms)^(-1)
    output:
        (Lz_tot, Lz_sr): total angular momentum and self-rotation angular momentum # ???
        (omega_tot, omega_sr): total angular velocity and self-rotation angular velocity # ms^-1
    '''
    Lz_tot = m * cp.sum(X*Fy - Y*Fx) * (dx*dy) * (m/hbar) * Num * x0**2
    Lz_sr = m * cp.sum(X*Fy1 - Y*Fx1)* (dx*dy) * (m/hbar) * Num * x0**2
    omega_tot = Lz_tot / Iz_tot / t0 if Iz_tot != 0 else cp.float32(0)
    omega_sr = Lz_sr / Iz_sr / t0 if Iz_sr != 0 else cp.float32(0)
    return (Lz_tot, Lz_sr, omega_tot, omega_sr)

def eigenaxis_angle(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray,
                    ) -> tuple[cp.float32, cp.float32]:
    '''
    functionality:
        find the polar angle of eigenaxis of the wavepacket's inertia tensor
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
    output:
        ang_l: the polar angle of the eigenaxis whose eigenvalue is larger # rad
        ang_s: the polar angle of the eigenaxis whose eigenvalue is smaller # rad
    '''
    ...