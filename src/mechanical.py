import numpy as np
import cupy as cp

def COM(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray) -> tuple[np.float32, np.float32]:
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
        cp.asnumpy(cp.sum(X * cp.abs(psi)**2) / cp.sum(cp.abs(psi)**2)),
        cp.asnumpy(cp.sum(Y * cp.abs(psi)**2) / cp.sum(cp.abs(psi)**2))
    )

def Iz(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, dx:np.float32, dy:np.float32) -> np.float32:
    '''
    functionality:
        get the moment of inertia around z axis
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        dx: grid spacing in x direction # μm
        dy: grid spacing in y direction # μm
    output:
        Iz: moment of inertia around z axis
    '''
    return cp.sum(cp.abs(psi)**2 * (X**2 + Y**2)) * (dx*dy)

def Iz_c(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, dx:np.float32, dy:np.float32) -> np.float32:
    '''
    functionality:
        get the moment of inertia around COM along z axis
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        dx: grid spacing in x direction # μm
        dy: grid spacing in y direction # μm
    output:
        Iz_c: moment of inertia around COM along z axis
    '''
    (Cx, Cy) = COM(psi, X, Y)
    return cp.sum(cp.abs(psi)**2 * ((X-Cx)**2 + (Y-Cy)**2)) * (dx*dy)
    ...

def wo_COM(ps:cp.ndarrayi, X:cp.ndarray, Y:cp.ndarray, Kx:cp.ndarray, Ky:cp.ndarray, Nx:np.int32, Ny:np.int32):
    '''
    functionality:
        get the wavepacket without COM motion by cleaning the plane wave component
    input:
    output:
    '''
    ...

def flow_fidle(psi:cp.ndarray, psi1:cp.ndarray, Kx:cp.ndarray, Ky:cp.ndarray, Nx:np.int32, Ny:np.int32):
    '''
    functionality:
    input:
    output:
    '''
    ...

def rotate(psi:cp.ndarray, Fx:cp.ndarray, Fy:cp.ndarray, Fx1:cp.ndarray, Fy1:cp.ndarray, X:cp.ndarray, Y:cp.ndarray):
    '''
    functionality:
    input:
    output:
    '''
    ...

def eigenaxis(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray,):
    '''
    functionality:
    input:
    output:
    '''
    ...