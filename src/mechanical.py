import numpy as np
import cupy as cp

def COM(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray
        ) -> tuple[np.float32, np.float32]:
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

def Iz(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, dx:np.float32, dy:np.float32
       ) -> np.float32:
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

def Iz_c(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, dx:np.float32, dy:np.float32
         ) -> np.float32:
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

def wo_COM(ps:cp.ndarrayi, X:cp.ndarray, Y:cp.ndarray, 
           Kx:cp.ndarray, Ky:cp.ndarray, Nx:np.int32, Ny:np.int32) -> cp.ndarray:
    '''
    functionality:
        get the wavepacket without COM motion by cleaning the plane wave component
    input:
    output:
        psi1: wavefunction without COM motion, shape (Ny, Nx)
    '''
    ...

def flow_fidle(psi:cp.ndarray, psi1:cp.ndarray, Kx:cp.ndarray, Ky:cp.ndarray, 
               Nx:np.int32, Ny:np.int32
               ) -> tuple[cp.ndarray, cp.ndarray, cp.ndarray, cp.ndarray]:
    '''
    functionality:
        get the flow field (velocity field times density field) of the wavepacket
    input:
    output:
        Fx: x component of the flow field, shape (Ny, Nx)
        Fy: y component of the flow field, shape (Ny, Nx)
        Fx1: x component of the flow field without COM motion, shape (Ny, Nx)
        Fy1: y component of the flow field without COM motion, shape (Ny, Nx)
    '''
    ...

def rotate(psi:cp.ndarray, Fx:cp.ndarray, Fy:cp.ndarray, 
           Fx1:cp.ndarray, Fy1:cp.ndarray, X:cp.ndarray, Y:cp.ndarray
           ) -> tuple[np.float32, np.float32, np.float32, np.float32]:
    '''
    functionality:
        get the angular momentum and angular velocity, for self-rotation and total-rotation
    input:
        psi: wavefunction, shape (Ny, Nx)
    output:
        (Lz_tot, Lz_sr): total angular momentum and self-rotation angular momentum # ???
        (w_tot, w_sr): total angular velocity and self-rotation angular velocity # ???
    '''
    ...

def eigenaxis_angle(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray,
                    ) -> tuple[np.float32, np.float32]:
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