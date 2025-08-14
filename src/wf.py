import numpy as np
import cupy as cp

def grid(x_range:list, y_range:list, Nx:int, Ny:int
                 ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    '''
    functionality: 
        generate a series of 2D grid, including X, Y in real space and Kx, Ky in momentum space
    input:
        x_range: (xmin, xmax) # μm
        y_range: (ymin, ymax) # μm
        Nx: number of points in x direction
        Ny: number of points in y direction
    output:
        vx: x coordinates, shape (Nx,)
        vy: y coordinates, shape (Ny,)
        [X, Y]: meshgrid of x and y coordinates, shape (Ny, Nx)
    '''
    vx = np.linspace(x_range[0], x_range[1], Nx)
    vy = np.linspace(y_range[0], y_range[1], Ny)
    [X, Y] = np.meshgrid(vx,vy)
    dx = vx[1] - vx[0]
    dy = vy[1] - vy[0]
    kx = np.fft.fftshift(np.linspace(-np.pi/dx, np.pi/dx, Nx))
    ky = np.fft.fftshift(np.linspace(-np.pi/dy, np.pi/dy, Ny))
    [Kx, Ky] = np.meshgrid(kx, ky)
    return (X, Y, Kx, Ky)

def operator(X:np.ndarray, Y:np.ndarray, Kx:np.ndarray, Ky:np.ndarray, 
             m:np.float32, omega_0:np.float32, trap_center:np.ndarray, beta:np.float32, 
             r_0:np.float32, imaginary_time:np.bool) -> tuple[cp.ndarray1, cp.ndarray]:
    '''
    functionality: 
        generate a series of 2D grid, V the potential energy operator and T the kinetic energy operator
    input:
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        Kx: kx coordinates meshgrid, shape (Ny, Nx) # μm^-1
        Ky: ky coordinates meshgrid, shape (Ny, Nx) # μm^-1
        m: mass of the particle # kg
        omega_0: trap frequency # ms^-1
        trap_center: center of the trap, shape (2,) # μm
        beta: anisotropy parameter # dimensionless
        r_0: a legnth scale, to keep beta dimensionless # μm
        imaginary_time: if True, perform imaginary time evolution
        
    output:
        V: potential energy operator, shape (Ny, Nx) # dimensionless
        T: kinetic energy operator, shape (Ny, Nx) # dimensionless
    '''
    loss = 1 - 0.4j if imaginary_time else 1
    T = cp.asarray(
        np.exp(-1j * loss * 
        (Kx**2 + Ky**2) / 2
            ), dtype=cp.complex32
                )
    V = cp.asarray(
        np.exp(-1j * loss * 
        0.5 * m * omega_0**2 * (((1-2*beta) * ((X-trap_center[0])**2+(Y-trap_center[1])**2)) + (beta/r_0) * ((X-trap_center[0])**2+(Y-trap_center[1])**2)**2)
            ), dtype=cp.complex32
                )
    return (V, T)
    ...

def wavefunction():
    ...