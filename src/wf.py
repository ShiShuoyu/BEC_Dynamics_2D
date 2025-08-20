import numpy as np
import cupy as cp

from constants import hbar, m, e0, t0

def grid(x_range:np.ndarray, y_range:np.ndarray, Nx:np.int32, Ny:np.int32
                 ) -> tuple[cp.ndarray, cp.ndarray, cp.ndarray, cp.ndarray, np.float32, np.float32]:
    '''
    functionality: 
        generate a series of 2D grid, including X, Y in real space and Kx, Ky in momentum space
    input:
        x_range: (xmin, xmax) # μm
        y_range: (ymin, ymax) # μm
        Nx: number of points in x direction
        Ny: number of points in y direction
    output:
        [X, Y]: meshgrid of x and y coordinates, shape (Ny, Nx)
        [Kx, Ky]: meshgrid and fftshift of kx and ky momentums, shape (Ny, Nx)
        [dx, dy]: grid spacing in x and y directions # μm
    '''
    vx = np.linspace(x_range[0], x_range[1], Nx)
    vy = np.linspace(y_range[0], y_range[1], Ny)
    [X, Y] = np.meshgrid(vx,vy)
    dx = vx[1] - vx[0]
    dy = vy[1] - vy[0]
    kx = np.fft.fftshift(np.linspace(-np.pi/dx, np.pi/dx, Nx))
    ky = np.fft.fftshift(np.linspace(-np.pi/dy, np.pi/dy, Ny))
    [Kx, Ky] = np.meshgrid(kx, ky)
    X, Y = cp.asarray(X, dtype=cp.float32), cp.asarray(Y, dtype=cp.float32)
    Kx, Ky = cp.asarray(Kx, dtype=cp.float32), cp.asarray(Ky, dtype=cp.float32)
    return (X, Y, Kx, Ky, dx, dy)

def operator(X:cp.ndarray, Y:cp.ndarray, Kx:cp.ndarray, Ky:cp.ndarray, 
             omega:np.float32, trap_center:np.ndarray, beta:np.float32, 
             r_0:np.float32, imaginary_time:np.bool, dt:np.float32) -> tuple[cp.ndarray, cp.ndarray]:
    '''
    functionality: 
        generate a series of 2D grid, V the potential energy operator and T the kinetic energy operator
    input:
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        Kx: kx coordinates meshgrid, shape (Ny, Nx) # μm^-1
        Ky: ky coordinates meshgrid, shape (Ny, Nx) # μm^-1
        omega: trap frequency # ms^-1
        trap_center: center of the trap, shape (2,) # μm
        beta: anisotropy parameter # dimensionless
        r_0: a legnth scale, to keep beta dimensionless # μm
        imaginary_time: if True, perform imaginary time evolution
        dt: time step # ms
    output:
        U: potential energy, shape (Ny, Nx) # dimensionless
        V_sqrt: potential energy operator over dt/2, shape (Ny, Nx) # dimensionless
        T: kinetic energy operator over dt, shape (Ny, Nx) # dimensionless
    '''
    loss = 1 - 0.4j if imaginary_time else 1

    U = cp.asarray(0.5 * m * omega**2 * (((1-2*beta) * ((X-trap_center[0])**2+(Y-trap_center[1])**2)) + (beta/r_0**2) * ((X-trap_center[0])**2+(Y-trap_center[1])**2)**2) / e0, dtype=cp.complex64)

    # pesudo potential to avoid the wavepacket to escape the trap during imaginary time evolution
    # in the region out of the trap, we set the minimal potential to 0.2 * U_max
    R_squared = (X - trap_center[0])**2 + (Y - trap_center[1])**2
    U_max = cp.max(cp.abs(U))
    max_index = cp.argmax(cp.abs(U))
    r2_threshold = R_squared.flatten()[max_index]
    condition = (R_squared > r2_threshold) & (cp.abs(U) < 0.2 * U_max)
    U = cp.where(condition, 0.2 * U_max, U)

    V_sqrt = cp.exp(-1j * loss * (dt/2) * U)

    T = cp.asarray(
        np.exp(-1j * loss * dt * 
        (Kx**2 + Ky**2) / 2
            ), dtype=cp.complex64
                ) # Here K.E. = p^2/2m / e0 = (hbar^2/2m)*k^2 / e0 = (1/2) * x0^2 * k^2 -> dimensionless
    return (U, V_sqrt, T)

def Norm(psi:cp.ndarray, dx:np.float32, dy:np.float32) -> np.float32:
    '''
    functionality:
        calculate the norm of the wavefunction
    input:
        psi: wavefunction, shape (Ny, Nx)
        dx: grid spacing in x direction # μm
        dy: grid spacing in y direction # μm
    output:
        N: norm of the wavefunction # dimensionless
    '''
    return cp.sum(cp.abs(psi)**2) * (dx*dy)

def wf_Gaussian(X:cp.ndarray, Y:cp.ndarray, BEC_center:np.ndarray, omega:np.float32, 
                dx:np.float32, dy:np.float32) -> cp.ndarray:
    '''
    functionality:
        generate a Gaussian wavefunction for the BEC
    input:
        X: x coordinates meshgrid, shape (Ny, Nx) # μm  
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        BEC_center: center of the BEC, shape (2,) # μm
        omega: the frequency of the harmonic trap whose ground state is this gaussian wavepacket # ms^-1
        dx: grid spacing in x direction # μm
        dy: grid spacing in y direction # μm
    output:
        psi: wavefunction, shape (Ny, Nx) # normalized
    '''
    psi = cp.asarray(
        np.exp(
            -(m*omega/2/hbar) * ((X-BEC_center[0])**2 + (Y-BEC_center[1])**2)
            ), dtype=cp.complex64
                    )
    return psi / np.sqrt(Norm(psi, dx, dy))

def wf_ThomasFermi():
    ...

def boost(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, vx:np.float32, vy:np.float32
          ) -> cp.ndarray:
    '''
    functionality:
        boost the wavepacket by implying a plane wave exp(1j*(kx*X + ky*Y))
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        (vx, vy): the initial velocity of the wavepacket # μm/ms
    output:
        psi: boosted wavefunction, shape (Ny, Nx)
    '''
    # psi = psi e^{i\frac{m\vec{v}}{\hbar}\cdot\vec{x}}
    return psi * cp.exp(1j * (vx*X + vy*Y) * t0)