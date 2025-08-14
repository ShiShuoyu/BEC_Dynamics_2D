import numpy as np
import cupy as cp

def time_evolution(psi:np.ndarray, V:np.ndarray, T:np.ndarray, dt:np.float32, 
                   imaginary_time:np.bool, g:np.float32) -> np.ndarray:
    '''
    input:
        psi: wave function, shape (N, N)
        V: potential energy operator, shape (N, N) 
        T: kinetic energy operator, shape (N, N)
        dt: time step
        imaginary_time: if True, perform imaginary time evolution
        g: interaction strength, scalar
    output:
        psi: evolved wave function, shape (N, N)
    '''
    ...
    return psi