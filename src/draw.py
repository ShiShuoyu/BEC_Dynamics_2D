import numpy as np
import cupy as cp
import matplotlib.pyplot as plt

def camera(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, colormap:str, 
           xlabel:str, ylabel:str, title:str, fontsize:np.float32) -> None:
    '''
    functionality:
        plot the density distribution of the wavepacket
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        colormap: colormap for the density plot
    '''
    psi = cp.asnumpy(psi)
    xmin = cp.asnumpy(cp.min(X))
    xmax = cp.asnumpy(cp.max(X))
    ymin = cp.asnumpy(cp.min(Y))
    ymax = cp.asnumpy(cp.max(Y))
    plt.imshow(np.abs(psi)**2, extent=(xmin,xmax,ymin,ymax), cmap=colormap, xlabel=xlabel, ylabel=ylabel, title=title, fontsize=fontsize)
    plt.show()
    return

def flow(Fx:cp.ndarray, Fy:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, 
         color:str, width:np.float32) -> None:
    Fx = cp.asnumpy(Fx)
    Fy = cp.asnumpy(Fy)
    X = cp.asnumpy(X)
    Y = cp.asnumpy(Y)
    plt.quiver(X, Y, Fx, Fy, scale=5, angles='xy', color=color, width=width)
    return