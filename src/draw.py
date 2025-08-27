import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from skimage.measure import block_reduce

def camera(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, colormap:str, 
           xlabel:str, ylabel:str, title:str, fontsize:float, file_name:str) -> None:
    '''
    functionality:
        plot the density distribution of the wavepacket
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        colormap: colormap for the density plot
    '''
    xmin = cp.asnumpy(cp.min(X))
    xmax = cp.asnumpy(cp.max(X))
    ymin = cp.asnumpy(cp.min(Y))
    ymax = cp.asnumpy(cp.max(Y))
    plt.imshow(cp.asnumpy(cp.abs(psi)**2), extent=(xmin,xmax,ymin,ymax), cmap=colormap)
    
    ax = plt.gca()
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_aspect('equal')

    plt.savefig(file_name, dpi=1200)
    plt.clf()

    return

def flow(Fx:cp.ndarray, Fy:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, 
         color:str, width:float, xlabel:str, ylabel:str, title:str, fontsize:float, reduce_exponent:int, file_name:str) -> None:
    '''
    '''
    n = 2**reduce_exponent
    X_reduced = block_reduce(cp.asnumpy(X), block_size=(n,n), func=np.mean)
    Y_reduced = block_reduce(cp.asnumpy(Y), block_size=(n,n), func=np.mean)
    Fx_reduced = block_reduce(cp.asnumpy(Fx), block_size=(n,n), func=np.mean)
    Fy_reduced = block_reduce(cp.asnumpy(Fy), block_size=(n,n), func=np.mean)
    dx = cp.asnumpy(X[0,1]-X[0,0])

    # to avoid overcrowding
    F2_max = np.max(Fx_reduced**2 + Fy_reduced**2)
    plt.quiver(X_reduced, Y_reduced, Fx_reduced, Fy_reduced, angles='xy', color=color, width=width, pivot='mid', scale=F2_max**0.5*550/n)

    ax = plt.gca()
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_aspect('equal')

    plt.savefig(file_name, dpi=1200)
    plt.clf()

    return