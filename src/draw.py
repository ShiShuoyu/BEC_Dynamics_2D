import numpy as np
import cupy as cp
import matplotlib.pyplot as plt
from skimage.measure import block_reduce

def camera_video(psi:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, colormap:str, 
           xlabel:str, ylabel:str, title:str, fontsize:float) -> None:
    '''
    functionality:
        plot the density distribution of the wavepacket for video
    input:
        psi: wavefunction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        colormap: colormap for the density plot
        xlabel: label for x axis
        ylabel: label for y axis
        title: title of the figure
        fontsize: fontsize for the labels and title
    '''
    xmin = cp.asnumpy(cp.min(X))
    xmax = cp.asnumpy(cp.max(X))
    ymin = cp.asnumpy(cp.min(Y))
    ymax = cp.asnumpy(cp.max(Y))
    plt.imshow(cp.asnumpy(cp.abs(psi)**2), extent=(xmin,xmax,ymin,ymax), cmap=colormap)
    
    ax = plt.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_aspect('equal')

    return

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
        xlabel: label for x axis
        ylabel: label for y axis
        title: title of the figure
        fontsize: fontsize for the labels and title
        file_name: name of the output file
    '''
    xmin = cp.asnumpy(cp.min(X))
    xmax = cp.asnumpy(cp.max(X))
    ymin = cp.asnumpy(cp.min(Y))
    ymax = cp.asnumpy(cp.max(Y))
    plt.imshow(cp.asnumpy(cp.abs(psi)**2), extent=(xmin,xmax,ymin,ymax), cmap=colormap)
    
    ax = plt.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_aspect('equal')

    plt.savefig(file_name, dpi=600)
    plt.clf()

    return

def flow(Fx:cp.ndarray, Fy:cp.ndarray, X:cp.ndarray, Y:cp.ndarray, 
         color:str, width:float, xlabel:str, ylabel:str, title:str, fontsize:float, reduce_exponent:int, file_name:str) -> None:
    '''
    functionality:
        plot the flow field of the wavepacket
    input:
        Fx: flow field in x direction, shape (Ny, Nx)
        Fy: flow field in y direction, shape (Ny, Nx)
        X: x coordinates meshgrid, shape (Ny, Nx) # μm
        Y: y coordinates meshgrid, shape (Ny, Nx) # μm
        color: color of the arrows
        width: width of the arrows
        xlabel: label for x axis
        ylabel: label for y axis
        title: title of the figure
        fontsize: fontsize for the labels and title
        reduce_exponent: reduce the number of arrows by a factor of 2**reduce_exponent
        file_name: name of the output file
    '''
    xmin = cp.asnumpy(cp.min(X))
    xmax = cp.asnumpy(cp.max(X))
    ymin = cp.asnumpy(cp.min(Y))
    ymax = cp.asnumpy(cp.max(Y))
    n = 2**reduce_exponent
    X_reduced = block_reduce(cp.asnumpy(X), block_size=(n,n), func=np.mean)
    Y_reduced = block_reduce(cp.asnumpy(Y), block_size=(n,n), func=np.mean)
    Fx_reduced = block_reduce(cp.asnumpy(Fx), block_size=(n,n), func=np.mean)
    Fy_reduced = block_reduce(cp.asnumpy(Fy), block_size=(n,n), func=np.mean)
    dx = float(X[0,1] - X[0,0])

    # to avoid overcrowding
    F2_max = np.max(Fx_reduced**2 + Fy_reduced**2)
    plt.quiver(X_reduced, Y_reduced, Fx_reduced, Fy_reduced, angles='xy', color=color, width=width, pivot='mid', scale=F2_max**0.5*100/n/dx)

    ax = plt.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)
    ax.set_aspect('equal')

    plt.savefig(file_name, dpi=600)
    plt.clf()

    return

def combine():
    ...

def quantity(time:cp.ndarray, quantity:cp.ndarray, xlabel:str, ylabel:str, title:str, fontsize:float, 
             file_name:str) -> None:
    '''
    functionality:
        plot the time variation of physical quantities
    input:
        time: time of sampling points # ms
        quantity: physical quantity at sampling points # units
        xlabel: label for x axis
        ylabel: label for y axis
        title: title of the figure
        fontsize: fontsize for the labels and title
        file_name: file name of the output figure
    '''
    xmin = cp.asnumpy(cp.min(time) - 0.05*(cp.max(time)-cp.min(time)))
    xmax = cp.asnumpy(cp.max(time) + 0.05*(cp.max(time)-cp.min(time)))
    ymin = cp.asnumpy(cp.min(quantity) - 0.05*(cp.max(quantity)-cp.min(quantity)))
    ymax = cp.asnumpy(cp.max(quantity) + 0.05*(cp.max(quantity)-cp.min(quantity)))

    plt.plot(cp.asnumpy(time), cp.asnumpy(quantity), linewidth=2)
    ax = plt.gca()
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.set_title(title, fontsize=fontsize)

    plt.savefig(file_name, dpi=600)
    plt.clf()
    return