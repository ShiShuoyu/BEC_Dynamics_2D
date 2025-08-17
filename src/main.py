import numpy as np
import matplotlib.pyplot as plt
import cupy as cp
import tkinter as tk
import argparse

from constants import *
import wf
import evolve

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description="Simulating and visualizing the machanical motion of 2D Bose-Einstein condensates (BECs)"
        )
    parser.add_argument('--duration', type=np.float32, default=200, 
                        help='Totle duration of simulation in ms')
    parser.add_argument('--dt', type=np.float32, default=0.2, 
                        help='Time step in ms')
    parser.add_argument('--sampling_interval', type=np.int32, default=20, 
                        help='After how many time steps to sample the wave function')
    parser.add_argument('--radius_xy', type=np.float32, nargs=2, default=(50, 50), 
                        help='Radius of the real space simulation region in μm (e.g. ''(50, 40)'' means -50 to 50 μm in x and -40 to 40 μm in y)')
    parser.add_argument('--number_xy', type=np.int32, nargs=2, default=(256, 256),
                        help='Number of points in x and y directions (e.g. ''(256, 256)'' means 256 points in both x and y directions) (When N is a power of 2, the FFT is faster)')
    parser.add_argument('--omega_trap', type=np.float32, default=0.01*np.pi, 
                        help='Trapping frequency in ms^-1')
    parser.add_argument('--center_trap', type=np.float32, nargs=2, default=(0, 0), 
                        help='Center of the trap in μm (e.g. ''(0, 0)'' means the center is at (0, 0) μm)')
    parser.add_argument('--beta', type=np.float32, default=-0.1, 
                        help='Anharmonic parameter (dimensionless) (V = 0.5*m*omega^2*((1-2*beta)*(x^2+y^2) + (beta/r_0)*(x^2+y^2)^2)), where r_0 is the length scale to keep beta dimensionless)')
    parser.add_argument('--center_bec', type=np.float32, nargs=2, default=(30, 0), 
                        help='Center of the BEC in μm (e.g. ''25, 25'' means the center is at (25, 25) μm, and r_0 will be set to sqrt(2)*25 μm)')
    parser.add_argument('--omega_bec', type=np.float32, default=0.01*np.pi, 
                        help='The frequency of the harmonic trap whose ground state is this gaussian wavepacket')
    parser.add_argument('--imaginary_time', action='store_true', 
                        help='If set, perform imaginary time evolution')
    parser.add_argument('--velocity', type=np.float32, nargs=2, default=(0, 0),
                        help='Velocity of the initial wavepacket in μm/ms (e.g. ''(0, 0)'' means no initial velocity)')
    parser.add_argument('--video', action='store_true', 
                        help='If set, save the simulation results as a video')
    args = parser.parse_args()
    print("parameters:")
    print(f"duration: {args.duration} ms")
    print(f"time step: {args.dt} ms")
    print(f"sampling interval: {args.sampling_interval}")
    print(f"radius in x: {args.radius_xy[0]} μm, radius in y: {args.radius_xy[1]} μm")
    print(f"number in x: {args.output[0]}, number in y: {args.output[1]}")
    print(f"trapping frequency ω: {args.omega} ms^-1")
    print(f"center of the trap: ({args.center_trap[0]}, {args.center_trap[1]}) μm")
    print(f"anharmonic parameter β: {args.beta}")
    print(f"center of the BEC: ({args.center_bec[0]}, {args.center_bec[1]}) μm")
    print(f"imaginary time evolution: {'on' if args.imaginary_time else 'off'}")
    print(f"initial velocity: ({args.velocity[0]}, {args.velocity[1]}) μm/ms")

    # Generate the grid, operators, and initial wavefunction
    (X, Y, Kx, Ky, dx, dy) = wf.grid(x_range=(-args.radius_xy[0], args.radius_xy[0]),
                                      y_range=(-args.radius_xy[1], args.radius_xy[1]),
                                        Nx=args.number_xy[0], Ny=args.number_xy[1])
    (U, V_sqrt, T) = wf.operator(X=X, Y=Y, Kx=Kx, Ky=Ky, m=m, omega=args.omega_trap,
                                 trap_center=np.array(args.center_trap), beta=args.beta)
    psi = wf.wf_Gaussian(X=X, Y=Y, BEC_center=np.array(args.center_bec), omega=args.omega_bec, dx=dx, dy=dy)
    psi = wf.boost(psi=psi, X=X, Y=Y, vx=args.velocity[0], vy=args.velocity[1])

    # Generate the time structure
    (dt, duration, n_steps, n_samples) = evolve.time_step(dt=args.dt, duration=args.duration, sampling_interval=args.sampling_interval)

    # Prepare the output arrays
    time = 0
    ...

    # Time evolution loop
    for step in range(n_steps):
        # Evolve the wavefunction
        psi = evolve.time_evolution(psi=psi, U=U, V_sqrt=V_sqrt, T=T, dt=dt, g=0)
        time = time + dt

        # Sample the wavefunction
        if step % (args.sampling_interval) == 0:
            if args.video:
                # Store the wavefunction for visualization
                ...
            # Store the mechanical quantities
            ...

if __name__ == "__main__":
    main()