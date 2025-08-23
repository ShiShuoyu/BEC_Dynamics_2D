import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import cupy as cp
import tkinter as tk
import json
from argparse import ArgumentParser

import arguments as agms
import wf
import evolve as ev
import mechanical as mc
import draw

def main():
    # parse command line arguments
    parser = agms.input_args()
    args = parser.parse_args()
    # display the simulation parameters
    agms.display_args(args)    

    # Update the constants.json file with the new atomic parameters
    with open("constants.json", "r") as f:
        constants = json.load(f)
    constants["Ar"] = args.relative_atomic_mass
    constants["a"] = args.scattering_length
    with open("constants.json", "w") as f:
        json.dump(constants, f, indent=4)

    # Generate the grid, operators, and initial wavefunction
    (X, Y, Kx, Ky, dx, dy) = wf.grid(x_range=(-args.radius_xy[0], args.radius_xy[0]),
                                      y_range=(-args.radius_xy[1], args.radius_xy[1]),
                                        Nx=args.number_xy[0], Ny=args.number_xy[1])
    (U, V_sqrt, T) = wf.operator(X=X, Y=Y, Kx=Kx, Ky=Ky, omega=args.omega_trap,
                                 trap_center=np.array(args.center_trap), beta=args.beta, 
                                 r_0=args.r_0, imaginary_time=args.imaginary_time, dt=args.dt)
    psi = wf.wf_Gaussian(X=X, Y=Y, BEC_center=np.array(args.center_bec), omega=args.omega_bec, dx=dx, dy=dy)
    psi = wf.boost(psi=psi, X=X, Y=Y, vx=args.velocity[0], vy=args.velocity[1])

    # Generate the time structure
    (dt, duration, n_steps, n_samples) = ev.time_step(dt=args.dt, duration=args.duration, sampling_interval=args.sampling_interval)

    # Prepare the output arrays
    time = 0
    ...
    # Prepare the output video
    if args.video:
        FFMpegWriter = animation.writers['ffmpeg']
        metadata = dict(title='BEC_Dynamics_2D', artist='matplotlib',
                        comment="Split step method")
        writer = FFMpegWriter(fps=30, metadata=metadata)

    fig, ax = plt.subplots()

    with writer.saving(fig, 'BEC_2D.mp4', 100):
        # Time evolution loop
        for step in range(n_steps):
            # Evolve the wavefunction
            psi = ev.time_evolution(psi=psi, U=U, V_sqrt=V_sqrt, T=T, dt=dt, Num=args.atom_number, omega_z=args.omega_trap_z)
            time = time + dt

            # Sample the wavefunction
            if step % (args.sampling_interval) == 0:
                if args.video:
                    # Store the wavefunction for visualization
                    draw.camera(psi=psi, X=X, Y=Y, colormap='hot', xlabel='x/μm', ylabel='y/μm', title=['time = ', str(time), 'ms'], fontsize=20)
                    writer.grab_frame()
                    # Store the mechanical quantities
                    ...

if __name__ == "__main__":
    main()