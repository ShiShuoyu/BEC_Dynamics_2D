import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import cupy as cp
import tkinter as tk
import json
from argparse import ArgumentParser
from tqdm import tqdm

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
    constants["Ar"] = int(args.relative_atomic_mass)
    constants["a"] = float(args.scattering_length)
    with open("constants.json", "w") as f:
        json.dump(constants, f, indent=4, ensure_ascii=False)

    # Time dimensionless
    from constants import t0, hbar, m, x0, e0, a, ab
    args.omega_trap = args.omega_trap / t0
    args.omega_bec = args.omega_bec / t0
    args.omega_trap_z = args.omega_trap_z / t0

    # interacting strength #
    g = (4*cp.pi*args.atom_number * hbar**2 * a*ab / m / e0 / x0**3
         ) * cp.sqrt(m*args.omega_trap_z/2/cp.pi/hbar * x0
                     ) # correction for squeezing a 3D wavepacket (along z direction) into a x0 μm layer

    # Generate the time structure
    (dt, duration, n_steps, n_samples) = ev.time_step(dt=args.dt, duration=args.duration, sampling_interval=args.sampling_interval)

    # Generate the grid, operators, and initial wavefunction
    (X, Y, Kx, Ky, dx, dy) = wf.grid(x_range=(-args.radius_xy[0], args.radius_xy[0]),
                                      y_range=(-args.radius_xy[1], args.radius_xy[1]),
                                        Nx=args.number_xy[0], Ny=args.number_xy[1])
    (U, V_sqrt, T) = wf.operator(X=X, Y=Y, Kx=Kx, Ky=Ky, omega=args.omega_trap,
                                 trap_center=np.array(args.center_trap), beta=args.beta, 
                                 r_0=args.r_0, imaginary_time=args.imaginary_time, dt=args.dt)
    psi = wf.wf_Gaussian(X=X, Y=Y, BEC_center=np.array(args.center_bec), omega=args.omega_bec, l=args.angular_momentum_bec[0], lz=args.angular_momentum_bec[1], dx=dx, dy=dy)
    psi = wf.boost(psi=psi, X=X, Y=Y, vx=args.velocity[0], vy=args.velocity[1])

    # Prepare the output arrays
    n_steps = n_steps + 1 # include the last step
    time = 0 # ms
    idx_sampling = 0 # index for sampling
    time_list = cp.zeros(n_samples, dtype=cp.float32)
    Iz_tot_list = cp.zeros(n_samples, dtype=cp.float32)
    Iz_sr_list = cp.zeros(n_samples, dtype=cp.float32)
    Lz_tot_list = cp.zeros(n_samples, dtype=cp.float32)
    Lz_sr_list = cp.zeros(n_samples, dtype=cp.float32)
    omega_tot_list = cp.zeros(n_samples, dtype=cp.float32)
    omega_sr_list = cp.zeros(n_samples, dtype=cp.float32)
    ang_c_list = cp.zeros(n_samples, dtype=cp.float32)
    ang_l_list = cp.zeros(n_samples, dtype=cp.float32)
    ang_s_list = cp.zeros(n_samples, dtype=cp.float32)

    if args.video and args.mechanics:
        print("Error: --video and --mechanics cannot be set at the same time.")
        return
    # time evolution loop
    if args.video:
        FFMpegWriter = animation.writers['ffmpeg']
        metadata = dict(title='BEC_Dynamics_2D', artist='matplotlib',
                        comment="Split step method")
        writer = FFMpegWriter(fps=25, codec='h264_nvenc', metadata=metadata,  extra_args=[
                '-b:v', '2000k', # bitrate
                '-crf', '20'
            ])

        fig, ax = plt.subplots()

        with writer.saving(fig, 'output/BEC_2D.mp4', dpi=150):
            # Time evolution loop
            for step in tqdm(range(n_steps)):
                # Sample the wavefunction
                if step % (args.sampling_interval) == 0:
                    # Store the wavefunction for visualization
                    draw.camera_video(Z=(cp.abs(psi)**2)*(args.atom_number), X=X, Y=Y, colormap='plasma', xlabel='x (μm)', ylabel='y (μm)', title=f'time = {time:.2f}ms', fontsize=16)
                    writer.grab_frame()

                # Evolve the wavefunction
                psi = ev.time_evolution(psi=psi, U=U, V_sqrt=V_sqrt, T=T, dt=args.dt, g=g, imaginary_time=args.imaginary_time)
                time = time + dt*t0 # physical time
    elif args.mechanics:
        for step in tqdm(range(n_steps)):
            if step % (args.sampling_interval) == 0:
                # Mechanical quantities
                psi1 = mc.wo_COM(psi=psi, X=X, Y=Y, Kx=Kx, Ky=Ky, dx=dx, dy=dy)
                Iz_tot = mc.Iz(psi=psi, X=X, Y=Y, dx=dx, dy=dy, Num=args.atom_number)
                Iz_sr = mc.Iz_c(psi=psi, X=X, Y=Y, dx=dx, dy=dy, Num=args.atom_number)
                (Fx,Fy,Fx1,Fy1) = mc.flow_field(psi=psi, psi1 = psi1, Kx=Kx, Ky=Ky)
                (Lz_tot, Lz_sr, omega_tot, omega_sr) = mc.rotate(Fx=Fx, Fy=Fy, Fx1=Fx1, Fy1=Fy1, dx=dx, dy=dy, X=X, Y=Y, Num=args.atom_number, Iz_tot=Iz_tot, Iz_sr=Iz_sr)
                (ang_c, ang_l, ang_s) = mc.eigenaxis_angle(psi=psi, X=X, Y=Y)
                # Store the mechanical quantities
                time_list[idx_sampling] = time
                Iz_tot_list[idx_sampling] = Iz_tot
                Iz_sr_list[idx_sampling] = Iz_sr
                Lz_tot_list[idx_sampling] = Lz_tot
                Lz_sr_list[idx_sampling] = Lz_sr
                omega_tot_list[idx_sampling] = omega_tot
                omega_sr_list[idx_sampling] = omega_sr
                ang_c_list[idx_sampling] = ang_c
                ang_l_list[idx_sampling] = ang_l
                ang_s_list[idx_sampling] = ang_s
                idx_sampling = idx_sampling + 1
            # Evolve the wavefunction
            psi = ev.time_evolution(psi=psi, U=U, V_sqrt=V_sqrt, T=T, dt=args.dt, g=g, imaginary_time=args.imaginary_time)
            time = time + dt*t0
    else:
        for step in tqdm(range(n_steps)):
            # Evolve the wavefunction
            psi = ev.time_evolution(psi=psi, U=U, V_sqrt=V_sqrt, T=T, dt=args.dt, g=g, imaginary_time=args.imaginary_time)
            time = time + dt*t0
    if args.figure: # output the density profile and flow field of the final state
        print('\nsaving figures ...')
        # density profile - save with both time-specific and consistent names
        draw.camera(Z=(cp.abs(psi)**2)*(args.atom_number), X=X, Y=Y, colormap='plasma', xlabel='x (μm)', ylabel='y (μm)', title=f'time = {time:.2f}ms', fontsize=16, file_name=f'output/density_t{time:.0f}ms.png')
        draw.camera(Z=(cp.abs(psi)**2)*(args.atom_number), X=X, Y=Y, colormap='plasma', xlabel='x (μm)', ylabel='y (μm)', title=f'time = {time:.2f}ms', fontsize=16, file_name='output/density_final.png')  # consistent name
        
        draw.camera(Z=(cp.real(psi)**2)*(args.atom_number), X=X, Y=Y, colormap='hot', xlabel='x (μm)', ylabel='y (μm)', title=f'time = {time:.2f}ms', fontsize=16, file_name=f'output/real_t{time:.0f}ms.png')
        draw.camera(Z=(cp.real(psi)**2)*(args.atom_number), X=X, Y=Y, colormap='hot', xlabel='x (μm)', ylabel='y (μm)', title=f'time = {time:.2f}ms', fontsize=16, file_name='output/real_final.png')  # consistent name
        
        draw.camera(Z=cp.angle(psi), X=X, Y=Y, colormap='viridis', xlabel='x (μm)', ylabel='y (μm)', title=f'time = {time:.2f}ms', fontsize=16, file_name=f'output/phase_t{time:.0f}ms.png')
        draw.camera(Z=cp.angle(psi), X=X, Y=Y, colormap='viridis', xlabel='x (μm)', ylabel='y (μm)', title=f'time = {time:.2f}ms', fontsize=16, file_name='output/phase_final.png')  # consistent name

        # flow field - save with both time-specific and consistent names
        psi1 = mc.wo_COM(psi=psi, X=X, Y=Y, Kx=Kx, Ky=Ky, dx=dx, dy=dy)
        (Fx,Fy,Fx1,Fy1) = mc.flow_field(psi=psi, psi1=psi1, Kx=Kx, Ky=Ky)
        draw.flow(Fx=Fx, Fy=Fy, X=X, Y=Y, color='blue', width=0.001, xlabel='x (μm)', ylabel='y (μm)', title=f'flow field at t = {time:.2f}ms', fontsize=16, reduce_exponent=3, file_name=f'output/flow_t{time:.0f}ms.png')
        draw.flow(Fx=Fx, Fy=Fy, X=X, Y=Y, color='blue', width=0.001, xlabel='x (μm)', ylabel='y (μm)', title=f'flow field at t = {time:.2f}ms', fontsize=16, reduce_exponent=3, file_name='output/flow_final.png')  # consistent name
        
        draw.flow(Fx=Fx1, Fy=Fy1, X=X, Y=Y, color='blue', width=0.001, xlabel='x (μm)', ylabel='y (μm)', title=f'flow field w/o COM at t = {time:.2f}ms', fontsize=16, reduce_exponent=3, file_name=f'output/flow1_t{time:.0f}ms.png')
        draw.flow(Fx=Fx1, Fy=Fy1, X=X, Y=Y, color='blue', width=0.001, xlabel='x (μm)', ylabel='y (μm)', title=f'flow field w/o COM at t = {time:.2f}ms', fontsize=16, reduce_exponent=3, file_name='output/flow1_final.png')  # consistent name
    if args.mechanics: # output the time evolution of mechanical quantities
        print('\nsaving mechanical quantities ...')
        draw.quantity(time=time_list, quantity=Iz_tot_list, xlabel='time (ms)', ylabel='Iz_tot (kg*μm^2)', title='Moment of inertia around z', fontsize=16, file_name='output/Iz_tot.png')
        draw.quantity(time=time_list, quantity=Iz_sr_list, xlabel='time (ms)', ylabel='Iz_sr (kg*μm^2)', title='Moment of inertia around COM', fontsize=16, file_name='output/Iz_sr.png')
        draw.quantity(time=time_list, quantity=Lz_tot_list, xlabel='time (ms)', ylabel='Lz_tot (kg*μm^2/ms)', title='Total angular momentum', fontsize=16, file_name='output/Lz_tot.png')
        draw.quantity(time=time_list, quantity=Lz_sr_list, xlabel='time (ms)', ylabel='Lz_sr (kg*μm^2/ms)', title='Intrinsic angular momentum', fontsize=16, file_name='output/Lz_sr.png')
        draw.quantity(time=time_list, quantity=omega_tot_list, xlabel='time (ms)', ylabel='omega_tot (ms^-1)', title='Total angular velocity', fontsize=16, file_name='output/omega_tot.png')
        draw.quantity(time=time_list, quantity=omega_sr_list, xlabel='time (ms)', ylabel='omega_sr (ms^-1)', title='Intrinsic angular velocity', fontsize=16, file_name='output/omega_sr.png')
        draw.quantity(time=time_list, quantity=ang_c_list, xlabel='time (ms)', ylabel='ang_com (rad)', title='Polar angle of COM', fontsize=16, file_name='output/ang_com.png')
        draw.quantity(time=time_list, quantity=ang_l_list, xlabel='time (ms)', ylabel='ang_major (rad)', title='Polar angle of major axis', fontsize=16, file_name='output/ang_l.png')
        draw.quantity(time=time_list, quantity=ang_s_list, xlabel='time (ms)', ylabel='ang_minor (rad)', title='Polar angle of minor axis', fontsize=16, file_name='output/ang_s.png')

        mc.save_quantities(time_list=time_list, Iz_tot_list=Iz_tot_list, Iz_sr_list=Iz_sr_list, Lz_tot_list=Lz_tot_list, Lz_sr_list=Lz_sr_list, omega_sr_list=omega_sr_list, omega_tot_list=omega_tot_list, ang_c_list=ang_c_list, ang_l_list=ang_l_list, ang_s_list=ang_s_list)
    print('\nDONE!')
    return

if __name__ == "__main__":
    main()