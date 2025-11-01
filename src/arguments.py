from argparse import ArgumentParser
import numpy as np
import cupy as cp

def input_args() -> ArgumentParser:
    '''
    functionality:
        parse the command line arguments
    output:
        parser: ArgumentParser object
    '''
    parser = ArgumentParser(
        description="Simulating and visualizing the machanical motion of 2D Bose-Einstein condensates (BECs)"
        )
    parser.add_argument('--relative_atomic_mass', type=cp.int32, default=87, 
                        help='Atomic mass of the particle in amu (e.g. ''87'' means Rb-87)')
    parser.add_argument('--scattering_length', type=cp.float32, default=10,
                        help='Scattering length in Bohr radius (e.g. ''10'' means a_s = 10 a_Bohr)')
    parser.add_argument('--atom_number', type=cp.int32, default=2500, 
                        help='Number of atoms')
    parser.add_argument('--duration', type=cp.float32, default=200, 
                        help='Totle duration of simulation in ms')
    parser.add_argument('--dt', type=cp.float32, default=0.5, 
                        help='Time step in ms')
    parser.add_argument('--sampling_interval', type=cp.int32, default=10, 
                        help='After how many time steps to sample the wave function')
    parser.add_argument('--radius_xy', type=cp.float32, nargs=2, default=(50, 50), 
                        help='Radius of the real space simulation region in μm (e.g. ''(50, 40)'' means -50 to 50 μm in x and -40 to 40 μm in y)')
    parser.add_argument('--number_xy', type=cp.int32, nargs=2, default=(256, 256),
                        help='Number of points in x and y directions (e.g. ''(256, 256)'' means 256 points in both x and y directions) (When N is a power of 2, the FFT is faster)')
    parser.add_argument('--omega_trap', type=cp.float32, default=0.01*cp.pi, 
                        help='Trapping frequency in ms^-1')
    parser.add_argument('--center_trap', type=cp.float32, nargs=2, default=(0, 0), 
                        help='Center of the trap in μm (e.g. ''(0, 0)'' means the center is at (0, 0) μm)')
    parser.add_argument('--beta', type=cp.float32, default=-0.1, 
                        help='Anharmonic parameter (dimensionless) (V = 0.5*m*omega^2*((1-2*beta)*(x^2+y^2) + (beta/r_0)*(x^2+y^2)^2)), where r_0 is the length scale to keep beta dimensionless)')
    parser.add_argument('--omega_trap_z', type=cp.float32, default=0.005*cp.pi, 
                        help='Trapping frequency along z in ms^-1')
    parser.add_argument('--r_0', type=cp.float32, default=30, 
                        help='Length scale to keep beta dimensionless (e.g. ''1'' means the length scale is 1 μm)')
    parser.add_argument('--center_bec', type=cp.float32, nargs=2, default=(30, 0), 
                        help='Center of the BEC in μm (e.g. ''25, 25'' means the center is at (25, 25) μm, and r_0 will be set to sqrt(2)*25 μm)')
    parser.add_argument('--omega_bec', type=cp.float32, default=0.01*cp.pi, 
                        help='The frequency of the harmonic trap whose ground state is this gaussian wavepacket')
    parser.add_argument('--angular_momentum_bec', type=cp.int32, nargs=2, default=(0, 0),
                        help='The angular momentum eigenstate (l, m) (e.g. ''(2, 1)'' means l = 2, m = 1)')
    parser.add_argument('--imaginary_time', action='store_true', 
                        help='If set, perform imaginary time evolution')
    parser.add_argument('--velocity', type=cp.float32, nargs=2, default=(0, 3*cp.pi/10),
                        help='Velocity of the initial wavepacket in μm/ms (e.g. ''(0, 0)'' means no initial velocity)')
    parser.add_argument('--video', action='store_true', 
                        help='If set, save the simulation process as a video')
    parser.add_argument('--figure', action='store_true', 
                        help='If set, save the simulation results as a figure')
    parser.add_argument('--mechanics', action='store_true', 
                        help='If set, save the mechanical quantities during the simulation')
    return parser

def display_args(parser:ArgumentParser) -> None:
    '''
    functionality:
        display the parsed arguments
    '''
    print("parameters:")
    print("\n## Atomic parameters ##")
    print(f"relative atomic mass: {parser.relative_atomic_mass} amu")
    print(f"scattering length: {parser.scattering_length} a_Bohr")
    print(f"number of atoms: {parser.atom_number}")
    print("\n## time structure ##")
    print(f"duration: {parser.duration} ms")
    print(f"time step: {parser.dt} ms")
    print(f"sampling interval: {parser.sampling_interval}")
    print("\n## grid structure ##")
    print(f"radius in x: {parser.radius_xy[0]} μm, radius in y: {parser.radius_xy[1]} μm")
    print(f"number in x: {parser.number_xy[0]}, number in y: {parser.number_xy[1]}")
    print("\n## potential trap ##")
    print(f"trapping frequency ω: {parser.omega_trap} ms^-1")
    print(f"center of the trap: ({parser.center_trap[0]}, {parser.center_trap[1]}) μm")
    print(f"anharmonic parameter β: {parser.beta}")
    print(f"length scale r_0: {parser.r_0} μm")
    print(f"trapping frequency along z direction ω_z: {parser.omega_trap_z} ms^-1")
    print("\n## initial wavepacket ##")
    print(f"center of the BEC: ({parser.center_bec[0]}, {parser.center_bec[1]}) μm")
    print(f"trapping frequency of BEC ω_BEC: {parser.omega_bec} ms^-1")
    print(f"initial velocity: ({parser.velocity[0]}, {parser.velocity[1]}) μm/ms")
    print(f"initial angular momentum eigenstate (l, m): ({parser.angular_momentum_bec[0]}, {parser.angular_momentum_bec[1]})")
    print("\n## overall simulation settings ##")
    print(f"imaginary time evolution: {'on' if parser.imaginary_time else 'off'}")
    print(f"video record: {'on' if parser.video else 'off'}")
    print(f"mechanical quantities: {'on' if parser.mechanics else 'off'}")
    print(f"figure output: {'on' if parser.figure else 'off'}\n")
    return None