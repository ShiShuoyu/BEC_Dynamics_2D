from argparse import ArgumentParser
import numpy as np

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
    parser.add_argument('--relative_atomic_mass', type=np.int32, default=87, 
                        help='Atomic mass of the particle in amu (e.g. ''87'' means Rb-87)')
    parser.add_argument('--scattering_length', type=np.float32, default=50,
                        help='Scattering length in Bohr radius (e.g. ''50'' means a_s = 50 a_Bohr)')
    parser.add_argument('--atom_number', type=np.int32, default=2500, 
                        help='Number of atoms')
    parser.add_argument('--duration', type=np.float32, default=200, 
                        help='Totle duration of simulation in ms')
    parser.add_argument('--dt', type=np.float32, default=1, 
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
    parser.add_argument('--omega_trap_z', type=np.float32, default=0.005*np.pi, 
                        help='Trapping frequency along z in ms^-1')
    parser.add_argument('--r_0', type=np.float32, default=30, 
                        help='Length scale to keep beta dimensionless (e.g. ''1'' means the length scale is 1 μm)')
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
    return parser

def display_args(parser:ArgumentParser) -> None:
    '''
    functionality:
        display the parsed arguments
    '''
    args = parser.parse_args()
    print("parameters:")
    print("\n## Atomic parameters ##")
    print(f"relative atomic mass: {args.relative_atomic_mass} amu")
    print(f"scattering length: {args.scattering_length} a_Bohr")
    print(f"number of atoms: {args.atom_number}")
    print("\n## time structure ##")
    print(f"duration: {args.duration} ms")
    print(f"time step: {args.dt} ms")
    print(f"sampling interval: {args.sampling_interval}")
    print("\n## grid structure ##")
    print(f"radius in x: {args.radius_xy[0]} μm, radius in y: {args.radius_xy[1]} μm")
    print(f"number in x: {args.number_xy[0]}, number in y: {args.number_xy[1]}")
    print("\n## potential trap ##")
    print(f"trapping frequency ω: {args.omega_trap} ms^-1")
    print(f"center of the trap: ({args.center_trap[0]}, {args.center_trap[1]}) μm")
    print(f"anharmonic parameter β: {args.beta}")
    print(f"length scale r_0: {args.r_0} μm")
    print(f"trapping frequency along z direction ω_z: {args.omega_trap_z} ms^-1")
    print("\n## initial wavepacket ##")
    print(f"center of the BEC: ({args.center_bec[0]}, {args.center_bec[1]}) μm")
    print(f"trapping frequency of BEC ω_BEC: {args.omega_bec} ms^-1")
    print(f"initial velocity: ({args.velocity[0]}, {args.velocity[1]}) μm/ms")
    print("\n## overall simulation settings ##")
    print(f"imaginary time evolution: {'on' if args.imaginary_time else 'off'}")
    print(f"video record: {'on' if args.video else 'off'}")
    return
    ...