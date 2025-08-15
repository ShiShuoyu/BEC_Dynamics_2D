import numpy as np
import matplotlib.pyplot as plt
import cupy as cp
import tkinter as tk
import argparse

from constants import *
import wf
import evolve

def main():
    # Parse command line arguments # 命令行参数解析
    parser = argparse.ArgumentParser(
        description="Simulating and visualizing the machanical motion of 2D Bose-Einstein condensates (BECs)"
        )
    parser.add_argument('--radius_xy', type=np.float32, nargs=2, default=(50, 50), 
                        help='Radius of the real space simulation region in μm (e.g. ''(50, 40)'' means -50 to 50 μm in x and -40 to 40 μm in y)')
    parser.add_argument('--number_xy', type=int, nargs=2, default=(256, 256),
                        help='Number of points in x and y directions (e.g. ''(256, 256)'' means 256 points in both x and y directions) (When N is a power of 2, the FFT is faster)')
    parser.add_argument('--omega', type=np.float32, default=0.01*np.pi, 
                        help='Trapping frequency in ms^-1')
    parser.add_argument('--center_trap', type=np.float32, nargs=2, default=(0, 0), 
                        help='Center of the trap in μm (e.g. ''(0, 0)'' means the center is at (0, 0) μm)')
    parser.add_argument('--beta', type=np.float32, default=-0.1, 
                        help='Anharmonic parameter (dimensionless) (V = 0.5*m*omega^2*((1-2*beta)*(x^2+y^2) + (beta/r_0)*(x^2+y^2)^2)), where r_0 is the length scale to keep beta dimensionless)')
    parser.add_argument('--center_bec', type=np.float32, nargs=2, default=(30, 0), 
                        help='Center of the BEC in μm (e.g. ''25, 25'' means the center is at (25, 25) μm, and r_0 will be set to sqrt(2)*25 μm)')
    parser.add_argument('--velocity', type=np.float32, nargs=2, default=(0, 0),
                        help='Velocity of the initial wavepacket in μm/ms (e.g. ''(0, 0)'' means no initial velocity)')
    args = parser.parse_args()
    print("parameters:")
    print(f"radius in x: {args.radius_xy[0]} μm, radius in y: {args.radius_xy[1]} μm")
    print(f"number in x: {args.output[0]}, number in y: {args.output[1]}")
    print(f"trapping frequency ω: {args.omega} ms^-1")
    print(f"center of the trap: ({args.center_trap[0]}, {args.center_trap[1]}) μm")
    print(f"anharmonic parameter β: {args.beta}")
    print(f"center of the BEC: ({args.center_bec[0]}, {args.center_bec[1]}) μm")
    print(f"initial velocity: ({args.velocity[0]}, {args.velocity[1]}) μm/ms")


if __name__ == "__main__":
    main()