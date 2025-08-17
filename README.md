# BEC_Dynamics_2D
Simulating and visualizing the mechanical motion of 2D Bose-Einstein condensates (BECs)

(This repository is under construction...)
### Physical Settings
The 2D BEC experiences a tunable anharmonic potential trap, and the inter-atomic interaction is treated by mean-field approximation.

### Numerical Method
Solving the time dependent Gross-Pitaevskii equation by split step method using FFT.

### Why 2D?
a) 2D simulation is much more faster than its 3D counterpart, and much more interesting than 1D case.  
b) Our eyes can only handle 2D images.

### Highlights
a) Simulation with units.  
b) Sampling the mechanical quantity (e.g. the angular momentum of the wavepacket) during time evolution.
c) GPU acceleration by CUPY.