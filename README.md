# BEC_Dynamics_2D
Simulating and visualizing the mechanical motion of 2D Bose-Einstein condensates (BECs)

(This repository is under construction...)
### Physical Settings
The 2D BEC experiences an anharmonic potential trap
$$
V(\vec{r}) = \frac{1}{2}m\omega^2\left((1-2\beta)\vec{r}^2+\frac{\beta}{r_0}\vec{r}^4\right).
$$  
And the inter-atomic interaction is treated by mean-field approximation, so the Hamiltonian reads
$$
H = \frac{\vec{p}^2}{2m} + \frac{1}{2}m\omega^2\left((1-2\beta)\vec{r}^2+\frac{\beta}{r_0}\vec{r}^4\right) + g|\psi|^2,
$$  
where $g$ is the interaction parameter.