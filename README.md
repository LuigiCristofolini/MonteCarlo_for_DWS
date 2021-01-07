# MonteCarlo_for_DWS


This is a template script, and you are free to adapt it your needs!

More details on the script, and on its usage, can be found in:

[0] V. Lorusso et al. Advances in Colloid and Interface Science 288 (2021) 102341
https://www.sciencedirect.com/science/article/pii/S0001868620306102

To know more about the physics behind DWS:

[1] Weitz, Pine, "Diffusing wave spectroscopy", in "Dynamic Light Scattering", Oxford University Press (1993)

[2] Durian Physical Review E 51, 3350 (1995)

[3] Pine, Weitz et al. Phys Rev Lett 60 1134 (1988)

MATLAB computing environment is chosen because  it is quite efficient by allowing implicit parallelization if the simulation is articulated in photon bunches.
Optimal bunch size is of the order of nParallel=10^5 photons with a typical Win10 64-bit running on 8 GB RAM.
 This number is optimized on my PC (Matlab 202ob on a W10 machine with 8Gb RAM), but please experience other values on you machine

Many bunches are successively and independently generated to build up a statistical ensemble.

The philosophy of the simulation is to follow photon propagation by a random walk process, keeping trace of the fate of each photon.
We choose a reference frame in which X is the horizontal direction of the impinging laser beam, Y is the other horizontal direction and Z is vertical.
At the initialization step, we define the cuvette size and two detection areas, one for backscattering and one for transmission experiments.
Initial positions for each photon are generated inside the sample, at a depth l^*, with lateral distribution in the YZ plane either uniform or gaussian to account for these two possible beam profiles.
Then the core of the Monte Carlo algorithm consists in generating the directions of the steps uniformly distributed over the sphere, while the step lengths,
following the approach of Durian et al. [2], follows a lognormal distribution with average value l^*, rather than being identically fixed to l^*.
This choice has no effect in the results for transmission geometry but yields better reproduction of the analytically exact formulas for the correlation functions simulated for the infinite slab in backscattering geometry.
At each step, the fate of each photon is checked: if it is still inside the sample volume, the counter of the path length is incremented by one unit and the photon is retained for the next steps of the simulation.
If on the contrary the photon has come out of the sample volume, simulation for this photon stops, and depending on whether this happened outside the measurement areas or inside one of them, the corresponding counter of the steps is discarded or retained in the appropriate (backscattering or transmission) pool.
This corresponds to the absorbing confinement condition, which we choose as the most realistic reproduction of the real experiment.
When enough photons have been simulated, from the pool of path lengths recorded for the backscattering and the transmission geometries, the corresponding path length distributions P(s) are generated and saved to a file.

The experimental geometry is described, in S.I. units, by
X is the direction of the laser (horizontal) impinging  beam
Y is the other horizontal
Z is vertical
(x0, y0, z0) initial positions
(xx, yy, zz) at each step, the instantaneous positions of the photons in their random walk

The shape of impinging laser beam can be chosen to be gaussian or homogeneous (if sigma_beam=Inf)
