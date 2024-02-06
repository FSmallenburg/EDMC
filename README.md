# Event-driven Monte Carlo for the WCA potential


This is an event-driven Monte Carlo (EDMD) code for particles interacting via the Weeks-Chandler-Andersen (WCA) potential. It is based on this [hard-sphere EDMD](https://github.com/FSmallenburg/EDMD) code (in particular the "Single" variant), adapted to implement the rejection-free Monte Carlo method introduced by [Peters and de With](https://doi.org/10.1103/PhysRevE.85.026703).

Contains an implementation of Xoshiro256+, see [Ref.](https://doi.org/10.1145/3460772), and [this page](https://prng.di.unimi.it/).

Code published along with the paper: (LINK TO BE ADDED).




## Compilation details

A simple makefile is included, as well as a sample initial configuration. All simulation parameters are defined as global variables near the top of the main code file. See the comments there for details.



## Simulation units and output

We simulate systems of $N$ particles in a constant volume $V$ in three dimensions. Each particle $i$ has a position $\mathbf{r}_i$, a diameter $\sigma_i$, a mass $m_i$ (currently set equal for all particles), and a radius $R_i$. 

The simulation is currently set up to handle either monodisperse or bidisperse systems. Other size distributions should be possible with minor modifications.

The simulation code operates in the following units:
-  Lengths are measured in units of the maximum particle size $\sigma$.
-  Mass is measured in units of a reference mass $m$ (currently just the mass of any particle, since they are all equal).
-  Time is measured in units of $\tau = \sqrt{\beta m \sigma^2}$. Here, $\beta = 1/k_B T$ with $k_B$ Boltzmann's constant.

Note that the simulation assumes that no particles with a diameter greater than $1 \sigma$ exist in the simulation, and uses this assumption in the creation of the cell list. Hence, all particles necessarily have a diameter $\sigma_i \leq \sigma$. By default, the simulation code sets the mass of all particles to be equal to $m$ when loading an initial configuration, but mass is explicitly taken into account when determining the effect of collisions. Hence, other choices for the mass can readily be implemented by adapting the initialization functions.

The simulation code measures the pressure $P$ during the simulation, and outputs it in the form of a reduced pressure $P^* = \beta P \sigma^3$. The average pressure in a given time interval $[t_a, t_b]$ is measured via the virial expression 
$P = \rho k_B T + \frac{1}{3V(t_b-t_a)} \sum  \mathbf{\delta p}_i,$ 
where  $\rho = N/V$ is the number density with $V$ the system volume and $N$ the number of particles, and the sum in the last term is taken over all collisions in the time interval $[t_a, t_b]$. For each collision, $\mathbf{r}_{ij}$ is the center-to-center vector connecting the two colliding particles $i$ and $j$, and $\mathbf{\delta p}_i$ is the momentum change of particle $i$ due to the collision. 

Similarly, the code outputs the potential energy $U/N\epsilon$. Both the pressure and potential energy are written to an output file (filename starting with "press").

    
## Snapshot file format

The configuration files from the simulation are written in a simple text-based format, which can contain multiple snapshots per file. For each frame, the format consists of $N+2$ lines (with $N$ the number of particles), as follows:
- One line containing just the number of particles
- One line containing the box size, specifying the box length $L_x$, $L_y$, and $L_z$ along the three axes, separated by whitespace.
- One line per particle containing: a letter indicating particle type, three numbers indicating the real-space particle coordinates, and one number indicating the particle radius. 

The movie files that are created by the simulation code include multiple of these frames consecutively in a single text file. Note that although the code assumes periodic boundary conditions, coordinates of particles that leave the box during the simulation will be printed as being outside of the simulation box, to allow for analysis of long-time dynamics. Hence, any structural analysis or visualization should apply periodic boundary conditions explicitly. 

The simulation codes can read in snapshots in this format as initial configurations (depending on the choice of "initialconfig" in the code). Periodic boundaries will be applied to the snapshot at the start of the simulation.  Note again that the simulation code assumes that all box lengths, positions, and radii are given in units of $\sigma$, which is the largest possible particle diameter. Hence, the radius of a particle in the initial configuration should never be given as a number larger than 0.5.

Adaptation to different configuration file formats can be done via modification of the ``loadparticles``, ``write``, and ``outputsnapshot`` functions.

