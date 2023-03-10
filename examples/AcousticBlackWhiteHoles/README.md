In the following examples, an acoustic black hole with a sonic horizon radius of rh = 1.5 m (directory <ABH_rh1.5>) and an acoustic white hole with the same radius (directory <AWH_rh1.5>) are simulated. Upon successful compilation of the source code in the build directory </path/to/some/arbitrary/build/directory/>, the simulation can be run by executing the following command in the corresponding case directory:

/path/to/some/arbitrary/build/directory/WaveDNA

In each of the case directories, the following post-processing scripts are available:

<spaceTimePlot.py>: This script generates a space-time plot for the acoustic pressure p1.
<plotWaveProfile.py>: This script generates a plot of one or several instantaneous wave profiles p1(r).
<plotHorizon.py>: This script generates a plot of the temporal acoustic pressure evolution p1(t) at the position rh of the sonic horizon.

See the documentation for more details.
