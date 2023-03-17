The following examples involve oscillating wave-emitting boundaries in spherical symmetry (the cases can just as well be run in one-dimensional Cartesian coordinates). In the sub-directory _FixedBoundary_3d_R0.015_, the wave-emitting spherical boundary is, for reference, fixed and has a radius of $R_0 = 0.015 \ \mathrm{m}$ (one wavelength) in order to reproduce the classical $1/r$-decay of a spherical wave. In the sub-directory _MovingBoundary_3d_R0.015_inducedFlow_, the oscillating wave-emitting boundary displaces the surrounding fluid and thus induces a spherically symmetric transient and spatially non-uniform flow field. In the sub-directory _MovingBoundary_3d_R0.015_quiescent_, the wave-emitting boundary moves in precisely the same fashion as in the previous case, however in a quiescent background medium (the boundary is in relative motion to the medium).

It is important to note that the specified excitation pressure amplitude is enforced as a boundary condition. Hence, the Mach number dependent amplitude shift in the classical case of convective amplification/attenuation would have to be considered in the boundary condition as given by the excitation pressure amplitude in the options file.

In spherical symmetry, the velocity amplitude of the sinusoidal boundary motion must be chosen in such a away that the coordinate of the boundary does not become zero or negative at any time instant (singularity of the Laplacian).

Upon successful compilation of the source code in the build directory ````</path/to/some/arbitrary/build/directory/>````, the simulation can be run by executing the following command in the corresponding case directory:

````/path/to/some/arbitrary/build/directory/WaveDNA````

Using the above command, the options file is expected to have the name ````run.DNA```` and be located in the directory in which Wave-DNA is executed. If the options file has a different name ````optionsFileName````, the name must be specified as a command line option as follows:

````/path/to/some/arbitrary/build/directory/WaveDNA -options optionsFileName````

In each of the case directories, the following post-processing scripts are available:
- ````spaceTimePlot.py````: This script generates a space-time plot for the acoustic pressure $p_1$, indicating the boundary motion.
- ````plotWaveProfile.py````: This script generates a plot of one or several instantaneous wave profiles $p_1(r)$, also indicating the slope of $1/r$-decay for a spherical wave.

See the documentation for more details.
