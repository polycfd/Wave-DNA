In the following examples, the Doppler effect for a constant motion of the velocity field and/or the emitter can be reproduced in one-dimensional Cartesian coordinates. The case directory _InducedFlow_ involves an example in which the wave emitting boundary moves and induces a uniform flow field by displacing the fluid. This case can be seen as the classical Doppler shift, where the fixed Eulerian receiver is in relative motion to the sourc.

The case directory _RelativeEmitterMotion_ involves an example in which the wave emitting boundary is again moving, but this time relative to a quiescent background medium. It is important to note that the specified excitation pressure amplitude is enforced as a boundary condition. Hence, the Mach number dependent amplitude shift in the classical case of convective amplification/attenuation would have to be considered in the boundary condition as given by the excitation pressure amplitude in the options file.

The case directory _PossibleCombinations_ includes eight sub-directories illustrating the possible combinations of moving and wave-emitting boundaries. The physical domain is sketched below, indicating the West and East boundaries. In the options file, the initial boundary positions are specified under ````FixedBoundaryPosition```` and ````InitialMovingBoundaryPosition````. The one with the smaller coordinate is identified as the West boundary.


        West o----o----o----o----o----o----o----o----o East


Any of the two bondaries may be the moving one. Independent of that, any of the two boundaries may be the wave-emitting one, specified by the excitation node. The moving boundary can either induce a flow field (````BackgroundMotionMode coupledToMovingBoundary````) or move relative to a quiescent background medium (````BackgroundMotionMode quiescent````). This gives the eight combinations indicated by the names of the sub-directories. The boundary conditions are assigned to the boundaries identified as West and East, specified under ````BoundaryConditionWest```` and ````BoundaryConditionEast````. For the wave-emitting boundary, the boundary condition is immaterial as it is overwritten by the excitation function. The other boundary can either be ````absorbing```` to let the wave pass the domain, or ````scattering```` to reflect the wave.

Furthermore, one could further run these configurations in spherical symmetry and/or with an oscillating emitter, so that the flow field becomes non-uniform. However, the coordinate range of the physical domain must be larger than zero then as the equation becomes singular and/or ill-posed for negative coordinates.

Important notes regarding the emission node:

1.) The node counting starts from 0. This means that if the domain consists of N points, spcified under ````NPoints````, the boundary opposite to the one associated with node 0 is associated with node N-1.

2.) The numerical algorithm is implemented in such a way that the moving boundary boundary is always associated with node 0, regardless of whether it is the West or East boundary. This means that if the wave shall be excited at the moving boundary, <ExcitationNode> must be set to 0. If the wave shall be excited at the boundary opposite to the moving one,<ExcitationNode> must be set to N-1. This is illustrated below.

3.) In principle, <ExcitationNode> can be set to any int value in the range [0, N-1], but it is emphasized that the simulation tool is designed to solve boundary value problems.

            West                                       East
        Moving o----o----o----o----o----o----o----o----o Fixed
               |                                       |
               |                                       |
            node 0                                    node N-1


            West                                       East
         Fixed o----o----o----o----o----o----o----o----o Moving
               |                                       |
               |                                       |
            node N-1                                  node 0

Upon successful compilation of the source code in the build directory ````/path/to/some/arbitrary/build/directory/````, the simulation can be run by executing the following command in the corresponding case directory:

````/path/to/some/arbitrary/build/directory/WaveDNA````

Using the above command, the options file is expected to have the name ````run.DNA```` and be located in the directory in which Wave-DNA is executed. If the options file has a different name ````optionsFileName````, the name must be specified as a command line option as follows:

````/path/to/some/arbitrary/build/directory/WaveDNA -options optionsFileName````

In each of the case directories, the following post-processing scripts are available:
- ````plotWaveProfile.py````: This script generates a plot of one or several instantaneous wave profiles $p_1(r)$.

See the documentation for more details.
