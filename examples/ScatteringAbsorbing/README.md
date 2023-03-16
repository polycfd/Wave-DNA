This is a quick example demonstrating that the domain boundary can either be scattering (waves are reflecting) or absorbing (waves can pass the domain). Here, the set-up is as follows:

            West                                       East
        Moving o----o----o----o----o----o----o----o----o Fixed
               |                                       |
               |                                       |
            node 0                                    node N-1

The velocity of the moving boundary is zero so that the domain is stationary. A Gaussian pulse is emitted at the West boundary and traveling to the right. The boundary condition at the East is specified by setting <BoundaryConditionEast> to either <absorbing> or <scattering>. The condition for the West boundary ( <BoundaryConditionWest>) is immaterial in this case, as it is overwritten by the excitation function.

Upon successful compilation of the source code in the build directory </path/to/some/arbitrary/build/directory/>, the simulation can be run by executing the following command in the corresponding case directory:

/path/to/some/arbitrary/build/directory/WaveDNA

Using the above command, the options file is expected to have the name <run.DNA> (default). If the options file has a different name <optionsFileName>, the name must be specified as a command line option as follows:

/path/to/some/arbitrary/build/directory/WaveDNA -options optionsFileName

In each of the case directories, the following post-processing scripts are available:

<plotWaveProfile.py>: This script generates a plot of two instantaneous wave profiles p1(r), one before and one after the emitted waves has reached the domain boundary.

See the documentation for more details.
