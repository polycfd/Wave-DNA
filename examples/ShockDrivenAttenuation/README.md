This example involves a progressively deforming, shock-forming and decaying wave in one-dimensional Cartesian coordinates, where the attenuation of the waveform is driven by the dissipation of energy in the shock front. Despite the absence of any physical attenuation term in the present release of Wave-DNA, the attenuation of higher wave frequencies can be achieved by the help of the predictor-corrector method.

The progressively deforming wave is obtained by setting <FluidNonLinearity 3.5>, and the predictor-corrector method is enabled by setting:

<nCorrectors 1>
<correctorWeight 1.5>

The <correctorWeight> determines the magnitude of the damping that acts predominantly on the higher wave frequencies, and <nCorrectors> is the number of corrector steps per time step, where the results hardly change after the first iteration.

Upon successful compilation of the source code in the build directory </path/to/some/arbitrary/build/directory/>, the simulation can be run by executing the following command in the corresponding case directory:

/path/to/some/arbitrary/build/directory/WaveDNA

Using the above command, the options file is expected to have the name <run.DNA> (default). If the options file has a different name <optionsFileName>, the name must be specified as a command line option as follows:

/path/to/some/arbitrary/build/directory/WaveDNA -options optionsFileName

In each of the case directories, the following post-processing scripts are available:

<plotWaveProfile.py>: This script generates a plot of an instantaneous wave profile p1(r). In addition to that, the envelope of the Fay solution based on the nominal shock formation is indicated, which can be found in:

David T. Blackstock. Connection between the Fay and Fubini solutions for plane sound waves of finite amplitude. The Journal of the Acoustical Society of America, 39(6):1019–1026, 1966. doi:10.1121/1.1909986.

See the documentation for more details.
