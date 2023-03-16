This example can be run with an empty options file, where the solver uses the default settings.

Upon successful compilation of the source code in the build directory </path/to/some/arbitrary/build/directory/>, the simulation can be run by executing the following command in the corresponding case directory:

/path/to/some/arbitrary/build/directory/WaveDNA

Using the above command, the options file is expected to have the name <run.DNA> (default). If the options file has a different name <optionsFileName>, the name must be specified as a command line option as follows:

/path/to/some/arbitrary/build/directory/WaveDNA -options optionsFileName

The following post-processing scripts are available:

<plotWaveProfile.py>: This script generates a plot of one or several instantaneous wave profiles p1(r).

See the documentation for more details.
