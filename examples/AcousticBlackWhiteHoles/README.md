In the following examples, an acoustic black hole (directory _BlackHole_) and an acoustic white hole (directory _WhiteHole_) are simulated. By default, the radius of the sonic horizon radius is $r_h = 1.5 \ \mathrm{m}$ in both cases.

Upon successful compilation of the source code in the build directory ````/path/to/some/arbitrary/build/directory/````, the simulation can be run by executing the following command in the corresponding case directory:

````/path/to/some/arbitrary/build/directory/WaveDNA````

Using the above command, the options file is expected to have the name ````run.DNA```` and be located in the directory in which Wave-DNA is executed. If the options file has a different name ````optionsFileName````, the name must be specified as a command line option as follows:

````/path/to/some/arbitrary/build/directory/WaveDNA -options optionsFileName````

In each of the case directories, the following post-processing scripts are available:
- ````spaceTimePlot.py````: This script generates a space-time plot for the acoustic pressure $p_1$.
- ````plotWaveProfile.py````: This script generates a plot of one or several instantaneous wave profiles $p_1(r)$.
- ````plotHorizon.py````: This script generates a plot of the temporal acoustic pressure evolution $p_1(t)$ at the position $r_h$ of the sonic horizon.

See the documentation for more details.
