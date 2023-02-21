# Wave-DNA
?

Wave-DNA is a software toolbox/program (?) to ??. The acronym Wave-DNA stands for ""??"".

Key features of Wave-DNA are:
- ??
- ??

## Developers
- [Sören Schenke](mailto:soeren.schenke@ovgu.de) (principal developer)
- ??
- [Fabian Denner](mailto:fabian.denner@ovgu.de) (maintainer)

## License and Copyright
Wave-DNS is under the copyright of its developers and made available as open-source software under the terms of the [MIT License](LICENSE).

## Repository Structure
The Wave-DNA repository is structured as follows:
- The [documentation](/documentation/) folder contains a short [pdf](/documentation/WDNA-Documentation.pdf) documentation of Wave-DNA. The documentation discusses the theory behind Wave-DNA, explains the code structure and how to use Wave-DNA.
- The [examples](/examples/) folder contains representative examples of how to use Wave-DNA and to demonstrate the most important features of Wave-DNA. A short explanation on how to run the examples is given in the [Quick Start Guide](#quick-start-guide) below.
- The [src](/src/) folder contains all source files of Wave-DNA.
- The [.clang-format](.clang-format) file, which defines the formatting rules for the source code.
- The [.gitignore](.gitignore) file telling _git_ which folders and files to ignore.
- The [LICENSE](LICENSE) file containing the MIT License text.
- The [README.md](README.md) file is the file you are currently reading.

## Quick Start Guide
Getting started with Wave-DNA is easy. After downloading Wave-DNA in the directory ````<path to Wave-DNA>````, define the following environment variables:
- ````WDNA_DIR```` to the directory in which Wave-DNA is located. Using bash, for instance, simply execute the command ````export WDNA_DIR=<path to Wave-DNA>```` or, even better, add this command to your bash profile.
- ````USRLIB_DIR```` to the directory in which libm.a or libm.dylib (the standard _math_ library) is located. This may, for instance, be ````/usr/lib64/```` on Linux systems or ````/usr/lib/```` on MacOS systems.

## Acknowledgements
The development of Wave-DNA has directly benefitted from research funding provided by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation), grant number 441063377.
