# Wave-DNA

<p align="left">
  <a href="https://doi.org/10.5281/zenodo.7249297">
    <img src="https://img.shields.io/badge/DOI-10.5281/zenodo.7249297-blue" alt="Latest version">
  </a>
</p>
<p align="left">
  <a href="https://doi.org/10.1016/j.jsv.2022.116814">
    <img src="https://img.shields.io/badge/Finite--difference_method-green" alt="FDM">
  </a>
  <a href="https://doi.org/10.1121/10.0020561">
    <img src="https://img.shields.io/badge/Convective_Kuznetsov_equation-green" alt="JASA">
  </a>
  <a href="https://doi.org/10.1063/5.0106145">
    <img src="https://img.shields.io/badge/Acoustic_black_hole_analogue-green" alt="ABH">
  </a>
</p>

Wave-DNA is a software tool for one-dimensional and spherically-symmetric nonlinear acoustic waves in transient and spatially variable background flow fields. The motion of the background medium is accounted for by considering a convective form of the lossless Kuznetsov wave equation, derived from first principles based on perturbations of the continuity equation and the transient Bernoulli equation, which is solved using a finite-difference method. In principle, the background flow field may be obtained from analytical solutions, numerical simulations, or experimental measurements. The acronym Wave-DNA stands for "Wave Doppler effects in Nonlinear Acoustics".

## Developers
- [Sören Schenke](mailto:soeren.schenke@ovgu.de) (principal developer)
- [Fabian Sewerin](mailto:fabian.sewerin@ovgu.de)
- [Berend van Wachem](mailto:berend.vanwachem@ovgu.de)
- [Fabian Denner](mailto:fabian.denner@polymtl.ca) (maintainer)

## License and Copyright
Wave-DNA is under the copyright of its developers and made available as open-source software under the terms of the [MIT License](LICENSE).

## Repository Structure
The Wave-DNA repository is structured as follows:
- The [documentation](/documentation/) folder contains a short [pdf](/documentation/WaveDNA-Documentation.pdf) documentation of Wave-DNA. The documentation discusses the theory behind Wave-DNA and the code structure, and explains how to install and use Wave-DNA.
- The [examples](/examples/) folder contains representative examples that show how to use Wave-DNA and demonstrate the most important features of Wave-DNA. 
- The [src](/src/) folder contains all source files of Wave-DNA.
- The [.clang-format](.clang-format) file, which defines the formatting rules for the source code.
- The [.gitignore](.gitignore) file telling _git_ which folders and files to ignore.
- The [LICENSE](LICENSE) file containing the MIT License text.
- The [README.md](README.md) file is the file you are currently reading.

## Acknowledgements
The development of Wave-DNA has directly benefitted from research funding provided by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation), grant number 441063377.
