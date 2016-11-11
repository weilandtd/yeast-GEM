# Genome-Scale Model of _S. cerevisiae_ - Aerobic Branch

This repository contains the current genome-scale metabolic model of _Saccharomyces cerevisiae_ used in the [@SysBioChalmers](https://github.com/SysBioChalmers) group. It's an improved version of [the consensus metabolic model, version 7.6](https://sourceforge.net/projects/yeast/), FBCv2 compliant, with subSystem and EC number information connected to reactions, and with several issues corrected. The usage of this model is exactly the same as any other COBRA developed model.

This repository is administered by Benjamín J. Sánchez ([@BenjaSanchez](https://github.com/benjasanchez)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

Last update: 2016-11-11

## Installation

### Required Software:
* You need a functional Matlab installation (MATLAB 7.3 or higher)
* The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox). An up-to-date version from COBRA GitHub repository is strongly recommended . Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com).

### Dependencies - Recommended Software:
* libSBML MATLAB API ([version 5.13.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/)  is recommended).
* Gurobi Optimizer for MATLAB ([version 6.5.2](http://www.gurobi.com/registration/download-reg) is recommended). 

### Installation Instructions
* Clone aerobic branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers/yeast-metabolic-network-7.6)
* Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)

## Contributors
* [Benjamín J. Sánchez](https://www.chalmers.se/en/staff/Pages/bensan.aspx), Chalmers University of Technology, Gothenburg Sweden
