# YeastMetabolicNetwork-GEM - Aerobic Branch

* Brief Model Description:

This repository contains the current genome-scale metabolic model of _Saccharomyces cerevisiae_ used in the [@SysBioChalmers](https://github.com/SysBioChalmers) group. It is an improved version of [the consensus metabolic model, version 7.6](https://sourceforge.net/projects/yeast/).

* Main Improvements to Original Model:

  * Format changes:
    * FBCv2 compliant.
    * Compatible with latest COBRA and RAVEN parsers.
  * Added information:
    * `subSystems` and `rxnECnumbers` added to reactions based on [KEGG](http://www.genome.jp/kegg/) & [Swissprot](http://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Saccharomyces+cerevisiae+%28strain+ATCC+204508+%2F+S288c%29+%28Baker%27s+yeast%29+%5B559292%5D%22+AND+reviewed%3Ayes) data.
    * `geneNames` added to genes based on [KEGG](http://www.genome.jp/kegg/) data.
    * Boundary metabolites tracked (available in [`ComplementaryScripts`](https://github.com/SysBioChalmers/yeast-metabolic-network-7.6/blob/master/ComplementaryScripts/boundaryMets.txt)).
  * Simulation improvements:
    * Glucan composition fixed in biomass pseudo-rxn.
    * Proton balance in membrane restored.
    * Ox.Pho. stoichiometry fixed.
    * NGAM rxn introduced.
    * GAM in biomass pseudo-rxn fixed and refitted to chemostat data.

* Model KeyWords:

**GEM Category:** Species; **Utilisation:** maximising growth; **Field:** metabolic-network reconstruction; **Type of Model:** curated, reconstruction; **Model Source:** [Yeast 7.6](https://sourceforge.net/projects/yeast/); **Taxonomy:** _Saccharomyces cerevisiae_; **Metabolic System:** General Metabolism; **Bioreactor**

* Last update: 2017-10-30

* The model:

|Taxonomy | Template Model | Reactions | Metabolites| Genes |
|:-------:|:--------------:|:---------:|:----------:|:-----:|
|_Saccharomyces cerevisiae_|[Yeast 7.6](https://sourceforge.net/projects/yeast/)|3494|2220|909|


This repository is administered by Benjamín J. Sánchez ([@BenjaSanchez](https://github.com/benjasanchez)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.


## Installation

### Required Software:

* A functional Matlab installation (MATLAB 7.3 or higher)
* One of the 2 following:
  * The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox).
  * The [RAVEN toolbox for MATLAB](https://github.com/SysBioChalmers/RAVEN). An up-to-date version from COBRA GitHub repository is strongly recommended . Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com).

### Dependencies - Recommended Software:
* libSBML MATLAB API (version [5.15.0](https://sourceforge.net/projects/sbml/files/libsbml/5.15.0/stable/MATLAB%20interface/) is recommended).
* The SBML Toolbox (version [4.1.0](https://sourceforge.net/projects/sbml/files/SBMLToolbox/4.1.0/)  is recommended).
* Gurobi Optimizer for MATLAB (version [6.5.2](http://www.gurobi.com/registration/download-reg) is recommended). 

### Installation Instructions
* Clone aerobic branch from [SysBioChalmers GitHub](https://github.com/SysBioChalmers/yeast-metabolic-network-7.6).
* Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com).


## Contributors
- [Benjamín J. Sánchez](https://www.chalmers.se/en/staff/Pages/bensan.aspx) ([@BenjaSanchez](https://github.com/benjasanchez)), Chalmers University of Technology, Gothenburg Sweden
- [Simonas Marcišauskas](https://www.chalmers.se/en/Staff/Pages/simmarc.aspx) ([@simas232](https://github.com/simas232)), Chalmers University of Technology, Gothenburg Sweden
- [Hongzhong Lu](https://www.chalmers.se/en/Staff/Pages/luho.aspx) ([@hongzhonglu](https://github.com/hongzhonglu)), Chalmers University of Technology, Gothenburg Sweden