yeast7scripts
========================

A collection of scripts for handling the [consensus genome scale model of yeast (Yeast7.6)](https://github.com/SysBioChalmers/yeast-metabolic-network-7.6).
* `convertYmn2FBC2.m`: Converts yeast metabolic network 7.6 from COBRA-compatible SBML2 to FBC v2, thereby adding the missing annotation data, which could not be retained with the older COBRA versions.
* `modelCorrections.m`: Corrects various issues in yeast7 (biomass composition, proton balance, Ox.Pho., GAM and NGAM).
* `changeBiomass.m`: Rescales the biomass composition for varying protein content in anaerobic case. Also changes GAM and NGAM.
* `anaerobicModel.m`: Transforms the model to anaerobic conditions.
* `saveYeastModel.m`: Saves yeast model as an `sbml` and `txt` file, and creates `boundaryMets.txt`.
* `boundaryMets.txt`: Contains a list of all boundary metabolites in model, listing the id and name.
* `missingFields`: Folder with functions for adding missing fields to the model.
   * `addGeneNames.m`: Adds the field `geneNames` by extracting the data from KEGG. 
   * `getMissingFields.m`: Retrieves missing information (`rxnECNumbers` and `subSystems`) from KEGG & Swissprot. It uses `changeRules.m` for properly reading the gene-reaction rules, and `findInDB.m`, `getAllPath.m` and `findSubSystem.m` for reading the databases. The latter 3 functions are adapted versions of functions from the [GECKO toolbox](https://github.com/SysBioChalmers/gecko).
   * `ProtDatabase.mat`: Contains the relevant data from Swissprot and KEGG.

__Contributors:__ Benjamín J. Sánchez ([@BenjaSanchez](https://github.com/benjasanchez)) & Simonas Marcišauskas ([@simas232](https://github.com/simas232)).

Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology.

Last update: 2017-10-27

========================

## Required Software:

* MATLAB (7.5 or higher)
* The [COBRA toolbox for MATLAB](https://github.com/opencobra/cobratoolbox). Note that specific versions of COBRA must be used to properly convert to FBC v2 without losing functionality. For more details please read _convertYmn2FBC2.m_

========================
