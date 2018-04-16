# History

### yeast 8.0.1:
* `.yml` format included for easier visualization of model changes
* Empty notes removed from model
* Issue and PR templates included
* `README.md` updated to comply with new repo's name

### yeast 8.0.0:
First version of the yeast8 model, to separate it from previous versions:

* Manual curation project:
  * All metabolite information manually curated (names, charges, kegg IDs, chebi IDs)
  * Reaction gene rules updated with curation from [the iSce926 model](http://www.maranasgroup.com/submission_models/iSce926.htm). 13 genes added in this process
* Format changes:
  * Folder `ComplementaryData` introduced
  * All data is stored in `.tsv` format now (can be navigated in Github)
  * Releases now come in `.xlsx` as well
* Other new features:
  * Added `loadYeastModel.m`
  * A much smarter `increaseVersion.m`
  * Lots of refactoring

### yeast 7.8.3:
* curated tRNA's formulas
* started tracking COBRA and RAVEN versions
* dropped SBML toolbox as requirement
* reorganized `complementaryScripts`
* switched to a CC-BY-4.0 license

### yeast 7.8.2:
* fixed subSystems bug: now they are saved as individual groups
* solved inter-OS issues
* remade license to follow GitHub format
* added `history.md` and made it a requirement to update when increasing version

### yeast 7.8.1:
* started following dependencies
* started keeping track of the version in the repo (`version.txt`)
* included `.gitignore`
* dropped `.mat` storage for `devel` + feature branches (but kept it in `master`)

### yeast 7.8.0:
First release of the yeast model

* Format changes:
  * FBCv2 compliant
  * Compatible with latest COBRA and RAVEN parsers
  * Biomass clustered by 5 main groups: protein, carbohydrate, lipid, RNA and DNA
* Added information:
  * `subSystems` and `rxnECnumbers` added to reactions based on [KEGG](http://www.genome.jp/kegg/) & [Swissprot](http://www.uniprot.org/uniprot/?query=*&fil=organism%3A%22Saccharomyces+cerevisiae+%28strain+ATCC+204508+%2F+S288c%29+%28Baker%27s+yeast%29+%5B559292%5D%22+AND+reviewed%3Ayes) data
  * `geneNames` added to genes based on [KEGG](http://www.genome.jp/kegg/) data
  * `rxnKEGGID` added from old version
  * `rxnNotes` enriched with Pubmed ids (`pmid`) from old version
  * `rxnConfidenceScores` added based on [automatic script](https://github.com/SysBioChalmers/YeastMetabolicNetwork-GEM/blob/f7870589d16c08e18057a8f6cc880466373b77a7/ComplementaryScripts/getConfidenceScores.m)
  * Boundary metabolites tracked (available in [`ComplementaryScripts`](https://github.com/SysBioChalmers/yeast-metabolic-network-7.6/blob/master/ComplementaryScripts/boundaryMets.txt))
* Simulation improvements:
  * Glucan composition fixed in biomass pseudo-rxn
  * Proton balance in membrane restored
  * Ox.Pho. stoichiometry fixed
  * NGAM rxn introduced
  * GAM in biomass pseudo-rxn fixed and refitted to chemostat data