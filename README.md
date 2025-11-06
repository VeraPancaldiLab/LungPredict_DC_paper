# LUNG PREDICT : Stade 4 NSCLC (Non-Small Cell Lung Cancer)

Dendritic cells spatial organization shapes tumor microenvironment and impact immunotherapy response in Non-Small-Cell lung cancer

Mains objective: we investigated the spatial organization of DC subsets in NSCLC and how this organization is able to shape the surrounding tumor microenvironment and transcription factor signaling, driving the immunotherapy response.

<p align="center">
<img src="framework.png?raw=true"/>
</p>
<p align="center">
<em>Figure 1. A schematic overview of the paper analysis</em>
</p>

## Data used in the analysis
- 116 samples with bulk RNAseq
- 68 samples with spatial
- 60 samples with bulk RNAseq and spatial
- 58 samples analyzed 
  
## Project organization
**input**: Files used for the analysis
- `Deconvolution_Gobbini`: Cell type deconvolution matrix  
- `MegaClusters_area`: DC Megaclusters area information
- `MegaClusters_density`: DC subsets density information
- `MegaClusters_entropy`: DC Megaclusters CD8 enrichment information
- `Spatial_densities.RData`: DC subsets densities after preprocessing
- `signatures/`: DC signatures from scRNAseq used for cell type deconvolution
  
**output**: Output files
- `Patients_groups.RData`: Patients groups assignations 
- `group_marker_genesets.rds`: TFs sets fromm differential analysis across groups
  
**Scripts**: Codes used for analysis.
- `Spatial_densities`: Preprocessing of spatial data 
- `Cell_densities_analysis`: Analysis of DC subsets densities across response information and survival
- `Megaclusters_analysis`: Analysis of megaclusters density and CD8 enrichment across responders and patient stratification 
- `RNAseq_analysis`: Bulk RNAseq analysis of patient groups to characterize them with transcriptomic features (pathways, TFs modules, immunescores, chemokines, deconvolution)
- `Signature_projection`: Validation of TF signature in independent cohort
  
**Results**: Figures used in the paper.

## Specifications
Analysis was done using R version 4.3.1 with the OS Ubuntu 22.04.3 LTS.

## Reproducibility
If you would like to reproduce the analysis done here, we invite you to use our provided r-environment. Setting it up will install all the neccessary packages, along with their specific versions in an isolated environment.  

For this, open the project `LP_MOSAIC.Rproj` inside the `scripts/` folder and in the R console run:

```r
# Download renv package (if not installed)
install.packages('renv')
# To activate the R environment 
renv::activate()
# To download and install all the require libraries and packages 
renv::restore() 
```

Note that this is an **once-step** only when running the repository for the **first time**. For the following times, you will only need to open the `LP_MOSAIC.Rproj` and you are ready to go!

Once all packages have been installed, you can start reproducing the analysis using the scripts inside the `scripts/` folder.

Make sure to run `renv::deactivate()` when finishing, to avoid conflicts whenever you start a different R project.

For more information about how R-environments work, visit the main page of the tool [renv](https://rstudio.github.io/renv/articles/renv.html).

## Citation
Soon available

## Contributing
If you are interested or have questions about the analysis done in this project, we invite you to open an issue in https://github.com/VeraPancaldiLab/LungPredict_DC_paper/issues or contact Marcelo Hurtado (marcelo.hurtado@inserm.fr) for more information.

## Authors
- [Marcelo Hurtado](https://github.com/mhurtado13)
- [Vera Pancaldi](https://github.com/VeraPancaldi)
- Elisa Gobbini
- Jenny Valladeau
- Pierre Duplouye
