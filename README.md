# LUNG PREDICT : Stade 4 NSCLC (Non-Small Cell Lung Cancer)

Dendritic cells spatial organization shapes tumor microenvironment and impact immunotherapy response in Non-Small-Cell lung cancer

Mains objective: we investigated the spatial organization of DC subsets in NSCLC and how this organization is able to shape the surrounding tumor microenvironment and transcription factor signaling, driving the immunotherapy response.

<p align="center">
<img src="framework.png?raw=true"/>
</p>
<p align="center">
<em>Figure 1. A schematic overview of the paper analysis</em>
</p>

Scripts
- Spatial_densities: Preprocessing of spatial data 
- Cell_densities_analysis: Analysis of DC subsets densities across response information and survival
- Megaclusters_analysis: Analysis of megaclusters density and CD8 enrichment across responders and patient stratification 
- RNAseq_analysis: Bulk RNAseq analysis of patient groups to characterize them with transcriptomic features (pathways, TFs modules, immunescores, chemokines, deconvolution)
- Signature_projection: Validation of TF signature in independent cohort
