# Rlv3481-lifestyles
Scripts associated with the manuscript "Genome-scale metabolic modelling of lifestyle changes in Rhizobium leguminosarum" by C. C. M. Schulte, V. K. Ramachandran, A. Papachristodoulou and P. S. Poole

All MATLAB scripts were run with the COBRA Toolbox v.3.0 ([Heirendt *et al.*, 2019](https://www.nature.com/articles/s41596-018-0098-2)).


```
|- README
|
|- iCS1224.mat                   # MATLAB structure of iCS1224, a metabolic model for Rhizobium leguminosarum bv. viciae 3841 
|- iCS1224.xlsx                  # Excel file for iCS1224, a metabolic model for Rhizobium leguminosarum
|                                bv. viciae 3841 
|
|- Scripts/                     
|
| |- UMSsimulations.m            # MATLAB script to compare model flux predictions with 13C flux data
| |- bacteroidBody_riptide.m     # MATLAB script to prepare bacteroid model for use with RIPTiDe 
| |- noduleBacteria_riptide.m    # MATLAB script to prepare nodule bacteria model for use with RIPTiDe 
| |- rhizosphere_riptide.m       # MATLAB script to prepare rhizosphere bacteria for use with RIPTiDe 
| |- BiologComparison.m          # MATLAB script to compare Biolog data with model predictions
| |- geneEssentiality.m          # MATLAB script to compare gene essentiality data (INSeq) with model predictions
| |- geneEssentiality.m          # MATLAB script to determine reporter metabolites in the rhizosphere of pea, alfalfa, sugar beet
| |- RIPTiDe_bacteroidBody.py    # Python script to run RIPTiDe for bacteroids
| |- RIPTiDe_noduleBacteria.py   # Python script to run RIPTiDe for nodule bacteria 
| |- RIPTiDe_rhizosphere.py      # Python script to run RIPTiDe for rhizosphere bacteria 
|
|- Data/             
|
| |- bacteroidBody_rpkm          # dRNA-Seq data for nodule body for use with RIPTiDe
| |- bacteroidTips_rpkm          # dRNA-Seq data for nodule tips for use with RIPTiDe
| |- rhizosphere_WT_rpkm         # RNA-Seq data for rhizosphere bacteria for use with RIPTiDe
| |- bacteroidCpds               # Compounds available to bacteroids for use with RIPTiDe
| |- noduleTipCpds               # Compounds available to nodule bacteria for use with RIPTiDe
| |- rhizosphereCpds             # Compounds available to rhizosphere bacteria for use with RIPTiDe
| |- bacteroidGenes              # Genes associated with proteins upregulated in bacteroids for use with RIPTiDe
| |- noduleBacteriaGenes         # Genes essential for nodule bacteria for use with RIPTiDe
| |- rhizosphereGenes            # Genes essential for rhizosphere bacteria for use with RIPTiDe
| |- rhizosphereCpds             # Compounds available to rhizosphere bacteria for use with RIPTiDe
| |- BiologData                  # Results of phenotype microarray experiment
| |- C13data                     # 13C flux data for R. leguminosarum
| |- essentialGenes              # Genes essential for growth in minimal media (INSeq) to compare with model predictions
| |- MA_rhizosphere_pea          # Microarray data for pea rhizosphere vs free-living growth
| |- MA_rhizosphere_alfalfa      # Microarray data for alfalfa rhizosphere vs free-living growth
| |- MA_rhizosphere_sugarbeet    # Microarray data for sugar beet rhizosphere vs free-living growth
|
|- Results/             
|
| |- bacteroidBody/              # RIPTiDe output for bacteroids
| |- bacteroidBody_o2lim/        # RIPTiDe output for oxygen-limited bacteroids
| |- noduleBacteria/             # RIPTiDe output for nodule bacteria
| |- rhizosphere/                # RIPTiDe output for rhizosphere bacteria
| |- reporterMets_pea            # Reporter metabolites for pea rhizosphere
| |- reporterMets_alfalfa        # Reporter metabolites for alfalfa rhizosphere
| |- reporterMets_sugarbeet      # Reporter metabolites for sugar beet rhizosphere
|
```
