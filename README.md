---
title: "InflaMix"
author: "Sandeep Raj"
date: "2023-07-08"
---

## InflaMix Manuscript Code

This github repository contains code and models that accompany the InflaMix manuscript while it is under review. 

After publication, this git will be updated with:
 - a link to the paper
 - guidance on how to use InflaMix
 - additional code to further facilitate direct implementation of InflaMix for research with new data sets

A description of the two data inputs for this project and a description of their column values is given in the comments at the beginning of scaling_sfig3.R. 

Scripts should be run in the following order which contains code for the corresponding figures:

0 . Please place Dataset 1 (deriv_cohort_d0_labs.csv) into the data folder. 

1. scaling_sfig3.R 
  - Supplementary Figure 3
  
2. modelgen_efig2h_i.R 
  - Extended Figure 2h, i
  
3. icl_bic_efig2a.R 
  - Extended Figure 2a
  
4. efig2b_g.R 
  - Extended Figure 2b-g
  
5. deriv_cohort_properties_fig1a_g_sfig2.R 
  - Figure 1a-g, Supplementary Figure 2
  
6. deriv_outcomes_fig1h_l_table2.R
 - Figure 1h-l
 
7. partial_lab_clustering_fig2ab_fig2a_c.R
 - Figure 2a-b, Figure 3a-c
 
8. validation_fig2c_k_fig3d_l_efig4_table2.R
 - Figure 2c-k, Figure 3d-l, Extended Figure 4a-f
 
9. cohort_heatmaps_efig3a_d.R
 - Extended Figure 3a-d
 
10. exploratory_analyses_MTV_table2.R

11. cluster_transitions_figure6a_e.R
 - Figure 6a-e

Packages and version information for most dependencies are included in the renv.lock lockfile. All required libraries are listed at the top of each script. 

The InflaMix model features are in the "model" directory. 

## Data availability

Data requests for patient-related laboratory measurements or clinical outcomes will be reviewed by the corresponding author in consultation with coauthors from Hackensack University Medical Center and Sheba Medical Center. Any data and materials that can be shared will be released via data transfer agreement. Laboratory values and their corresponding upper limits of normal for the model-derivation cohort of MSK patients are provided in the Supplementary Data (Data Set 1). With this data, the following figures can be reproduced:
 - Figure 1 a-c, g
 - Extended Figure 2 a-i
 - Extended Figure 3 a
 - Supplementary Figure 2
 - Supplementary Figure 3
