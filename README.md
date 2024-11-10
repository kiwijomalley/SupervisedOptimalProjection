# SupervisedOptimalProjection
Supplementary Materials for the paper "Methodology for Supervised Optimization of the Construction of Physician Shared-Patient Networks" that is currently under review. This includes R code used in the analyses described in the paper and the Appendix to be published with the paper,

## R Scripts
The computing scripts relied upon by this paper are twofold:
1) The scripts used to construct the unipartite shared-patient physician networks under the settings of the Continuity, Revisit and Multplicity factors as described in the Background section of the paper are written in Python. These are available at the repository: https://github.com/kiwijomalley/OptimalBipartiteProjection. For additional description about these factors, readers are referred to O'Malley et al (2022) at: https://doi.org/10.6339/22-JDS1064
2) The R scripts used to compute the various measures of diagnostic accuracy are contained in the file xyz.R. These receive as input a multi-edgelist containing the edges used as the gold standard (e.g., those obtained directly from surveying physicians) and the edge weights constructed using one or more projections defined by the three factors in 1) or otherwise. The data wrangling of the data needed to form the shared-patient weights is described in the GitHub page in 1) and so is not reproduced here.

The 5 R scripts and their functions follow:
1) ConcordDataBuildNew.r
   - Step 1 of data wrangling process to perform concordance analysis of NPO survey nominations and patient-sharing in Medicare across the 20 projections.
   - Saves the generated networks for a summarization program (ConcordStatsAllNew.r) to then perform the concordance analysis.

2) ConcordDataWrangle.r
   - Step 2 of data wrangling process
   - Standardizes form of relational data sets across the different ways of vonstructing the network
   - Assumes that all NPO physicians at a hospital know each other

3) ConcordStatsAllNew.r
   - Step 1 of data analyses
   - Performs concordance analysis of NPO survey nominations and patient-sharing in Medicare
   - Makes plots of networks formed using different projections
   - Runs on the 40 different ways of forming weighted physician networks (including the 20 undirected projections)
  
4) ConcordStatsComposite.r
   - Step 2 of data analyses
   - Estimates hierarchical p2 and p2 models for relational data to try and establish whether the direction of patient-sharing that has the greater weight is the more likely to correspond to a true nomination.
   - Also evaluates descriptive statistics

5) p2est.r
   - #R function for setting up and estimating p2 and hierarchical p2 models, respectively.
   - Use type=1 or 2 for the hierarchical p2 model
   - Use type=3 or 4 for the p2 model
   - The difference between type=1 vs type=2 and type=3 vs type=4 is immaterial once inside this function. It reflects whether a single model is being estimated or multiple are (i.e., this procedure is called multiple times by an outside function)

## Appendix
The Appendix accompanying the paper is contained in the file: Supervised_Optimal_Projection_Revision_AppendixUnblinded.pdf. 

Section A of the Appendix includes additional dyadic diagnostic measures and conditional dyadic analyses. Section B contains additional empirical results relevant to the application: the demonstration of the minimal risk of overfitting from not forming separate training and test datasets when estimating and comparing the performance of different network construction methods, additional plots of the networks, associations of the Multiplicity factor levels with edge-level diagnostic accuracy, and the fitted Hierarchical P2 model used to evaluate the overall association of the within-hospital shared-patient physician networks with the survey nomination physician networks.
