# Supervised_Optimal_Projection
Supplementary Materials include R code for paper under review

The computing scripts relied upon by this paper are twofold:
1) The scripts used to construct the unipartite shared-patient physician networks under the settings of the Continuity, Revisit and Multplicity factors as described in the Background section of the paper are written in Python. These are available at the repository: https://github.com/kiwijomalley/OptimalBipartiteProjection. For additional description about these factors, readers are referred to O'Malley et al (2022) at: https://doi.org/10.6339/22-JDS1064
2) The R scripts used to compute the various measures of diagnostic accuracy are contained in the file xyz.R. These receive as input a multi-edgelist containing the edges used as the gold standard (e.g., those obtained directly from surveying physicians) and the edge weights constructed using one or more projections defined by the three factors in 1) or otherwise. The data wrangling of the data needed to form the shared-patient weights is described in the GitHub page in 1) and so is not reproduced here.
