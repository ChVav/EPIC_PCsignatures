
This version of the code accompanies a manuscript submitted to bioRxiv: <br>
Contradictory Results: Vavourakis, C. D., Herzog, C. M., & Widschwendter, M. (2023). <br>
Devising reliable and accurate epigenetic predictors: choosing the optimal computational solution (submitted). <br>
https://doi.org/10.1101/2023.10.13.562187

## Scope

* We evaluated the precision and accuracy of a PC-version of a clock predicting chronological age from DNA methylation array data, <br>
an approach described by Higgins-Chen et al (https://doi.org/10.1038/s43587-022-00248-2) for bolstering reliability of epigenetic clocks. <br>
* We retrained the Hannum age clock and devised a multi-class cancer signature using the autoML framework Autogluon both with and without reducing data dimensions by PCA. <br>
* There was a clear loss in predictive accuracy in the PC versions of the clock, which is a (simple) penalized regression model, and the (richer) cancer prediction model, which is an ensemble of several deep learning models. <br>
* We also modeled how training size and dimension reduction by PCA may affect the age prediction model accuracy (RMSE/slope) and repeatability (ICC) in external test sets, for Elastic net regression and 10K-fold cross validation. 

## Data availability
All microarray data is available freely at NCBI GEO or at the European Genome-phenome Archive (EGA) database under restricted access.
For accession IDs please refer to /0-data/tables/SupplementaryTable1.csv.

