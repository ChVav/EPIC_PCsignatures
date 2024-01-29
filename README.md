
This version of the code accompanies a pending submission to BioData Mining, Brief Report: Vavourakis, C. D., Herzog, C. M., & Widschwendter, M. (2023). Devising reliable and accurate epigenetic predictors: choosing the optimal computational solution (submitted). <br>

## Scope and data availability

* I retrained the Hannum age clock and devised a multi-class cancer signature using the autoML framework Autogluon to test the effect of reducing data dimensions by PCA on the performance of the respective simple and rich prediction models. 
For PC versions of the clock and cancer signature I used the approach by Higgins-Chen et al (https://doi.org/10.1038/s43587-022-00248-2) for bolstering reliability of epigenetic clocks predicting age from methylation data.
There was a clear loss in predictive accuracy in the PC versions of the clock and cancer prediction model compared to the models trained in the original feature space. <br>
* I also modeled how training size and dimension reduction by PCA may affect the age prediction model accuracy (RMSE/slope) and repeatability (ICC) in external test sets, for Elastic net regression and 10K-fold cross validation. 
This was done both for 450K data derived from blood samples and EPIC v1.0. data derived from cervical smears. <br>
* To investigate precision and repeatability of trained prediction models for EPIC v.1.0 DNA methylation data, a test data set with technical quadruplicates derived from blood and cervical smear samples from four different women was generated. 
Fresh, and frozen blood were treated as separate tissue types, so in total 4 * 4 * 3 =  48 microarrays were analyzed. <br>
The original raw microarray data will be made available at the European Genome-phenome Archive (EGA) database under restricted access (EGAS00001007184), once a corresponding manuscript is accepted <br>
All other training and test data are already available (under restricted access) in public repositories, see /0-data/tables/SupplementaryTable1.csv.

