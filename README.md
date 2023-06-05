
Code accompanying the draft paper "Devising reliable and accurate epigenetic clocks: choosing the optimal computational solution". <br>

## Scope and data availability
I retrained the Hannum age clock and devised a multi-class cancer signature using the autoML framework Autogluon to test the effect of reducing data dimensions by PCA on the performance of the respective simple and rich prediction models. For PC versions of the clock and cancer signature I used the approach by Higgins-Chen et al (https://doi.org/10.1038/s43587-022-00248-2) for bolstering reliability of epigenetic clocks predicting age from methylation data. <br>
I also modeled how the training size (n=32-1657) affects the age prediction model accuracy (RMSE) and repeatability (ICC) in an external test set, for a simple Elastic net regression model. <br>
To investigate precision and repeatability of trained prediction models, an EPIC v.1.0 DNA methylation dataset with technical quadruplicates derived from blood and cervical smear samples from four different women was generated. Fresh, and frozen blood were treated as separate tissue types, so in total 4 * 4 * 3 =  48 microarrays were analyzed. <br>
The original raw microarray data will be made available at the European Genome-phenome Archive (EGA) database under restricted access (EGAS00001007184). <br>
All other training and test data are already available (under restricted access) in public repositories, see /0-data/tables/SupplementaryTable1.csv.

