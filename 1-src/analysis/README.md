
Corresponding pheno-files describing the samples in each dataset are given in 0-data/pheno.

## Preprocessing in house generated data (EPIC)
These are example scripts of our in house pipeline for data preprocessing.
Raw IDATs can be downloaded using the EGA access ID's given in 0-data/tables/SupplementaryTable1.csv.
Scripts for downstream analysis assume the beta-matrices outputted by the preprocessing pipeline are stored in 0-data/beta folder:

* beta_3CDisc.Rds:
    - 1657 cervical smear samples
    - training clocks as a function of training set size (Figure 2)
    - training a multiclass predictor for breast, ovarian and endometrial cancer (Extended data item 3)
* beta_3CExtVal.Rds:
    - 449 cervical smear samples
    - testing the accuracy of PC versus non-PC epigenetic biomarkers as a function of training set size (Figure 2, Extended data item 3)
* beta_repeatability.Rds:
    - 4 x 4 technical replicates cervical smear samples
    - testing the repeatability of PC versus non-PC clocks as a function of training set size (Figure 2)

## Preparing external 450K datasets 
Following 450K datasets (preprocessed beta-matrices) were downloaded from GEO and saved as .Rds files (data frames):

* GSE40279 (beta_Hannum.Rds): 
    - 656 blood samples
    - training the Hannum_PC clock (Figure 1) 
    - testing the accuracy of PC versus non-PC clocks as a function of training set size (Figure 2)
* GSE55763 (beta_BloodRep_450K.Rds): 
    - 2639 blood samples
    - testing the accuracy of the original and PC version of the Hannum clock on blood samples (Figure 1)
    - training clocks as a function of training set size (Figure 2)
* GSE55763 (beta_BloodRep_450K.Rds):
    - 2x36 technical replicates blood samples
    - testing the repeatability of the clocks' predictions on blood samples (Figure 1,2)




