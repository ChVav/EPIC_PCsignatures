
Pheno-files describing the samples in each dataset are given in 0-data/pheno.

## Preprocessing in house generated data (EPIC)
These are example scripts of our in house pipeline for data preprocessing.
In addition, probes not on EPIC v2.0 (see 0-data/flagged_probes/notonEPICv2_names.csv) and probes with a consistently low intensity normalized intensity across in house test sets (see 0-data/flagged_probes/lowMI_names.csv) were removed.

Raw IDATs can be downloaded using the access ID's given in 0-data/tables/SupplementaryTable1.csv.

Scripts for downstream analysis assume the beta-matrices outputted by the preprocessing pipeline are stored in 0-data/beta folder.

## Preparing external datasets (Hannum, 450K)

### Training Hannum_PC
For training the Hannum_PC clock, the preprocessed beta-matrix was downloaded from GEO and extracted as follows in R:

```
gunzip("<path-to-download>/GSE40279/GSE40279_average_beta.txt.gz") # should be 656 samples, age range 19-101
beta <- read.table(file = "<path-to-download>/GSE40279/GSE40279_average_beta.txt", header =TRUE, sep = "\t")
saveRDS(beta,file = "/0-data/beta/beta_Hannum.Rds")
```

### Testing different versions Hannum clock

Also for comparing the performance of the original vs the PC version of the Hannum clock, the preprocessed beta-matrix was downloaded from GEO, in R:

```
gunzip("<path-to-download-folder>/GSE55763/GSE55763_normalized_betas.txt.gz")

# make list of samples needed for training
samples <- read.RDS("/0-data/pheno/pheno_BloodRep_450K.Rds") %>% pull(basename)
write.csv(samples, file="<path-to-download-folder>/GSE55763/samples.csv", quote = FALSE, row.names=FALSE)
```

The relevant samples were then extracted, so that a more reasonably sized dataset could be read in with R:

```
#!/bin/bash/
sed -i '1,1d' samples.csv  #remove first line from file
cut -f"$({ echo 1; grep -Fxn "$(<samples.csv)" < <(head -1 GSE55763_normalized_betas.txt | tr '\t' '\n'); } | cut -f1 -d: | paste -sd,)" GSE55763_normalized_betas.txt > blood_subset.txt
```

```
beta <- read.table(file = "<path-to-download-folder>/GSE55763/blood_subset.txt", header =TRUE, sep = "\t") # 473864 CpGs
colnames(beta) <- gsub("X","",colnames(beta)) # remove trailing X in columnames
beta <- beta %>% column_to_rownames(var = "ID_REF")

# check samples in same order as pheno and save
beta <- beta[,match(pheno$basename,colnames(beta))]
identical(pheno$basename, colnames(beta)) 

saveRDS(beta,file = "/0-data/beta/beta_BloodRep_450K.Rds")
```



