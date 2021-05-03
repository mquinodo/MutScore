# MutScore

This repository contains the scripts used to compute MutScore, do figures and analysis included in the manuscript and for the two Shiny apps: MutScore and mutScore-batch.

The Shiny apps can be found at: 
https://iob-genetic.shinyapps.io/mutscore/
https://iob-genetic.shinyapps.io/mutscore-batch/

## Important note

The scripts used from ANNOVAR and the variants from HGMD database are not present in this repository due to protection of the data by licences.

## Prerequisites
+ BCFTools [[Link](https://samtools.github.io/bcftools/howtos/install.html)] (>= v1.9)
+ BEDTools [[Link](https://bedtools.readthedocs.io/en/latest/content/installation.html)] (>= v2.25.0)
+ Perl [[Link](https://www.perl.org/get.html)] (>= v5.22.0)
+ R [[Link](https://cran.r-project.org/mirrors.html)] (>= v3.2.0)
+ Picard [[Link](https://broadinstitute.github.io/picard/)]

## Installation
The scripts do not require compilation.

## Scripts

+ 00_databases-preprocessing.sh

   -> Download and processing of databases (ClinVar, ANNOVAR, dbNSFP4.0, CONDEL, ClinPred, OMIM, gnomAD and pext scores).

+ 01_clinvar_v0.8-all.sh

   -> Annotation of ClinVar, dbNSFP4.0 and gnomAD variants and parsing of the files.

+ 02_mutscore-computation-plots.R

   -> Processing of the annotated variants: computation of the amino-acid change and positional scores.
   -> Computation of MutScore model and scores.
   -> Analysis and figures used in the manuscript.
   
+ 03_cross-validation.R

   -> Cross-validation on the training set with figure used in the manuscript.
   
+ 04_MutScore-annotation.sh

   -> Annotation of all possible missense with MutScore.
   
+ 05_LiftOver.sh

   -> Conversion of scores in hg19 to hg38.
   
+ 06_data-for-mutland.R

   -> Processing of the data for MutLand plots.
   
+ 07_Shiny_preprocesssing.R

   -> Preprocessing of data for the two Shiny apps.
   
