# MutScore

This repository contains the scripts used to compute MutScore, do figures and analysis included in the manuscript and for the two Shiny apps: MutScore and mutScore-batch.

The Shiny apps can be found at: 
https://iob-genetic.shinyapps.io/mutscore/
https://iob-genetic.shinyapps.io/mutscore-batch/

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
-> Bla

+ 02_mutscore-computation-plots.R
+ 03_cross-validation.R
+ 04_MutScore-annotation.sh
+ 05_LiftOver.sh
+ 06_data-for-mutland.R
+ 07_Shiny_preprocesssing.R

