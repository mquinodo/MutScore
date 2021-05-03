#!/bin/bash

### Made by Mathieu Quinodoz
### March 2021

hereMain=/home/mquinodo/SYNO/WES/Clinvar2

here=$hereMain/data-20201121

example=$hereMain/01_ClinVar-VCF/clinvar_20201121.vcf

gunzip $example.gz

grep -P "^#" $example > $here/temp-header.tsv

awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $5 "\t" $3 "\t" $4 "\t.\t.\tALLELEID=824438"}' $here/mutscore-hg19.tsv > $here/temp-body.tsv

cat $here/temp-header.tsv $here/temp-body.tsv > $here/temp.vcf 

java -jar /usr/local/bin/picard.jar LiftoverVcf \
I=$here/temp.vcf \
O=$here/temp-hg38.vcf \
CHAIN=$hereMain/00_scripts/databases-processing/hg19ToHg38.over.chain \
REJECT=$here/temp-rejected.vcf \
R=/home/mquinodo/GRCh38/GRCh38_full_analysis_set_plus_decoy_hla.fa

grep -v "#" $here/temp-hg38.vcf | sed -e 's/chr//g' | awk -F"\t" '{print $1 "\t" $2 "\t" $4 "\t" $5 "\t" $3}' > $here/mutscore-hg38.tsv

rm $here/temp-body.tsv $here/temp-header.tsv $here/temp.vcf $here/temp-rejected.vcf $here/temp-hg38.vcf $here/temp-hg38.vcf.idx

gzip $example

