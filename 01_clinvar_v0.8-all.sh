#!/bin/bash

### Made by Mathieu Quinodoz
### February 2021

############# take different type of variants (PLP, BLB, VUS and conflicting) from clinvar vcf

here=/home/mquinodo/SYNO/WES/Clinvar2

mkdir -p $here/01_ClinVar-VCF

# download ClinVar VCF file(s) from their website: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/
# and put it here: $here/01_ClinVar-VCF/

# for training set
batch=20201121

gunzip $here/01_ClinVar-VCF/clinvar_$batch.vcf.gz

mkdir -p $here/clinvar-$batch-PLP/variants $here/clinvar-$batch-BLB/variants $here/clinvar-$batch-VUS/variants $here/clinvar-$batch-CON/variants

grep -P "CLNSIG=Pathogenic;|CLNSIG=Likely_pathogenic;|CLNSIG=Pathogenic/Likely_pathogenic;" $here/01_ClinVar-VCF/clinvar_$batch.vcf | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $6 "\t" $8}' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=4.079;DP=24;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=70.00;MQRankSum=0.000;QD=15.85;ReadPosRankSum=0.810;SOR=0.481\tGT:AD:DP:GQ:PL\t0/1:13,10:23:99:372,0,445"}' > $here/clinvar-$batch-PLP/variants/clinvar-$batch-PLP.variants.haplotype.vcf 
grep -P "CLNSIG=Likely_benign;|CLNSIG=Benign;|CLNSIG=Benign/Likely_benign;" $here/01_ClinVar-VCF/clinvar_$batch.vcf | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $6 "\t" $8}' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=4.079;DP=24;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=70.00;MQRankSum=0.000;QD=15.85;ReadPosRankSum=0.810;SOR=0.481\tGT:AD:DP:GQ:PL\t0/1:13,10:23:99:372,0,445"}' > $here/clinvar-$batch-BLB/variants/clinvar-$batch-BLB.variants.haplotype.vcf 
grep -P "CLNSIG=Uncertain_significance;" $here/01_ClinVar-VCF/clinvar_$batch.vcf | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $6 "\t" $8}' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=4.079;DP=24;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=70.00;MQRankSum=0.000;QD=15.85;ReadPosRankSum=0.810;SOR=0.481\tGT:AD:DP:GQ:PL\t0/1:13,10:23:99:372,0,445"}' > $here/clinvar-$batch-VUS/variants/clinvar-$batch-VUS.variants.haplotype.vcf 
grep -P "CLNSIG=Conflicting_interpretations_of_pathogenicity;" $here/01_ClinVar-VCF/clinvar_$batch.vcf | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $6 "\t" $8}' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=4.079;DP=24;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=70.00;MQRankSum=0.000;QD=15.85;ReadPosRankSum=0.810;SOR=0.481\tGT:AD:DP:GQ:PL\t0/1:13,10:23:99:372,0,445"}' > $here/clinvar-$batch-CON/variants/clinvar-$batch-CON.variants.haplotype.vcf 

gzip $here/01_ClinVar-VCF/clinvar_$batch.vcf

# for testing set 1
batch=20210404

gunzip $here/01_ClinVar-VCF/clinvar_$batch.vcf.gz

mkdir -p $here/clinvar-$batch-PLP/variants $here/clinvar-$batch-BLB/variants $here/testing-set-1

grep -P "CLNSIG=Pathogenic;|CLNSIG=Likely_pathogenic;|CLNSIG=Pathogenic/Likely_pathogenic;" $here/01_ClinVar-VCF/clinvar_$batch.vcf | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $6 "\t" $8}' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=4.079;DP=24;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=70.00;MQRankSum=0.000;QD=15.85;ReadPosRankSum=0.810;SOR=0.481\tGT:AD:DP:GQ:PL\t0/1:13,10:23:99:372,0,445"}' > $here/clinvar-$batch-PLP/variants/clinvar-$batch-PLP.variants.haplotype.vcf 
grep -P "CLNSIG=Likely_benign;|CLNSIG=Benign;|CLNSIG=Benign/Likely_benign;" $here/01_ClinVar-VCF/clinvar_$batch.vcf | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $8 "\t" $6 "\t" $8}' | awk -F"\t" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=4.079;DP=24;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=70.00;MQRankSum=0.000;QD=15.85;ReadPosRankSum=0.810;SOR=0.481\tGT:AD:DP:GQ:PL\t0/1:13,10:23:99:372,0,445"}' > $here/clinvar-$batch-BLB/variants/clinvar-$batch-BLB.variants.haplotype.vcf 

awk -F"\t" '{print $1 "-" $2 "-" $4 "-" $5}' $here/clinvar-$batch-PLP/variants/clinvar-$batch-PLP.variants.haplotype.vcf | sed -e 's/chr//g' > $here/testing-set-1/clinvar-$batch-PLP.pos.tsv
awk -F"\t" '{print $1 "-" $2 "-" $4 "-" $5}' $here/clinvar-$batch-BLB/variants/clinvar-$batch-BLB.variants.haplotype.vcf | sed -e 's/chr//g' > $here/testing-set-1/clinvar-$batch-BLB.pos.tsv

gzip $here/01_ClinVar-VCF/clinvar_$batch.vcf

rm -rf $here/clinvar-$batch-PLP $here/clinvar-$batch-BLB


############# annotate ClinVar, gnomAD and dbNFSP with ANNOVAR

version=clinvar_annotation
script=$here/00_scripts/$version
batch=20201121

for file in clinvar-$batch-BLB clinvar-$batch-PLP clinvar-$batch-VUS clinvar-$batch-CON dbNFSP4.0 gnomAD-all
do
  outdir=$here/$file/$version
  mkdir -p $outdir
  vcf=$here/$file/variants/$file.variants.haplotype.vcf
  nice nohup bash $script/annotation_wo_mutscore.sh $here $vcf $file $script $outdir > $outdir/annotation.log.tsv &
done


############# extract variants after annotation by isoforms for ClinVar

batch=20201121

for type in clinvar-$batch-PLP clinvar-$batch-VUS clinvar-$batch-CON clinvar-$batch-BLB
do
  echo $type
  file=$here/$type/clinvar_annotation/$type.avinput.exonic_splicing.clean.annovar
  grep -v -P "^#|BRCA2:NM_000059.3:exon10:c.1888dupA:p.L629fs" $file.tsv | grep -P "^nonsynonymous|^frameshift|^stopgain|^synonymous" | cut -f1-8,10-54 | awk -F"\t" -v type="$type" '{printf $1; for (j = 2; j <= 9; j++) {printf "\t" $j}; printf "\t.\t.\t."; for (j = 10; j <= NF; j++) {printf "\t" $j}; print ""}' | awk -F"\t" -v type="$type" '{n=split($3,a,","); for (i = 1; i < n; i++) {split(a[i],b,":"); printf type "\t" b[2] "\t" b[4] "\t" b[5] "\t" $1 "\t" b[1] "\t" b[3]; for (j = 3; j <= NF; j++) {printf "\t" $j}; print ""}}' > $file.$type.isoform-info.temp.tsv
  grep -P "^splicing|^intronic|^UTR" $file.tsv | cut -f1-8,10-54 | awk -F"\t" -v type="$type" '{printf $1; for (j = 2; j <= 9; j++) {printf "\t" $j}; printf "\t.\t.\t."; for (j = 10; j <= NF; j++) {printf "\t" $j}; print ""}' | awk -F"\t" -v type="$type" '{printf type "\t.\t.\t.\t"; printf $1 "\t" $2 "\t."; for (j = 3; j <= NF; j++) {printf "\t" $j}; print ""}' > $file.$type.isoform-info.temp2.tsv
  cat $file.$type.isoform-info.temp.tsv $file.$type.isoform-info.temp2.tsv  > $file.$type.isoform-info.$batch.tsv
  rm $file.$type.isoform-info.temp2.tsv $file.$type.isoform-info.temp.tsv
done

# gene list of genes with PLP variants
mkdir -p $here/data-$batch
type=clinvar-$batch-PLP
cat $here/$type/clinvar_annotation/$type.avinput.exonic_splicing.clean.annovar.$type.isoform-info.$batch.tsv | cut -f6 | sort | uniq | awk -F"\t" '{print "\t" $1 "\t"}' > $here/data-$batch/gene-list.$batch.tsv


############# extract variants after annotation by isoforms for gnomAD and dbNFSP

batch=20201121

type=gnomAD-all
file=$here/$type/clinvar_annotation/$type.avinput.exonic_splicing.clean.annovar
grep -v "^#" $file.tsv | grep -P "^nonsynonymous|^frameshift|^stopgain|^synonymous" | cut -f1-8,10-54 | awk -F"\t" -v type="$type" '{printf $1; for (j = 2; j <= 9; j++) {printf "\t" $j}; printf "\t.\t.\t."; for (j = 10; j <= NF; j++) {printf "\t" $j}; print ""}' | awk -F"\t" -v type="$type" '{n=split($3,a,","); for (i = 1; i < n; i++) {split(a[i],b,":"); printf type "\t" b[2] "\t" b[4] "\t" b[5] "\t" $1 "\t" b[1] "\t" b[3]; for (j = 3; j <= NF; j++) {printf "\t" $j}; print ""}}' > $file.$type.isoform-info.temp.tsv
grep -P "^splicing|^intronic|^UTR" $file.tsv | cut -f1-8,10-54 | awk -F"\t" -v type="$type" '{printf $1; for (j = 2; j <= 9; j++) {printf "\t" $j}; printf "\t.\t.\t."; for (j = 10; j <= NF; j++) {printf "\t" $j}; print ""}' | awk -F"\t" -v type="$type" '{printf type "\t.\t.\t.\t"; printf $1 "\t" $2 "\t."; for (j = 3; j <= NF; j++) {printf "\t" $j}; print ""}' > $file.$type.isoform-info.temp2.tsv
cat $file.$type.isoform-info.temp.tsv $file.$type.isoform-info.temp2.tsv > $file.$type.isoform-info.$batch.tsv
rm $file.$type.isoform-info.temp2.tsv $file.$type.isoform-info.temp.tsv

type=dbNFSP4.0
file=$here/$type/clinvar_annotation/$type.avinput.exonic_splicing.clean.annovar
grep -v "^#" $file.tsv | grep -P "^nonsynonymous" | cut -f1-8,10-54 | awk -F"\t" -v type="$type" '{printf $1; for (j = 2; j <= 9; j++) {printf "\t" $j}; printf "\t.\t.\t."; for (j = 10; j <= NF; j++) {printf "\t" $j}; print ""}' | awk -F"\t" -v type="$type" '{n=split($3,a,","); for (i = 1; i < n; i++) {split(a[i],b,":"); printf type "\t" b[2] "\t" b[4] "\t" b[5] "\t" $1 "\t" b[1] "\t" b[3]; for (j = 3; j <= NF; j++) {printf "\t" $j}; print ""}}' > $file.$type.isoform-info.$batch.tsv


############# merge gnomAD and dbNFSP with ClinVar variants (after parsing of ClinVar information)

batch=20201121

mkdir -p $here/data-$batch
mkdir -p $here/plots-$batch

gene=ALL
cat $here/*/clinvar_annotation/*.isoform-info.$batch.tsv | grep -v wholegene | cut -f1-7,9- > $here/data-$batch/$gene.temp.tsv
cut -f1-12,14- $here/data-$batch/$gene.temp.tsv > $here/data-$batch/$gene.temp1.tsv
cut -f13 $here/data-$batch/$gene.temp.tsv | awk -F"\t" '{n=split($1,a,";"); c=0; for (i = 1; i <= n; i++) {if(a[i] ~ /^CLNREVSTAT/) {n=split(a[i],b,"="); print b[2]; c=1}  } if( c == 0) print "NA" }' > $here/data-$batch/$gene.REV.tsv
cut -f13 $here/data-$batch/$gene.temp.tsv | awk -F"\t" '{n=split($1,a,";"); c=0; for (i = 1; i <= n; i++) {if(a[i] ~ /^ORIGIN/) {n=split(a[i],b,"="); print b[2]; c=1}  } if( c == 0) print "NA" }' > $here/data-$batch/$gene.ORI.tsv
cut -f13 $here/data-$batch/$gene.temp.tsv | awk -F"\t" '{n=split($1,a,";"); c=0; for (i = 1; i <= n; i++) {if(a[i] ~ /^CLNDISDB/) {n=split(a[i],b,"="); print b[2]; c=1}  } if( c == 0) print "NA" }' > $here/data-$batch/$gene.temp3.tsv
cat $here/data-$batch/$gene.temp3.tsv | awk -F"\t" '{n=split($1,a,","); c=0; for (i = 1; i <= n; i++) {if(a[i] ~ /^OMIM/) {n=split(a[i],b,":"); if( c == 1) {printf "," b[2]}; if( c == 0) {printf b[2]}; c=1}  } if( c == 0) {printf "NA"} print ""}' > $here/data-$batch/$gene.temp4.tsv
cat $here/data-$batch/$gene.temp4.tsv | awk -F"\t" '{n=split($1,a,"|"); print a[1]}' > $here/data-$batch/$gene.OMIM.tsv
paste -d"\t" $here/data-$batch/$gene.temp1.tsv $here/data-$batch/$gene.REV.tsv $here/data-$batch/$gene.OMIM.tsv $here/data-$batch/$gene.ORI.tsv > $here/data-$batch/$gene.tsv
rm $here/data-$batch/$gene.temp*.tsv $here/data-$batch/$gene.REV.tsv $here/data-$batch/$gene.OMIM.tsv $here/data-$batch/$gene.ORI.tsv
head $here/data-$batch/$gene-$batch.tsv > $here/data-$batch/$gene.head.tsv
mkdir -p $here/Shiny


