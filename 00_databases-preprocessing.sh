#!/bin/bash

### Made by Mathieu Quinodoz
### March 2020

hereMain=/home/mquinodo/SYNO/WES/Clinvar3

#### ClinVar

here=$hereMain
# VCF with variants from ClinVar can be downloaded here: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/
# It has be put in the following directory:
mkdir -p $here/01_ClinVar-VCF


#### refGene and dbscsnv11

here=$hereMain/00_scripts

perl $here/clinvar_annotation/annotate_variation.pl -webfrom annovar -downdb refGeneWithVer -buildver hg19 $here/clinvar_annotation/humandb/

perl $here/clinvar_annotation/annotate_variation.pl -webfrom annovar -downdb dbscsnv11 -buildver hg19 $here/clinvar_annotation/humandb/


#### dbNSFP4.0

here=$hereMain/00_scripts/databases-processing

wget -P $here ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.0a.zip

unzip $here/dbNSFP4.0a.zip

for file in dbNSFP4.0a_variant.chr1 dbNSFP4.0a_variant.chr13  dbNSFP4.0a_variant.chr17  dbNSFP4.0a_variant.chr20  dbNSFP4.0a_variant.chr4  dbNSFP4.0a_variant.chr8  dbNSFP4.0a_variant.chrY dbNSFP4.0a_variant.chr10  dbNSFP4.0a_variant.chr14  dbNSFP4.0a_variant.chr18  dbNSFP4.0a_variant.chr21  dbNSFP4.0a_variant.chr5  dbNSFP4.0a_variant.chr9 dbNSFP4.0a_variant.chr11  dbNSFP4.0a_variant.chr15  dbNSFP4.0a_variant.chr19  dbNSFP4.0a_variant.chr22  dbNSFP4.0a_variant.chr6  dbNSFP4.0a_variant.chrM dbNSFP4.0a_variant.chr12  dbNSFP4.0a_variant.chr16  dbNSFP4.0a_variant.chr2   dbNSFP4.0a_variant.chr3   dbNSFP4.0a_variant.chr7  dbNSFP4.0a_variant.chrX
do
  gunzip $here/$file.gz
  cut -f1-9,37-84,86-155,229,302 $here/$file | awk -F "\t" '{ printf $1; for(i=2;i<=9;++i) printf "\t" $i; for(i=10;i<=NF;++i) {if($i ~ /;/) {split($i,a,";"); m=0; c=0; for (j in a) {if(a[j]!="."){m=m+a[j]; c=c+1;}}; max=m/(c+0.000000001); printf "\t" max} else {printf "\t" $i} } print ""}' | awk -F"\t" '{printf $8 "\t" $9 "\t" $3 "\t" $4 ;for(i=10;i<=NF;++i) printf "\t" $i; print ""}' | sed -e 's/.;//g' | cut -f1-5,8,11,14,17,21,26,29,32,35,37,40,44,47,49,53,55,57,60,71,72,74,78,81,84,87,89,103,104,106,108,110,112,114,116,119,123,124 | awk -F"\t" '{printf $1 "\t" $2 "\t" $2 ; for(i=3;i<=NF;++i) printf "\t" $i; print ""}' > $here/$file.sel.tsv
  rm $here/$file
done

cat $here/dbNSFP4.0a_variant.chr*.sel.tsv | grep -v -P "hg19_chr" | sort -k 1,1 -k 2,2n -k 3,3n > $here/dbNSFP4.0a.all-sel.tsv
cat $here/head-dbNSFP4.0.txt $here/dbNSFP4.0a.all-sel.tsv | awk -F"\t" '{if($1!=".") print $0}'  > $here/../clinvar_annotation/humandb/hg19_dbNSFP4.0a.txt
perl $here/../clinvar_annotation/annovar_idx.pl $here/../clinvar_annotation/humandb/hg19_dbNSFP4.0a.txt 100 > $here/../clinvar_annotation/humandb/hg19_dbNSFP4.0a.txt.idx
rm $here/dbNSFP4.0a_variant.chr*.sel.tsv $here/dbNSFP4.0a.all-sel.tsv $here/dbNSFP4.0a.all-sel.tsv

here=$hereMain

mkdir -p $here/dbNFSP4.0/variants
awk -F"\t" '{if($19!=".") print $0}' $here/00_scripts/clinvar_annotation/humandb/hg19_dbNSFP4.0a.txt | cut -f1-5,19 > $here/dbNFSP4.0/variants/temp.tsv
awk -F"\t" '{print $1 "\t" $2 "\t" "." "\t" $4 "\t" $5 "\t.\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=1.242;DP=8;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=28.81;MQRankSum=0.319;QD=21.45;ReadPosRankSum=-0.319;SOR=1.863\tGT:AD:DP:GQ:PL\t0/1:2,6:8:41:179,0,41"}' $here/dbNFSP4.0/variants/temp.tsv > $here/dbNFSP4.0/variants/dbNFSP4.0.variants.haplotype.vcf 
rm $here/dbNFSP4.0/variants/temp.tsv


#### CONDEL

here=$hereMain/00_scripts/databases-processing

wget -P $here http://bbglab.irbbarcelona.org/fannsdb/downloads/fannsdb.tsv.gz

gunzip $here/fannsdb.tsv.gz

tail -n+2 $here/fannsdb.tsv | awk -F"\t" '{if($15>(-1)){ print $1 "\t" $2 "\t" $2 "\t" $4 "\t" $5 "\t" $15}}' | uniq > $here/../clinvar_annotation/humandb/hg19_CONDEL.txt
perl $here/../clinvar_annotation/annovar_idx.pl $here/../clinvar_annotation/humandb/hg19_CONDEL.txt 100 > $here/../clinvar_annotation/humandb/hg19_CONDEL.txt.idx


#### ClinPred

# download the file ClinPred.txt here:
# https://drive.google.com/file/d/1Rh696VtfU4BBI33EltsW5BltOCDqc8Q4/view
# and put it in databases-processing folder

here=$hereMain/00_scripts/databases-processing

tail -n+2 $here/ClinPred.txt | awk -F"\t" '{ print $1 "\t" $2 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' > $here/../clinvar_annotation/humandb/hg19_ClinPred.txt
perl $here/../clinvar_annotation/annovar_idx.pl $here/../clinvar_annotation/humandb/hg19_ClinPred.txt 100 > $here/../clinvar_annotation/humandb/hg19_ClinPred.txt.idx


#### OMIM

here=$hereMain/00_scripts/databases-processing

# genemap2.txt has to be downloaded from OMIM after registration here: https://omim.org/downloads
# The file has to be put in databases-processing folder

file=$here/genemap2.txt

grep -v -P "^#|^chrX|^chrY" $file | awk -F"\t" '{if($3-$2<1000000) print $0}' | cut -f7,13 | awk -F"\t" '{if($2!=""){n=split($1,var,", "); for (i = 1; i <=n; i++) {print var[i] "\t" $2} }}' | \
awk -F"\t" '{n=split($2,var,"; "); for (i = 1; i <=n; i++) {print $1 "\t" var[i]} }' | \
grep -v -P "\t\?" | grep -P "Autosomal dominant|Autosomal recessive|somatic" > $here/gene-disease.tsv

grep -P "Autosomal dominant|somatic" $here/gene-disease.tsv | cut -f1 | sort | uniq > $here/gene-dom.tsv
grep -P "Autosomal recessive" $here/gene-disease.tsv | cut -f1 | sort | uniq > $here/gene-rec.tsv

comm -13 $here/gene-dom.tsv $here/gene-rec.tsv > $here/gene-rec-pure.tsv
comm -23 $here/gene-dom.tsv $here/gene-rec.tsv > $here/gene-dom-pure.tsv
comm -12 $here/gene-dom.tsv $here/gene-rec.tsv > $here/gene-both.tsv

rm $here/gene-disease.tsv $here/gene-dom.tsv $here/gene-rec.tsv


#### gnomAD 

here=$hereMain/00_scripts/databases-processing

wget -P $here https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/exomes/gnomad.exomes.r2.1.1.sites.vcf.bgz
wget -P $here https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz

mv $here/gnomad.exomes.r2.1.1.sites.vcf.bgz $here/gnomad.exomes.vcf.gz
mv $here/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz $here/gnomad.genomes.vcf.gz

gunzip $here/gnomad.exomes.vcf.gz
gunzip $here/gnomad.genomes.vcf.gz

for file in $here/gnomad.genomes $here/gnomad.exomes
do
	grep -P "\tPASS\t" $file.vcf | cut -f1-5 > $file.cut.vcf
	awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t100\t.\tAC=1;AF=0.500;AN=2;BaseQRankSum=1.242;DP=8;ExcessHet=3.0103;FS=0.000;MLEAC=1;MLEAF=0.500;MQ=28.81;MQRankSum=0.319;QD=21.45;ReadPosRankSum=-0.319;SOR=1.863\tGT:AD:DP:GQ:PL\t0/1:2,6:8:41:179,0,41"}' $file.cut.vcf | grep -P "^chr" > $file.cut.add.vcf
	grep -v "#" $file.cut.add.vcf > $file.cut.add.final.vcf
	rm $file.cut.add.vcf $file.cut.vcf
done

here=$hereMain

mkdir -p $here/gnomAD-all/variants
cat $here/00_scripts/databases-processing/gnomad.genomes.cut.add.vcf $here/00_scripts/databases-processing/gnomad.exomes.cut.add.vcf | sort | uniq > $here/gnomAD-all/variants/gnomAD-all.variants.haplotype.vcf
rm $here/00_scripts/databases-processing/gnomad.*.cut.add.final.vcf $here/gnomad.genomes.vcf $here/gnomad.exomes.vcf


#### pext scores

here=$hereMain/00_scripts

wget -P $here/databases-processing https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-tx-annotation/gnomad_browser/all.baselevel.021620.tsv.bgz

file=$here/databases-processing/all.baselevel.021620

mv $file.tsv.bgz $file.tsv.gz
gunzip $file.tsv.gz

sed -e 's/NaN/0/g' $file.tsv > $file.zero.tsv
awk '{printf $2 "\t"; split($2,a,":"); printf a[1] "\t" a[2] "\t" a[2] "\t" $57 "\t"; b=0; for (i=4;i<=56;i++){if($i>b){b=$i}} print b}' $file.zero.tsv > $file.short.tsv
sort -k1,1 -k5,5nr $file.short.tsv | sort -u -k1,1 | cut -f2-6 | sort -k1,1n -k2,2n > $file.short.unique.tsv
tail -n+2 $file.short.unique.tsv | awk -F"\t" '{print "chr" $1 "\t" $2 "\t" $3 "\t" sprintf("%.9f", $4) "\t" sprintf("%.5f", $5)}' | sort -V -k1,1 -k2,2 > $here/databases-processing/hg19_pext-temp.txt
awk -F"\t" '{nuc="A_T_C_G"; split(nuc,nucs,"_"); for (i = 1; i <=4; i++){for (j = 1; j <=4; j++) { print $1 "\t" $2 "\t" $3 "\t" nucs[i] "\t" nucs[j] "\t" $4 "\t" $5}}}' $here/databases-processing/hg19_pext-temp.txt > $here/databases-processing/hg19_pext2.txt
cat $here/databases-processing/head-pext.tsv $here/databases-processing/hg19_pext2.txt > $here/clinvar_annotation/humandb/hg19_pext.txt

perl $here/clinvar_annotation/annovar_idx.pl $here/clinvar_annotation/humandb/hg19_pext.txt 1000 > $here/clinvar_annotation/humandb/hg19_pext.txt.idx

rm $file.zero.tsv $file.short.tsv $file.short.unique.tsv $here/databases-processing/hg19_pext2.txt $here/databases-processing/hg19_pext-temp.txt

