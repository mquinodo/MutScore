#!/bin/bash

### Made by Mathieu Quinodoz
### January 2021

here=$1
vcf=$2
file=$3
script=$4
outdir=$5

######################################################################################################

path=$outdir

# Delete variants with multiple alternative alleles
grep -v "##" $vcf | egrep -v "1/2"  > $path/$file.tsv
name=$file

# convert to annovar input
perl $script/convert2annovar.pl -format vcf4 $path/$name.tsv -outfile $path/$name.avinput -include -withzyg

# Prediction software information
cd $script

perl $script/table_annovar.pl $path/$name.avinput $here/00_scripts/clinvar_annotation/humandb -intronhgvs 1000000 -protocol ClinPred,dbscsnv11,CONDEL,dbNSFP4.0a,pext -operation f,f,f,f,f -build hg19 -nastring . -remove -polish 

# treating data of annovar for correct annotation of indels (not same annotation)
tail -n+2 $path/$name.avinput.hg19_multianno.txt | cut -f6- > $path/$name.avinput.hg19_multianno.noH.txt
head -n1 $path/$name.avinput.hg19_multianno.txt > $path/$name.avinput.hg19_multianno.header.txt
awk -F"\t" '{print $9 "\t" $10 "\t" $10 "\t" $12 "\t" $13}' $path/$name.avinput > $path/$name.avinput.short
paste -d "\t" $path/$name.avinput.short $path/$name.avinput.hg19_multianno.noH.txt > $path/$name.avinput.temp
cat $path/$name.avinput.hg19_multianno.header.txt $path/$name.avinput.temp > $path/$name.avinput.hg19_multianno.coord.txt
rm $path/$name.avinput.hg19_multianno.noH.txt $path/$name.avinput.hg19_multianno.header.txt $path/$name.avinput.short $path/$name.avinput.temp
nameAnno=$file.avinput.hg19_multianno.coord.txt
name=$file.avinput

# refGene annotation
perl $script/annotate_variation.pl -buildver hg19 -hgvs -dbtype refGeneWithVer $path/$name $script/humandb/ 

# merge files
perl $script/merge_exonic_and_splicing_all.pl $path/$name.exonic_variant_function $path/$name.exonic_splicing.tsv $path/$name.variant_function
file_exon=$path/$name.exonic_splicing

# clean file
perl $script/clean_exonic_output.pl ${file_exon}.tsv ${file_exon}.clean.tsv
file_exon=${file_exon}.clean

# add annovar
perl $script/add_annovar.pl ${file_exon}.tsv ${file_exon}.annovar.tsv ${path}/$nameAnno
file_exon=${file_exon}.annovar

# remove temporary files
rm $path/$file.avinput $file.tsv $path/$file.avinput.exonic_splicing.clean.tsv $path/$file.avinput.exonic_splicing.tsv $path/$file.avinput.exonic_variant_function
rm $path/$file.avinput.hg19_multianno.coord.txt $path/$file.avinput.invalid_input $path/$file.avinput.log $path/$file.avinput.variant_function $path/$file.avinput.hg19_multianno.txt


