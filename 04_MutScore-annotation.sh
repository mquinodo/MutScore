#!/bin/bash

### Made by Mathieu Quinodoz
### February 2021

################################################

hereMain=/home/mquinodo/SYNO/WES/Clinvar2

# index for ANNOVAR file
here=$hereMain/00_scripts/clinvar_annotation
perl $here/annovar_idx.pl $here/humandb/hg19_mutscore.txt 1000 > $here/humandb/hg19_mutscore.txt.idx

####### DO THE ANNOTATION for MutScore  ########

here=$hereMain

version=clinvar_annotation
script=$here/00_scripts/$version
batch=20201121

for file in dbNFSP4.0
do
  outdir=$here/$file/${version}-mutscore
  mkdir -p $outdir
  tempdir=$here/$file/${version}-mutscore/temp
  mkdir -p $tempdir
  vcf=$here/$file/variants/$file.variants.haplotype.vcf
  nice nohup bash $script/annotation_with_mutscore.sh $here $vcf $file $script $outdir > $outdir/annotation.log.tsv &
done

##

here=$hereMain

batch=20201121

rm -rf $here/$file/${version}-mutscore/temp

type=dbNFSP4.0
file=$here/$type/clinvar_annotation-mutscore/$type.avinput.exonic_splicing.clean.annovar
grep -v "^#" $file.tsv | grep -P "^nonsynonymous" | cut -f1-8,10-57 | awk -F"\t" -v type="$type" '{printf $1; for (j = 2; j <= 9; j++) {printf "\t" $j}; printf "\t.\t.\t."; for (j = 10; j <= NF; j++) {printf "\t" $j}; print ""}' | awk -F"\t" -v type="$type" '{n=split($3,a,","); for (i = 1; i < n; i++) {split(a[i],b,":"); printf type "\t" b[2] "\t" b[4] "\t" b[5] "\t" $1 "\t" b[1] "\t" b[3]; for (j = 3; j <= NF; j++) {printf "\t" $j}; print ""}}' > $file.$type.isoform-info.$batch.tsv

##

here=$hereMain

batch=20201121

mkdir -p $here/data-$batch

gene=ALL
cat $here/dbNFSP4.0/clinvar_annotation-mutscore/*.isoform-info.$batch.tsv | grep -v wholegene | cut -f1-7,9- > $here/data-$batch/$gene.temp.tsv
cut -f1-12,14- $here/data-$batch/$gene.temp.tsv > $here/data-$batch/$gene.temp1.tsv
cut -f13 $here/data-$batch/$gene.temp.tsv | awk -F"\t" '{n=split($1,a,";"); c=0; for (i = 1; i <= n; i++) {if(a[i] ~ /^CLNREVSTAT/) {n=split(a[i],b,"="); print b[2]; c=1}  } if( c == 0) print "NA" }' > $here/data-$batch/$gene.REV.tsv
cut -f13 $here/data-$batch/$gene.temp.tsv | awk -F"\t" '{n=split($1,a,";"); c=0; for (i = 1; i <= n; i++) {if(a[i] ~ /^ORIGIN/) {n=split(a[i],b,"="); print b[2]; c=1}  } if( c == 0) print "NA" }' > $here/data-$batch/$gene.ORI.tsv
cut -f13 $here/data-$batch/$gene.temp.tsv | awk -F"\t" '{n=split($1,a,";"); c=0; for (i = 1; i <= n; i++) {if(a[i] ~ /^CLNDISDB/) {n=split(a[i],b,"="); print b[2]; c=1}  } if( c == 0) print "NA" }' > $here/data-$batch/$gene.temp3.tsv
cat $here/data-$batch/$gene.temp3.tsv | awk -F"\t" '{n=split($1,a,","); c=0; for (i = 1; i <= n; i++) {if(a[i] ~ /^OMIM/) {n=split(a[i],b,":"); if( c == 1) {printf "," b[2]}; if( c == 0) {printf b[2]}; c=1}  } if( c == 0) {printf "NA"} print ""}' > $here/data-$batch/$gene.temp4.tsv
cat $here/data-$batch/$gene.temp4.tsv | awk -F"\t" '{n=split($1,a,"|"); print a[1]}' > $here/data-$batch/$gene.OMIM.tsv
paste -d"\t" $here/data-$batch/$gene.temp1.tsv $here/data-$batch/$gene.REV.tsv $here/data-$batch/$gene.OMIM.tsv $here/data-$batch/$gene.ORI.tsv > $here/data-$batch/$gene-$batch-dbNFSP4.0.tsv
rm $here/data-$batch/$gene.temp*.tsv $here/data-$batch/$gene.REV.tsv $here/data-$batch/$gene.OMIM.tsv $here/data-$batch/$gene.ORI.tsv

