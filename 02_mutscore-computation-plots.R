
library(gridExtra)
library(MASS)
library(pROC)
library(dplyr)
library(randomForest)
require(caTools)
library(rpart)
library(rpart.plot)
library(stringr)

### TO BE CHANGED ACCORDINGLY
here="/home/mquinodo/SYNO/WES/Clinvar2/"
version="20201121"

###

dataLoc=paste("data-",version,"/",sep="")
plot=paste("plots-",version,"/",sep="")

###

data<-read.table(file=paste(here,dataLoc,"ALL.tsv",sep=""),header=FALSE,sep="\t")  # OK

save(data,file=paste(here,dataLoc,"ALL.RData",sep=""))

load(paste(here,dataLoc,"ALL.RData",sep=""))

# changing type of column 1 (PLP, BLB,...)
data[,1]=as.character(data[,1])

# setting gnomAD frequency as the maximum between gnomAD-exomes and gnomAD-genomes and 0.0000001 if 0
data[,56]=as.character(data[,56])
data[which(data[,56]=="."),56]=0.000001
data[,56]=as.numeric(data[,56])
data[which(data[,56]==0),56]=0.000001
data[,56]=log10(data[,56])
data[,57]=as.character(data[,57])
data[which(data[,57]=="."),57]=0.000001
data[,57]=as.numeric(data[,57])
data[which(data[,57]==0),57]=0.000001
data[,57]=log10(data[,57])
data[,56]=apply(data[,c(56,57)],1,max)

save(data, file = paste(here,dataLoc,"data1.RData",sep=""))

load(file = paste(here,dataLoc,"data1.RData",sep=""))

# transforming protein notation to position and store in column 63
data[,63]=0
temp=gsub("f","",gsub(".{1}$","",substring(gsub("p.","",data[,4]),2)))
temp2=temp
temp3=temp
for (i in 1:length(temp)){
  temp2[i]=strsplit(temp[i],"_")[[1]][1]
  temp3[i]=strsplit(temp2[i],"del")[[1]][1]
}
data[,63]=as.numeric(temp3)

save(data, file = paste(here,dataLoc,"data2.RData",sep=""))

load(file = paste(here,dataLoc,"data2.RData",sep=""))

# shortening of names in column 1
data[which(data[,1]==paste("clinvar-",version,"-PLP",sep="")),1]="PLP"
data[which(data[,1]==paste("clinvar-",version,"-BLB",sep="")),1]="BLB"
data[which(data[,1]==paste("clinvar-",version,"-VUS",sep="")),1]="VUS"
data[which(data[,1]==paste("clinvar-",version,"-CON",sep="")),1]="CON"

# removing TTN other isoforms because they are too many and all PLP are present in this one (NM_001267550.2)
data=data[-which(data[,6]=="TTN" & data[,2]!="NM_001267550.2"),]

# separating genes with at least one PLP (in data) and other genes (in others)
PLPgene=unique(data[which(data[,1]=="PLP"),6])
others=data[which(is.element(data[,6],PLPgene)==F),]
data=data[which(is.element(data[,6],PLPgene)),]

# put common gnomAD (AF> max of PLPs) as BLB and annotate max-AF-PLP by gene and number of PLP per gene
genelist=unique(data[,6])
data[,64]=0
data[,65]=0
others[,64]=0
others[,65]=0
for (i in 1:length(genelist)){
  print(i)
  PLP=which(data[,6]==genelist[i] & data[,1]=="PLP")
  ALL=which(data[,6]==genelist[i])
  maxPLP=max(data[PLP,56])
  data[ALL,64]=length(PLP)
  data[ALL,65]=maxPLP
  if(length(PLP)>9){
    data[ALL[which(data[ALL,1]=="gnomAD-all" & data[ALL,56]>maxPLP)],1]="OK"
  }
}
data[which(data[,1]=="OK"),57]=(-7) # common gnomAD-BLB have "-7" value in col 57
data[which(data[,1]=="OK"),1]="BLB"
data=data[-which(data[,1]=="gnomAD-all"),]

# removing missense affecting Met1 (misstart)
start=which(data[,63]==1)
if(length(start>0)){data=data[-start,]}

# taking only non-synonymous variants (missense)
data=data[which(data[,5]=="nonsynonymous SNV"),]
others=others[which(others[,5]=="nonsynonymous SNV"),]

save(data,others, file = paste(here,dataLoc,"data3.RData",sep=""))

load(file = paste(here,dataLoc,"data3.RData",sep=""))

# unique name with isoform in column 66 (e.g. isoform-chr-begin-end-ref-alt or "NM_198317.3-1-897009-897010-A-G")
data[,8]=gsub("chr","",data[,8])
data[,66]=""
data[,66]=gsub(" ","",apply(data[,c(2,8:12)],1,paste,collapse="-"))
others[,8]=gsub("chr","",others[,8])
others[,66]=""
others[,66]=gsub(" ","",apply(others[,c(2,8:12)],1,paste,collapse="-"))

# removing doubles between true BLBs and common gnomAD BLBs
BLB=which(data[,1]=="BLB")
DF <- as.matrix(data.frame(table(data[BLB,66])))
b=which(as.numeric(DF[,2])>1)
list=as.character(DF[b,1])
a=which(is.element(data[,66],list) & data[,57]==(-7) & data[,1]=="BLB")
if(length(a)>0) {data=data[-a,]}

# test if ok (length(b) should be 0)
BLB=which(data[,1]=="BLB")
DF <- as.matrix(data.frame(table(data[BLB,66])))
b=which(as.numeric(DF[,2])>1)
length(b)

save(data,others, file = paste(here,dataLoc,"data4.RData",sep=""))

load(file = paste(here,dataLoc,"data4.RData",sep=""))

# merging all genes
data=rbind(data,others)
rm(others)

# transform protein change notation to amino acid change (p.R154G -> RG)
temp=gsub("p.","",data[,4])
temp2=gsub('[[:digit:]]+', '', temp)
data[,14]=temp2

# separating training and other variants (passenger)
train=data[which(data[,1]=="PLP" | data[,1]=="BLB"),]
train[,1]=as.factor(train[,1])
passenger=data[which(data[,1]=="VUS" | data[,1]=="CON" | data[,1]=="dbNFSP4.0"),]
passenger[,1]=as.factor(passenger[,1])
rm(data)

# computing aa changes score from the training set
train2=train[,c(1,6,14,5)]
train2=unique(train2)
PLP=which(train2[,1]=="PLP")
BLB=which(train2[,1]=="BLB")
changes=unique(train2[,3])
freq<-matrix(nrow=length(changes),ncol=4)
freq[,1]=changes
for (i in 1:length(changes)){
  freq[i,2]=length(which(train2[PLP,3]==changes[i]))
  freq[i,3]=length(which(train2[BLB,3]==changes[i]))
}
freq[,2]=as.numeric(freq[,2])
freq[,3]=as.numeric(freq[,3])
freq[,4]=as.numeric(freq[,2])/(as.numeric(freq[,2])+as.numeric(freq[,3]))

# annotating aa change score for all variants in col 69
train[,67]=0
passenger[,67]=0
train[,68]=0
passenger[,68]=0
train[,69]=0
passenger[,69]=0
for (i in 1:dim(freq)[1]){
  print(i)
  train[which(train[,14]==freq[i,1]),69]=freq[i,4]
  passenger[which(passenger[,14]==freq[i,1]),69]=freq[i,4]
}

train[,69]=as.numeric(as.character(train[,69]))
passenger[,69]=as.numeric(as.character(passenger[,69]))

save(train,passenger, file = paste(here,dataLoc,"data5.RData",sep=""))

load(file = paste(here,dataLoc,"data5.RData",sep=""))

####### positional score
train[,6]=as.character(train[,6])
train[,1]=as.factor(as.character(train[,1]))
passenger[,6]=as.character(passenger[,6])
passenger[,1]=as.factor(as.character(passenger[,1]))

PLPgene=unique(train[which(train[,1]=="PLP"),6])
others=passenger[which(is.element(passenger[,6],PLPgene)==F),]
passenger2=passenger[which(is.element(passenger[,6],PLPgene)==T),]
rm(passenger)

genelist=unique(c(as.character(train[,2]),as.character(passenger2[,2])))

# 71 = nb of PLP for isoform; 70 = positional score; 72 = "LOW" if less than 10 PLPs or less than 1 BLB
train[,70]=0
train[,71]=0
train[,72]=0
passenger2[,70]=0
passenger2[,71]=0
passenger2[,72]=0
others[,70]=0
others[,71]=0
others[,72]=0
c=0
for (i in genelist){
  c=c+1
  print(c)
  seltrain=which(train[,2]==i)
  selpass=which(passenger2[,2]==i)
  PLP=which(train[,2]==i & train[,1]=="PLP")
  BLB=which(train[,2]==i & train[,1]=="BLB")
  train[seltrain,71]=length(PLP)
  passenger2[selpass,71]=length(PLP)
  if(length(PLP)>9 & length(BLB)>0){
    rf <- randomForest(
      V1 ~ V63,
      data=train[seltrain,], ntree = 500, importance = TRUE
    )
    train[seltrain,70]=rf$votes[,2]
    if(length(selpass)>0){pred=predict(rf, newdata=passenger2[selpass,],"prob")
    passenger2[selpass,70]=pred[,2]}
  } else {
    train[seltrain,72]="LOW"
    if(length(selpass)>0){passenger2[selpass,72]="LOW"}
  }
}

# putting positional score at 0 for isoforms with not sufficent nb of PLPs
i="LOW"
seltrain=which(train[,72]==i)
train[seltrain,70]=0
selpass=which(passenger2[,72]==i)
passenger2[selpass,70]=0

passenger=passenger2
rm(passenger2)

save(train,passenger,others, file = paste(here,dataLoc,"data6.RData",sep=""))

load(file = paste(here,dataLoc,"data6.RData",sep=""))

length(unique(c(train[,2],passenger[,2],others[,2])))

######################

# finding regions with enriched PLP variants

datareg=unique(passenger[which(passenger[,1]=="dbNFSP4.0"),c(2,63,70,8,9)]) 

save(train,datareg, file = paste(here,dataLoc,"data6-forregions.RData",sep=""))

load(file = paste(here,dataLoc,"data6-forregions.RData",sep=""))

iso=unique(train[which(train[,1]=="PLP"),2])
train=train[which(is.element(train[,2],iso)),]
datareg=datareg[which(is.element(datareg[,1],iso)),]

limh=0.05

allregions=c()
for (j in 1:length(iso)){
  print(j)
  # taking all changes for one isoform with: prot-pos / pos-score / chr / begin
  gene=datareg[which(datareg[,1]==iso[j]),c(2,3,4,5)]
  train2=train[which(train[,2]==iso[j]),]
  gene=gene[sort(gene[,1],index.return=T)$ix,] # sort by aa pos
  if (dim(gene)[1]>0){
	  regions=c()
	  inside=F
	  count=0
	  for (i in 1:dim(gene)[1]){
	    # if not in region and threshold is achivied
	    if(inside==F & gene[i,2]>=limh){
	      count=count+1
	      # begin region with: gene / isoform / prot-pos / prot-pos / chr / begin / begin
	      regions=rbind(regions,c(as.character(train2[1,6]),iso[j],gene[i,1],gene[i,1],gene[i,3],gene[i,4],gene[i,4]))
	      inside=T
	    }
	    # if in region but we are leaving it
	    if(inside==T & gene[i,2]<limh){
	      # change prot-pos to end and change begin pos to end pos
	      regions[count,4]=gene[i,1]
	      regions[count,7]=gene[i,4]
	      inside=F
	    }
	    # if we are inside region but reach the end of the gene
	    if(inside==T & i==dim(gene)[1]){
	      regions[count,4]=gene[i,1]
	      regions[count,7]=gene[i,4]
	      inside=F
	    }
	  }
	  if(count!=0){
	    regions=cbind(regions,regions)
	    for (i in 1:dim(regions)[1]){
	      # 8=nb-PLP-inside / 9=nb-BLB-inside / 10=perc-PLP-inside / 11=perc-region-size-overGene / 12=nb-tot-PLP-gene / 13=max-prot-pos / 14=region-size(aa)
	      regions[i,8]=length(which(train2[,63]>=as.numeric(regions[i,3]) & train2[,63]<=as.numeric(regions[i,4]) & train2[,1]=="PLP"))
	      regions[i,9]=length(which(train2[,63]>=as.numeric(regions[i,3]) & train2[,63]<=as.numeric(regions[i,4]) & train2[,1]=="BLB"))
	      regions[i,10]=as.numeric(regions[i,8])/length(which(train2[,1]=="PLP"))
	      regions[i,11]=(as.numeric(regions[i,4])-as.numeric(regions[i,3]))/max(gene[,1])
	      regions[i,12]=length(which(train2[,1]=="PLP"))
	      regions[i,13]=max(gene[,1])
	      regions[i,14]=abs(as.numeric(regions[i,4])-as.numeric(regions[i,3]))
	    }
	    allregions=rbind(allregions,regions[,1:14])
	  }
	}
}
allregions1=allregions

# same for BLB regions

liml=0.05

allregions=c()
for (j in 1:length(iso)){
  print(j)
  gene=datareg[which(datareg[,1]==iso[j]),c(2,3,4,5)]
  train2=train[which(train[,2]==iso[j]),]
  gene=gene[sort(gene[,1],index.return=T)$ix,]
  if (dim(gene)[1]>0){
	  regions=c()
	  inside=F
	  count=0
	  for (i in 1:dim(gene)[1]){
	    if(inside==F & gene[i,2]<=liml){
	      count=count+1
	      regions=rbind(regions,c(as.character(train2[1,6]),iso[j],gene[i,1],gene[i,1],gene[i,3],gene[i,4],gene[i,4]))
	      inside=T
	    }
	    if(inside==T & gene[i,2]>liml){
	      regions[count,4]=gene[i,1]
	      regions[count,7]=gene[i,4]
	      inside=F
	    }
	    if(inside==T & i==dim(gene)[1]){
	      regions[count,4]=gene[i,1]
	      regions[count,7]=gene[i,4]
	      inside=F
	    }
	  }
	  if(count!=0){
	  	regions=cbind(regions,regions)
	    for (i in 1:dim(regions)[1]){
	      regions[i,8]=length(which(train2[,63]>=as.numeric(regions[i,3]) & train2[,63]<=as.numeric(regions[i,4]) & train2[,1]=="PLP"))
	      regions[i,9]=length(which(train2[,63]>=as.numeric(regions[i,3]) & train2[,63]<=as.numeric(regions[i,4]) & train2[,1]=="BLB"))
	      regions[i,10]=as.numeric(regions[i,9])/length(which(train2[,1]=="BLB"))
	      regions[i,11]=(as.numeric(regions[i,4])-as.numeric(regions[i,3]))/max(gene[,1])
	      regions[i,12]=length(which(train2[,1]=="BLB"))
	      regions[i,13]=max(gene[,1])
	      regions[i,14]=abs(as.numeric(regions[i,4])-as.numeric(regions[i,3]))
	    }
	    allregions=rbind(allregions,regions[,1:14])
	  }
	}
}
allregions2=allregions

#### cluster PLP

# taking only regions with more than 5 PLPs
regionsp=allregions1[which(as.numeric(allregions1[,8])>=5),]

# regionsp
# 1=gene / 2=isoform / 3=begin-region / 4=end-region / 5=chr / 6=begin-pos / 7=end-pos / 8=nb-PLP / 9=nb-BLB
# 10=%-PLP-inside / 11=%-region-gene / 12=nb-PLP-tot / 13=size-gene / 14=size-region

iso=unique(regionsp[,2])
resi=cbind(iso,iso,iso,iso,iso,iso,iso,iso)
for (i in 1:length(iso)){
  print(i)
  a=which(regionsp[,2]==iso[i])
  # gene name
  resi[i,2]=regionsp[a[1],1]
  # perc-PLP inside regionS
  resi[i,3]=sum(as.numeric(regionsp[a,10]))
  # perc gene inside regionS
  resi[i,4]=sum(as.numeric(regionsp[a,11]))
  # minimum score
  resi[i,5]=min(as.numeric(resi[i,3]),1-as.numeric(resi[i,4]))
  # OR
  resi[i,6]=fisher.test(rbind(c(sum(as.numeric(regionsp[a,8])),as.numeric(regionsp[a[1],12])-sum(as.numeric(regionsp[a,8]))),c(sum(as.numeric(regionsp[a,14])),as.numeric(regionsp[a[1],13])-sum(as.numeric(regionsp[a,14])))))$estimate
  # p-value
  resi[i,7]=fisher.test(rbind(c(sum(as.numeric(regionsp[a,8])),as.numeric(regionsp[a[1],12])-sum(as.numeric(regionsp[a,8]))),c(sum(as.numeric(regionsp[a,14])),as.numeric(regionsp[a[1],13])-sum(as.numeric(regionsp[a,14])))))$p.value
}
resi[,8]=p.adjust(as.numeric(resi[,7]),method="fdr")
resi1=resi

# resi
# 1=isoform / 2=gene / 3=%PLP-regionS / 4=%gene-regionS / 5=score / 6=Fisher-estimate / 7=Fisher-p.value

#### cluster BLB

# taking only regions with more than 10 BLBs
regionsp=allregions2[which(as.numeric(allregions2[,9])>=10),]

iso=unique(regionsp[,2])
resi=cbind(iso,iso,iso,iso,iso,iso,iso,iso)
for (i in 1:length(iso)){
  print(i)
  a=which(regionsp[,2]==iso[i])
  resi[i,2]=regionsp[a[1],1]
  resi[i,3]=sum(as.numeric(regionsp[a,10]))
  resi[i,4]=sum(as.numeric(regionsp[a,11]))
  resi[i,5]=min(as.numeric(resi[i,3]),1-as.numeric(resi[i,4]))
  resi[i,6]=fisher.test(rbind(c(sum(as.numeric(regionsp[a,9])),as.numeric(regionsp[a[1],12])-sum(as.numeric(regionsp[a,9]))),c(sum(as.numeric(regionsp[a,14])),as.numeric(regionsp[a[1],13])-sum(as.numeric(regionsp[a,14])))))$estimate
  resi[i,7]=fisher.test(rbind(c(sum(as.numeric(regionsp[a,9])),as.numeric(regionsp[a[1],12])-sum(as.numeric(regionsp[a,9]))),c(sum(as.numeric(regionsp[a,14])),as.numeric(regionsp[a[1],13])-sum(as.numeric(regionsp[a,14])))))$p.value
}
resi[,8]=p.adjust(as.numeric(resi[,7]),method="fdr")
resi2=resi

regionsp=allregions1[which(as.numeric(allregions1[,8])>=5),]
regionsb=allregions2[which(as.numeric(allregions2[,9])>=10),]

colnames(resi1)=c("isoform","gene","%PLPin","%genein","score","OR","p-value","cor-p")
colnames(resi2)=c("isoform","gene","%PLPin","%genein","score","OR","p-value","cor-p")

write.table(resi1,file=paste(here,dataLoc,"regionsPLP.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(resi2,file=paste(here,dataLoc,"regionsBLB.tsv",sep=""),quote=F,sep="\t",row.names=F)

save(regionsb,regionsp,allregions1,allregions2,resi1,resi2, file = paste(here,dataLoc,"regions.RData",sep=""))
save(regionsb,regionsp,allregions1,allregions2,resi1,resi2, file = paste(here,"Shiny/regions.RData",sep=""))



####### positional score for permutation test p-values

train[,6]=as.character(train[,6])
train[,1]=as.factor(as.character(train[,1]))
passenger[,6]=as.character(passenger[,6])
passenger[,1]=as.factor(as.character(passenger[,1]))

PLPgene=unique(train[which(train[,1]=="PLP"),6])
passenger2=passenger[which(is.element(passenger[,6],PLPgene)==T),]
rm(passenger)

# selecting only transcripts with more than 10 PLPs and one BLB
PLP=table(train[which(train[,1]=="PLP"),2])
s1=names(PLP[which(as.numeric(PLP)>=10)])
BLB=table(train[which(train[,1]=="BLB"),2])
s2=names(BLB[which(as.numeric(BLB)>=0)])
s3=intersect(s1,s2)

train=train[which(is.element(train[,2],s3)==T),]
passenger2=passenger2[which(is.element(passenger2[,2],s3)==T),]

passenger2[,16:55]=0

save(train,passenger2, file = paste(here,dataLoc,"data5-review2.RData",sep=""))

###

load(file = paste(here,dataLoc,"data5-review2.RData",sep=""))

# 71 = nb of PLP for isoform; 70 = positional score; 72 = "LOW" if less than 10 PLPs or less than 1 BLB
train[,70]=0
train[,71]=0
train[,72]=0
passenger2[,70]=0
passenger2[,71]=0
passenger2[,72]=0

trainALL=train
passALL=passenger2

genelist=unique(c(as.character(trainALL[,2]),as.character(passALL[,2]))) # length=2578
pres1=matrix(nrow=length(genelist),ncol=2)
pres2=matrix(nrow=length(genelist),ncol=1000)
pres1[,1]=genelist

for (i in 1:length(genelist)){

	{

		train=trainALL[which(trainALL[,2]==genelist[i]),]
		passenger2=passALL[which(passALL[,2]==genelist[i]),]

	  print(i)
	  PLP=which(train[,1]=="PLP")
	  BLB=which(train[,1]=="BLB")
	  train[,71]=length(PLP)
	  print("step1 finished")
	  rf <- randomForest(
	    V1 ~ V63,
	    data=train, ntree = 500, importance = TRUE
	  )
	  train[,70]=rf$votes[,2]
	  pred=predict(rf, newdata=passenger2,"prob")
	  passenger2[,70]=pred[,2]

	  datareg=unique(passenger2[which(passenger2[,1]=="dbNFSP4.0"),c(2,63,70,8,9)]) 

	  limh=0.05
	  allregions=c()
	  # taking all changes for one isoform with: prot-pos / pos-score / chr / begin
	  gene=datareg[,c(2,3,4,5)]
	  train2=train
	  gene=gene[sort(gene[,1],index.return=T)$ix,] # sort by aa pos
	  if (dim(gene)[1]>0){
		  regions=c()
		  inside=F
		  count=0
		  for (k in 1:dim(gene)[1]){
		    # if not in region and threshold is achivied
		    if(inside==F & gene[k,2]>=limh){
		      count=count+1
		      # begin region with: gene / isoform / prot-pos / prot-pos / chr / begin / begin
		      regions=rbind(regions,c(as.character(train2[1,6]),genelist[i],gene[k,1],gene[k,1],gene[k,3],gene[k,4],gene[k,4]))
		      inside=T
		    }
		    # if in region but we are leaving it
		    if(inside==T & gene[k,2]<limh){
		      # change prot-pos to end and change begin pos to end pos
		      regions[count,4]=gene[k,1]
		      regions[count,7]=gene[k,4]
		      inside=F
		    }
		    # if we are inside region but reach the end of the gene
		    if(inside==T & k==dim(gene)[1]){
		      regions[count,4]=gene[k,1]
		      regions[count,7]=gene[k,4]
		      inside=F
		    }
		  }
		  if(count!=0){
		    regions=cbind(regions,regions)
		    for (k in 1:dim(regions)[1]){
		      # 8=nb-PLP-inside / 9=nb-BLB-inside / 10=perc-PLP-inside / 11=perc-region-size-overGene / 12=nb-tot-PLP-gene / 13=max-prot-pos / 14=region-size(aa)
		      regions[k,8]=length(which(train2[,63]>=as.numeric(regions[k,3]) & train2[,63]<=as.numeric(regions[k,4]) & train2[,1]=="PLP"))
		      regions[k,9]=length(which(train2[,63]>=as.numeric(regions[k,3]) & train2[,63]<=as.numeric(regions[k,4]) & train2[,1]=="BLB"))
		      regions[k,10]=as.numeric(regions[k,8])/length(which(train2[,1]=="PLP"))
		      regions[k,11]=(as.numeric(regions[k,4])-as.numeric(regions[k,3]))/max(gene[,1])
		      regions[k,12]=length(which(train2[,1]=="PLP"))
		      regions[k,13]=max(gene[,1])
		      regions[k,14]=abs(as.numeric(regions[k,4])-as.numeric(regions[k,3]))
		    }
		    allregions=rbind(allregions,regions[,1:14])
		  }
		}
		regionsp=allregions[which(as.numeric(allregions[,8])>=5),]

		#pres1[i,2]=sum(as.numeric(regionsp[,14]))
		if(length(regionsp)>14){
			t1=sum(as.numeric(regionsp[,10]))
			t2=sum(as.numeric(regionsp[,11]))
			pres1[i,2]=min(t1,1-t2)
		} else if(length(regionsp)==14){
			t1=sum(as.numeric(regionsp[10]))
			t2=sum(as.numeric(regionsp[11]))
			#pres1[i,2]=min(t1,1-t2)
			pres1[i,2]=t1+1-t2
		} else {
			pres1[i,2]=0
		}

		trainALL2=train
		passALL2=passenger2

		if(as.numeric(pres1[i,2])>0){
			for (m in 1:1000){
				#print(m)

				train=trainALL2
				passenger2=passALL2

				maximum=max(passenger2[,63])

				train[which(train[,1]=="PLP"),63]=sample(1:maximum,length(which(train[,1]=="PLP")),replace=T)
				train[which(train[,1]=="BLB"),63]=sample(1:maximum,length(which(train[,1]=="BLB")),replace=T)
			  train[,71]=length(PLP)

			  rf <- randomForest(
			    V1 ~ V63,
			    data=train, ntree = 500
			  )
			  train[,70]=rf$votes[,2]
			  pred=predict(rf, newdata=passenger2,"prob")
			  passenger2[,70]=pred[,2]
				
			  datareg=unique(passenger2[which(passenger2[,1]=="dbNFSP4.0"),c(2,63,70,8,9)]) 


			  limh=0.05
			  # taking all changes for one isoform with: prot-pos / pos-score / chr / begin
			  gene=datareg[,c(2,3,4,5)]
			  train2=train
			  gene=gene[sort(gene[,1],index.return=T)$ix,] # sort by aa pos
			  if (dim(gene)[1]>0){
				  regions=c()
				  inside=F
				  count=0
				  for (k in 1:dim(gene)[1]){
				    # if not in region and threshold is achivied
				    if(inside==F & gene[k,2]>=limh){
				      count=count+1
				      # begin region with: gene / isoform / prot-pos / prot-pos / chr / begin / begin
				      regions=rbind(regions,c(as.character(train2[1,6]),genelist[i],gene[k,1],gene[k,1],gene[k,3],gene[k,4],gene[k,4]))
				      inside=T
				    }
				    # if in region but we are leaving it
				    if(inside==T & gene[k,2]<limh){
				      # change prot-pos to end and change begin pos to end pos
				      regions[count,4]=gene[k,1]
				      regions[count,7]=gene[k,4]
				      inside=F
				    }
				    # if we are inside region but reach the end of the gene
				    if(inside==T & k==dim(gene)[1]){
				      regions[count,4]=gene[k,1]
				      regions[count,7]=gene[k,4]
				      inside=F
				    }
				  }
				  if(count!=0){
				    regions=cbind(regions,regions)
				    for (k in 1:dim(regions)[1]){
				      # 8=nb-PLP-inside / 9=nb-BLB-inside / 10=perc-PLP-inside / 11=perc-region-size-overGene / 12=nb-tot-PLP-gene / 13=max-prot-pos / 14=region-size(aa)
				      regions[k,8]=length(which(train2[,63]>=as.numeric(regions[k,3]) & train2[,63]<=as.numeric(regions[k,4]) & train2[,1]=="PLP"))
				      # regions[k,9]=length(which(train2[,63]>=as.numeric(regions[k,3]) & train2[,63]<=as.numeric(regions[k,4]) & train2[,1]=="BLB"))
				       regions[k,10]=as.numeric(regions[k,8])/length(which(train2[,1]=="PLP"))
				       regions[k,11]=(as.numeric(regions[k,4])-as.numeric(regions[k,3]))/max(gene[,1])
				      # regions[k,12]=length(which(train2[,1]=="PLP"))
				      # regions[k,13]=max(gene[,1])
				      regions[k,14]=abs(as.numeric(regions[k,4])-as.numeric(regions[k,3]))
				    }
				  }
				}
				regionsp=regions[which(as.numeric(regions[,8])>=5),]


				# regionsp
				# 1=gene / 2=isoform / 3=begin-region / 4=end-region / 5=chr / 6=begin-pos / 7=end-pos / 8=nb-PLP / 9=nb-BLB
				# 10=%-PLP-inside / 11=%-region-gene / 12=nb-PLP-tot / 13=size-gene / 14=size-region
				#pres1[i,2]=sum(as.numeric(regionsp[,14]))
				if(length(regionsp)>14){
					t1=sum(as.numeric(regionsp[,10]))
					t2=sum(as.numeric(regionsp[,11]))
					#pres2[i,m]=min(t1,1-t2)
					pres2[i,m]=t1+1-t2

				} else if(length(regionsp)==14){
					t1=sum(as.numeric(regionsp[10]))
					t2=sum(as.numeric(regionsp[11]))
					#pres2[i,m]=min(t1,1-t2)
					pres2[i,m]=t1+1-t2

				} else {
					pres2[i,m]=0
				}

			}
		}

	}

	save(pres1,pres2,file=paste(here,dataLoc,"data5-scores.RData",sep=""))

}

########
### BENIGN

load(file = paste(here,dataLoc,"data5-review2.RData",sep=""))

# 71 = nb of PLP for isoform; 70 = positional score; 72 = "LOW" if less than 10 PLPs or less than 1 BLB
train[,70]=0
train[,71]=0
train[,72]=0
passenger2[,70]=0
passenger2[,71]=0
passenger2[,72]=0

trainALL=train
passALL=passenger2

genelist=unique(c(as.character(trainALL[,2]),as.character(passALL[,2]))) # length=2578
pres1=matrix(nrow=length(genelist),ncol=2)
pres2=matrix(nrow=length(genelist),ncol=1000)
pres1[,1]=genelist

for (i in 1:length(genelist)){

	{

		train=trainALL[which(trainALL[,2]==genelist[i]),]
		passenger2=passALL[which(passALL[,2]==genelist[i]),]

	  print(i)
	  PLP=which(train[,1]=="PLP")
	  BLB=which(train[,1]=="BLB")
	  train[,71]=length(PLP)
	  print("step1 finished")
	  rf <- randomForest(
	    V1 ~ V63,
	    data=train, ntree = 500, importance = TRUE
	  )
	  train[,70]=rf$votes[,2]
	  pred=predict(rf, newdata=passenger2,"prob")
	  passenger2[,70]=pred[,2]

	  datareg=unique(passenger2[which(passenger2[,1]=="dbNFSP4.0"),c(2,63,70,8,9)]) 

	  liml=0.05
	  allregions=c()
	  # taking all changes for one isoform with: prot-pos / pos-score / chr / begin
	  gene=datareg[,c(2,3,4,5)]
	  train2=train
	  gene=gene[sort(gene[,1],index.return=T)$ix,] # sort by aa pos
	  if (dim(gene)[1]>0){
		  regions=c()
		  inside=F
		  count=0
		  for (k in 1:dim(gene)[1]){
		    # if not in region and threshold is achivied
		    if(inside==F & gene[k,2]<=liml){
		      count=count+1
		      # begin region with: gene / isoform / prot-pos / prot-pos / chr / begin / begin
		      regions=rbind(regions,c(as.character(train2[1,6]),genelist[i],gene[k,1],gene[k,1],gene[k,3],gene[k,4],gene[k,4]))
		      inside=T
		    }
		    # if in region but we are leaving it
		    if(inside==T & gene[k,2]>liml){
		      # change prot-pos to end and change begin pos to end pos
		      regions[count,4]=gene[k,1]
		      regions[count,7]=gene[k,4]
		      inside=F
		    }
		    # if we are inside region but reach the end of the gene
		    if(inside==T & k==dim(gene)[1]){
		      regions[count,4]=gene[k,1]
		      regions[count,7]=gene[k,4]
		      inside=F
		    }
		  }
		  if(count!=0){
		    regions=cbind(regions,regions)
		    for (k in 1:dim(regions)[1]){
		      # 8=nb-PLP-inside / 9=nb-BLB-inside / 10=perc-PLP-inside / 11=perc-region-size-overGene / 12=nb-tot-PLP-gene / 13=max-prot-pos / 14=region-size(aa)
		      regions[k,8]=length(which(train2[,63]>=as.numeric(regions[k,3]) & train2[,63]<=as.numeric(regions[k,4]) & train2[,1]=="PLP"))
		      regions[k,9]=length(which(train2[,63]>=as.numeric(regions[k,3]) & train2[,63]<=as.numeric(regions[k,4]) & train2[,1]=="BLB"))
		      regions[k,10]=as.numeric(regions[k,9])/length(which(train2[,1]=="BLB"))
		      regions[k,11]=(as.numeric(regions[k,4])-as.numeric(regions[k,3]))/max(gene[,1])
		      regions[k,12]=length(which(train2[,1]=="BLB"))
		      regions[k,13]=max(gene[,1])
		      regions[k,14]=abs(as.numeric(regions[k,4])-as.numeric(regions[k,3]))
		    }
		    allregions=rbind(allregions,regions[,1:14])
		  }
		}
		regionsp=allregions[which(as.numeric(allregions[,9])>=10),]

		if(length(regionsp)>14){
			t1=sum(as.numeric(regionsp[,10]))
			t2=sum(as.numeric(regionsp[,11]))
			pres1[i,2]=min(t1,1-t2)
		} else if(length(regionsp)==14){
			t1=sum(as.numeric(regionsp[10]))
			t2=sum(as.numeric(regionsp[11]))
			pres1[i,2]=min(t1,1-t2)

		} else {
			pres1[i,2]=0
		}

		trainALL2=train
		passALL2=passenger2

		if(as.numeric(pres1[i,2])>0){
			for (m in 1:1000){
				#print(m)

				train=trainALL2
				passenger2=passALL2

				maximum=max(passenger2[,63])

				train[which(train[,1]=="PLP"),63]=sample(1:maximum,length(which(train[,1]=="PLP")),replace=T)
				train[which(train[,1]=="BLB"),63]=sample(1:maximum,length(which(train[,1]=="BLB")),replace=T)
			  train[,71]=length(PLP)

			  rf <- randomForest(
			    V1 ~ V63,
			    data=train, ntree = 500
			  )
			  train[,70]=rf$votes[,2]
			  pred=predict(rf, newdata=passenger2,"prob")
			  passenger2[,70]=pred[,2]
				
			  datareg=unique(passenger2[which(passenger2[,1]=="dbNFSP4.0"),c(2,63,70,8,9)]) 
					  
			  liml=0.05
			  allregions=c()
			  # taking all changes for one isoform with: prot-pos / pos-score / chr / begin
			  gene=datareg[,c(2,3,4,5)]
			  train2=train
			  gene=gene[sort(gene[,1],index.return=T)$ix,] # sort by aa pos
			  if (dim(gene)[1]>0){
				  regions=c()
				  inside=F
				  count=0
				  for (k in 1:dim(gene)[1]){
				    # if not in region and threshold is achivied
				    if(inside==F & gene[k,2]<=liml){
				      count=count+1
				      # begin region with: gene / isoform / prot-pos / prot-pos / chr / begin / begin
				      regions=rbind(regions,c(as.character(train2[1,6]),genelist[i],gene[k,1],gene[k,1],gene[k,3],gene[k,4],gene[k,4]))
				      inside=T
				    }
				    # if in region but we are leaving it
				    if(inside==T & gene[k,2]>liml){
				      # change prot-pos to end and change begin pos to end pos
				      regions[count,4]=gene[k,1]
				      regions[count,7]=gene[k,4]
				      inside=F
				    }
				    # if we are inside region but reach the end of the gene
				    if(inside==T & k==dim(gene)[1]){
				      regions[count,4]=gene[k,1]
				      regions[count,7]=gene[k,4]
				      inside=F
				    }
				  }
				  if(count!=0){
				    regions=cbind(regions,regions)
				    for (k in 1:dim(regions)[1]){
				      # 8=nb-PLP-inside / 9=nb-BLB-inside / 10=perc-PLP-inside / 11=perc-region-size-overGene / 12=nb-tot-PLP-gene / 13=max-prot-pos / 14=region-size(aa)
				      regions[k,8]=length(which(train2[,63]>=as.numeric(regions[k,3]) & train2[,63]<=as.numeric(regions[k,4]) & train2[,1]=="PLP"))
				      regions[k,9]=length(which(train2[,63]>=as.numeric(regions[k,3]) & train2[,63]<=as.numeric(regions[k,4]) & train2[,1]=="BLB"))
				      regions[k,10]=as.numeric(regions[k,9])/length(which(train2[,1]=="BLB"))
				      regions[k,11]=(as.numeric(regions[k,4])-as.numeric(regions[k,3]))/max(gene[,1])
				      regions[k,12]=length(which(train2[,1]=="BLB"))
				      regions[k,13]=max(gene[,1])
				      regions[k,14]=abs(as.numeric(regions[k,4])-as.numeric(regions[k,3]))
				    }
				    allregions=rbind(allregions,regions[,1:14])
				  }
				}
				regionsp=allregions[which(as.numeric(allregions[,9])>=10),]


				# regionsp
				# 1=gene / 2=isoform / 3=begin-region / 4=end-region / 5=chr / 6=begin-pos / 7=end-pos / 8=nb-PLP / 9=nb-BLB
				# 10=%-PLP-inside / 11=%-region-gene / 12=nb-PLP-tot / 13=size-gene / 14=size-region
				#pres1[i,2]=sum(as.numeric(regionsp[,14]))
				if(length(regionsp)>14){
					t1=sum(as.numeric(regionsp[,10]))
					t2=sum(as.numeric(regionsp[,11]))
					pres2[i,m]=min(t1,1-t2)
				} else if(length(regionsp)==14){
					t1=sum(as.numeric(regionsp[10]))
					t2=sum(as.numeric(regionsp[11]))
					pres2[i,m]=min(t1,1-t2)
				} else {
					pres2[i,m]=0
				}

			}
		}

	}

	save(pres1,pres2,file=paste(here,dataLoc,"data5-scores-BLB.RData",sep=""))

}


load(file=paste(here,dataLoc,"data5-scores.RData",sep=""))

sel=which(is.na(pres1[,2])==F)
pres1=pres1[sel,]
pres2=pres2[sel,]

a1=as.numeric(pres1[,2])
a2=apply(pres2,1,max)
p=1:length(a1)
for (i in 1:length(a1)){
	t=c(a1[i],pres2[i,])
	p[i]=(1001-rank(t))/1000
}
p1=p.adjust(p,method="fdr",n=length(p))

load(file = paste(here,dataLoc,"regions.RData",sep=""))
resi1[,3]=NA
for (i in 1:length(p)){
	resi1[which(resi1[,1]==pres1[i,1]),3]=p[i]
	resi1[which(resi1[,1]==pres1[i,1]),4]=p1[i]
}

newPLP=resi1[,c(3,4)]

########


load(file=paste(here,dataLoc,"data5-scores-BLB.RData",sep=""))

sel=which(is.na(pres1[,2])==F)
pres1=pres1[sel,]
pres2=pres2[sel,]

a1=as.numeric(pres1[,2])
a2=apply(pres2,1,max)
p=1:length(a1)
for (i in 1:length(a1)){
	t=c(a1[i],pres2[i,])
	p[i]=(1001-rank(t))/1000
}
p1=p.adjust(p,method="fdr",n=length(p))


load(file = paste(here,dataLoc,"regions.RData",sep=""))
resi2[,3]=NA
for (i in 1:length(p)){
	resi2[which(resi2[,1]==pres1[i,1]),3]=p[i]
	resi2[which(resi2[,1]==pres1[i,1]),4]=p1[i]
}

newBLB=resi2[,c(3,4)]

save(newPLP,newBLB,file=paste(here,dataLoc,"regions-correction.RData",sep=""))


# review correction

load(file = paste(here,dataLoc,"regions.RData",sep=""))
load(file = paste(here,dataLoc,"regions-correction.RData",sep=""))

resi1[,c(7,8)]=newPLP
resi2[,c(7,8)]=newBLB

save(regionsb,regionsp,allregions1,allregions2,resi1,resi2, file = paste(here,dataLoc,"regions-review.RData",sep=""))
save(regionsb,regionsp,allregions1,allregions2,resi1,resi2, file = paste(here,"Shiny/regions-review.RData",sep=""))

write.table(resi1,file=paste(here,dataLoc,"regionsPLP-review.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(resi2,file=paste(here,dataLoc,"regionsBLB-review.tsv",sep=""),quote=F,sep="\t",row.names=F)


########################################

load(file = paste(here,dataLoc,"data6.RData",sep=""))
#train passenger others

gnomad=train[which(train[,57]==(-7)),]
train=train[-which(train[,57]==(-7)),]
gnomad[,1]="gnomAD-all"
passenger=rbind(passenger,gnomad)

# adding unique change in 73 (e.g. "1:   897009:   897010:A:G")
train[,73]=0
train[,73]=apply(train[,8:12],1,paste,collapse=":")
passenger[,73]=0
passenger[,73]=apply(passenger[,8:12],1,paste,collapse=":")
others[,73]=0
others[,73]=apply(others[,8:12],1,paste,collapse=":")

save(train, file = paste(here,dataLoc,"data7-train.RData",sep=""))
save(passenger, file = paste(here,dataLoc,"data7-passenger.RData",sep=""))
save(others, file = paste(here,dataLoc,"data7-others.RData",sep=""))

###############################

# this part is to select the isoform for variants present in multiple isoforms (the one affecting the isoform with the most PLPs)

load(file = paste(here,dataLoc,"data7-train.RData",sep=""))

# table with isoforms and count of PLPs
DF <- as.matrix(data.frame(table(train[,2])))
DF[,2]=as.numeric(DF[,2])+rnorm(n=length(DF[,2]),mean=0.5,sd=0.01)

genelist=unique(train[,6])
taken=as.numeric(train[,1])
taken[1:length(taken)]=1
# loop on genes
for (i in 1:length(genelist)){
  print(i)
  a=which(train[,6]==genelist[i])
  # list of isoforms from this gene
  isof=DF[which(is.element(DF[,1],train[a,2])),]
  # if more than 1 isoform
  if(length(unique(train[a,2]))>1){
    posi=unique(train[a,73])
    # loop on unique changes
    for (k in 1:length(posi)){
      a2=which(train[a,73]==posi[k])
      if(length(a2)>1){
        temp=0
        sele=0
        score=0
        for (j in 1:length(a2)){
          score=as.numeric(DF[which(DF[,1]==train[a[a2[j]],2]),2])
          if(score>temp){
            temp=score
            sele=j
          }
        }
        taken[a[a2[-sele]]]=0
      }
    }
  }
}

train=train[which(taken==1),]

save(train,DF, file = paste(here,dataLoc,"data8-train.RData",sep=""))

##############

load(file = paste(here,dataLoc,"data8-train.RData",sep=""))
#load(file = paste(here,dataLoc,"data7-passenger.RData",sep=""))

passenger1=passenger[which(passenger[,1]=="CON"),]
passenger2=passenger[which(passenger[,1]=="VUS"),]
passenger3=passenger[which(passenger[,1]=="gnomAD-all"),]
passenger4=passenger[which(passenger[,1]=="dbNFSP4.0"),]

rm(passenger)

genelist=unique(passenger1[,6])
taken=passenger1[,63]
taken[1:length(taken)]=1
for (i in 1:length(genelist)){
  print(i)
  a=which(passenger1[,6]==genelist[i])
  isof=DF[which(is.element(DF[,1],passenger1[a,2])),]
  if(length(unique(passenger1[a,2]))>1){
    posi=unique(passenger1[a,73])
    for (k in 1:length(posi)){
      a2=which(passenger1[a,73]==posi[k])
      if(length(a2)>1){
        temp=0
        sele=0
        score=0
        for (j in 1:length(a2)){
          if(is.element(passenger1[a[a2[j]],2],DF[,1])){
            score=as.numeric(DF[which(DF[,1]==passenger1[a[a2[j]],2]),2])
            if(score>temp){
              temp=score
              sele=j
            }
          }
        }
        taken[a[a2[-sele]]]=0
      }
    }
  }
}
passenger1=passenger1[which(taken==1),]

genelist=unique(passenger2[,6])
taken=passenger2[,63]
taken[1:length(taken)]=1
for (i in 1:length(genelist)){
  print(i)
  a=which(passenger2[,6]==genelist[i])
  isof=DF[which(is.element(DF[,1],passenger2[a,2])),]
  if(length(unique(passenger2[a,2]))>1){
    posi=unique(passenger2[a,73])
    for (k in 1:length(posi)){
      a2=which(passenger2[a,73]==posi[k])
      if(length(a2)>1){
        temp=0
        sele=0
        score=0
        for (j in 1:length(a2)){
          if(is.element(passenger2[a[a2[j]],2],DF[,1])){
            score=as.numeric(DF[which(DF[,1]==passenger2[a[a2[j]],2]),2])
            if(score>temp){
              temp=score
              sele=j
            }
          }
        }
        taken[a[a2[-sele]]]=0
      }
    }
  }
}
passenger2=passenger2[which(taken==1),]

genelist=unique(passenger3[,6])
taken=passenger3[,63]
taken[1:length(taken)]=1
for (i in 1:length(genelist)){
  print(i)
  a=which(passenger3[,6]==genelist[i])
  isof=DF[which(is.element(DF[,1],passenger3[a,2])),]
  if(length(unique(passenger3[a,2]))>1){
    posi=unique(passenger3[a,73])
    for (k in 1:length(posi)){
      a2=which(passenger3[a,73]==posi[k])
      if(length(a2)>1){
        temp=0
        sele=0
        score=0
        for (j in 1:length(a2)){
          if(is.element(passenger3[a[a2[j]],2],DF[,1])){
            score=as.numeric(DF[which(DF[,1]==passenger3[a[a2[j]],2]),2])
            if(score>temp){
              temp=score
              sele=j
            }
          }
        }
        taken[a[a2[-sele]]]=0
      }
    }
  }
}
passenger3=passenger3[which(taken==1),]

genelist=unique(passenger4[,6])
taken=passenger4[,63]
taken[1:length(taken)]=1
for (i in 1:length(genelist)){
  print(i)
  a=which(passenger4[,6]==genelist[i])
  isof=DF[which(is.element(DF[,1],passenger4[a,2])),]
  if(length(unique(passenger4[a,2]))>1){
    posi=unique(passenger4[a,73])
    for (k in 1:length(posi)){
      a2=which(passenger4[a,73]==posi[k])
      if(length(a2)>1){
        temp=0
        sele=0
        score=0
        for (j in 1:length(a2)){
          if(is.element(passenger4[a[a2[j]],2],DF[,1])){
            score=as.numeric(DF[which(DF[,1]==passenger4[a[a2[j]],2]),2])
            if(score>temp){
              temp=score
              sele=j
            }
          }
        }
        taken[a[a2[-sele]]]=0
      }
    }
  }
}
passenger4=passenger4[which(taken==1),]

passenger=rbind(passenger1,passenger2,passenger3,passenger4)

save(passenger, file = paste(here,dataLoc,"data8-passenger.RData",sep=""))

rm(passenger,passenger1,passenger2,passenger3,passenger4)

######

load(file = paste(here,dataLoc,"data7-other.RData",sep=""))

a=duplicated(others[,73])

others=others[-which(a==T),]

save(others, file = paste(here,dataLoc,"data8-others.RData",sep=""))

#####################################################

load(file = paste(here,dataLoc,"data8-others.RData",sep=""))
load(file = paste(here,dataLoc,"data8-passenger.RData",sep=""))
load(file = paste(here,dataLoc,"data8-train.RData",sep=""))

data=rbind(train,passenger,others)
rm(train,passenger,others)

# putting median where missing values for all scores excpet splicing and pext (zeroed)

a=duplicated(data[,73])

### PREDICTORS: SIFT and others
for (m in c(16,19:55)){
  print(m)
  data[,m]=as.character(data[,m])
  med=median(as.numeric(data[which(data[,m]!="." & data[,m]!="-" & is.na(data[,m])==FALSE && a==F),m]))
  data[which(data[,m]=="." | data[,m]=="-" | is.na(data[,m])==TRUE),m]=med
  data[,m]=as.numeric(data[,m])
}
# splicing and pext
for (m in c(17,18,58,59)){
  print(m)
  data[,m]=as.character(data[,m])
  data[which(data[,m]=="." | data[,m]=="-" | is.na(data[,m])==TRUE),m]=0
  data[,m]=as.numeric(data[,m])
}

train=data[which(data[,1]=="PLP" | data[,1]=="BLB"),]
VUSCON=data[which(data[,1]=="CON" | data[,1]=="VUS"),]
gnomad=data[which(data[,1]=="gnomAD-all"),]
dbnfsp=data[which(data[,1]=="dbNFSP4.0"),]

save(train,VUSCON,gnomad,dbnfsp, file = paste(here,dataLoc,"data9.RData",sep=""))

load(file = paste(here,dataLoc,"data9.RData",sep=""))


############################################

                 # MutScore #

############################################

# removing BLB variants from the training set in genes that have no PLP variants
table(train[,1])
genes=unique(train[which(train[,1]=="PLP"),6])
train=train[which(is.element(train[,6],genes)),]
table(train[,1])

# manual extraction of genes with no Mendelian phenotypes and wrong variants

a=duplicated(train[,73])

prob=train[which(is.element(train[,73],train[a,73])),]

write.table(prob,file=paste(here,dataLoc,"problematicPLP.tsv",sep=""),quote=F,sep="\t",row.names=F)

# genes names written after curation of problematicPLP.tsv file
wrong=c("FPGT-TNNI3K","UGT1A5","UGT1A8","UGT1A6","UGT1A7","UGT1A3","UGT1A9","UGT1A4","UGT1A10","ABHD14A-ACY1","SMN2","KAAG1","CNPY3-GNMT","PGBD3","INS-IGF2","NDUFC2-KCTD14","FXYD6-FXYD2","BIVM-ERCC5","ST20-MTHFS","CORO7-PAM16","CORO7-PAM16","KCNE1B","CBSL","U2AF1L5","CRYAA2","SIK1B","LOC102724788","OPN1MW3","OPN1MW2")
train=train[-which(is.element(train[,6],wrong)),]
train=train[-which(train[,6]=="MCFD2" & train[,4]=="p.S10T"),]
train[which(is.na(train[,70])),70]=0
gnomad[which(is.na(gnomad[,70])),70]=0

colnames(train)[1]="ClinVar"
train[,1]=as.character(train[,1])
train[,1]=as.factor(train[,1])

# compute RandomForest for MutScore

rf <- randomForest(
  ClinVar ~ V20 + V24 + V21 + V28 + V48 + V49 + V50 + V51 + V52 + V53 + V54 + V55 + V17 + V18 + V58 + V59 + V69 + V70,
  data=train, ntree = 1000, importance = TRUE
)
rf.roc<-roc(train$ClinVar,rf$votes[,2])
auc(rf.roc)

write.matrix(rf$importance,file=paste(here,dataLoc,"importance.txt",sep=""))
save(rf, file = paste(here,dataLoc,"rf.RData",sep=""))

# compute RandomForest for MutScore with REVEL

rf <- randomForest(
  ClinVar ~ V20 + V24 + V21 + V28 + V48 + V49 + V50 + V51 + V52 + V53 + V54 + V55 + V17 + V18 + V58 + V59 + V69 + V70 + V33,
  data=train, ntree = 1000, importance = TRUE
)
rf.roc<-roc(train$ClinVar,rf$votes[,2])
auc(rf.roc)

write.matrix(rf$importance,file=paste(here,dataLoc,"importance-REVEL.txt",sep=""))
save(rf, file = paste(here,dataLoc,"rf-revel.RData",sep=""))
train[,74]=rf$votes[,2]

auc(roc(train[,1],train[,74]))


rf <- randomForest(
  ClinVar ~ V20 + V24 + V21 + V28 + V48 + V49 + V50 + V51 + V52 + V53 + V54 + V55 + V17 + V18 + V58 + V59 + V69 + V70 + V29,
  data=train, ntree = 1000, importance = TRUE
)
rf.roc<-roc(train$ClinVar,rf$votes[,2])
auc(rf.roc)

write.matrix(rf$importance,file=paste(here,dataLoc,"importance-VEST4.txt",sep=""))
save(rf, file = paste(here,dataLoc,"rf-VEST4.RData",sep=""))
train[,74]=rf$votes[,2]

auc(roc(train[,1],train[,74]))



#######################

# computing MutScore for all missense (column 74)

load(file = paste(here,dataLoc,"rf.RData",sep=""))

train[,74]=rf$votes[,2]
pred2 = predict(rf, newdata=VUSCON,"prob")
VUSCON[,74]=pred2[,2]
pred2 = predict(rf, newdata=gnomad,"prob")
gnomad[,74]=pred2[,2]
pred2 = predict(rf, newdata=dbnfsp,"prob")
dbnfsp[,74]=pred2[,2]

# unique identifier in column 75 (e.g. "1-897009-A-G")
train[,75]=apply(train[,c(8,9,11,12)],1,paste,collapse="-")
train[,75]=gsub(" ","",train[,75])

VUSCON[,75]=apply(VUSCON[,c(8,9,11,12)],1,paste,collapse="-")
VUSCON[,75]=gsub(" ","",VUSCON[,75])

gnomad[,75]=apply(gnomad[,c(8,9,11,12)],1,paste,collapse="-")
gnomad[,75]=gsub(" ","",gnomad[,75])

dbnfsp[,75]=apply(dbnfsp[,c(8,9,11,12)],1,paste,collapse="-")
dbnfsp[,75]=gsub(" ","",dbnfsp[,75])

tools=matrix(0,ncol=2,nrow=45)
tools[,1]=c("ClinPred","dbscSNV-ADA","dbscSNV-RF","CONDEL","SIFT","SIFT4G","Polyphen2-HDIV","Polyphen2-HVAR","LRT","MutationTaster","MutationAssessor","FATHMM","PROVEAN","VEST4",
	"MetaSVM","MetaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL","fathmm-XF","Eigen","Eigen-PC","GenoCanyon","integrated-fitCons",
	"GERP++_NR","GERP++_RS","phyloP100way_vertebrate","phyloP30way_mammalian","phyloP17way_primate","phastCons100way_vertebrate","phastCons30way_mammalian","phastCons17way_primate",
	"SiPhy_29way","Pext-Mean","Pext-Max","aa-change","Positional-score","MutSore")
tools[,2]=c(16:55,58:59,69,70,74)

save(train,VUSCON,gnomad,dbnfsp,rf,tools, file = paste(here,dataLoc,"data10.RData",sep=""))

save(train,VUSCON,gnomad,dbnfsp,rf,tools, file = paste(here,dataLoc,"data10-freq.RData",sep=""))

############################################

                 # PLOTS #

############################################

load(file = paste(here,dataLoc,"data10.RData",sep=""))

# MutScore for PLP and BLB variants of the training set for the Shiny app
PLP=train[which(train[,1]=="PLP"),74]
BLB=train[which(train[,1]=="BLB"),74]
save(PLP,BLB, file = paste(here,"Shiny/PLP-BLB.RData",sep=""))

####### plot heatmap of correlation between scores #########

tools2=matrix(0,ncol=2,nrow=43)

tools2[,1]=c("ClinPred","dbscSNV-ADA","dbscSNV-RF","CONDEL","SIFT","SIFT4G","Polyphen2-HDIV","Polyphen2-HVAR","LRT","MutationTaster","MutationAssessor","FATHMM","PROVEAN","VEST4",
	"MetaSVM","MetaLR","REVEL","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL","fathmm-XF","Eigen","Eigen-PC","GenoCanyon","integrated-fitCons",
	"GERP++_NR","GERP++_RS","phyloP100way_vertebrate","phyloP30way_mammalian","phyloP17way_primate","phastCons100way_vertebrate","phastCons30way_mammalian","phastCons17way_primate",
	"SiPhy_29way","Pext-Mean","Pext-Max","aa-change","Positional-score","MutSore")
tools2[,2]=c(16:31,33,35:55,58:59,69,70,74)

train2=train
train2[,c(20,21,24,27,28)]=(-1)*train2[,c(20,21,24,27,28)]

colnames(train2)[as.numeric(tools2[,2])]=tools2[,1]
res<-cor(train2[,as.numeric(tools2[,2])],method="spearman")
library("gplots")
pdf(file = paste(here,plot,"heatmap-correlation.pdf",sep=""))
heatmap.2(res, scale = "none", col = bluered(100), 
          trace = "none", density.info = "none",cexRow=0.5,cexCol =0.5,na.rm=T)
dev.off()


#### plot importance of features in MutScore model ###########

pdf(file = paste(here,plot,"features-accuracy.pdf",sep=""),height=6,width=6)

rf$importance[,1]=c("SIFT","LRT","SIFT4G","PROVEAN","GERP++RS","phyloP100way","phyloP30way","phyloP17way","phastCons100way","phastCons30way","phastCons17way",
	"SiPhy29way","dbscSNV-ADA","dbscSNV-RF","pext-mean","pext-max","aa change ratio","Positional score")
par(mai=c(1,3,1,1))
a=sort(as.numeric(rf$importance[,3]),index.return=T)$ix
barplot(as.numeric(rf$importance[a,3]),horiz=T,names.arg=c(rf$importance[a,1]),las=1,cex.names=0.8,main="Mean decrease in accuracy")
dev.off()

pdf(file = paste(here,plot,"features-Gini.pdf",sep=""),height=6,width=6)
par(mai=c(1,3,1,1))
a=sort(as.numeric(rf$importance[,4]),index.return=T)$ix
barplot(as.numeric(rf$importance[a,4]),horiz=T,names.arg=c(rf$importance[a,1]),las=1,cex.names=0.8,main="Mean decrease in Gini index")
dev.off()

###############################

plotROC <- function(data,name,title) {

	pdf(file = paste(here,plot,"",name,".pdf",sep=""),width=10,height=10)

	nbgene=length(unique(data[,6]))
	par(mfrow=c(1,1))
	par(mai=c(1,1,1,1))
	plot(1,1,xaxs="i",yaxs="i",col=0,xlim=c(1.02,-0.02),ylim=c(-0.02,1.02),xlab="Specificity",ylab="Sensitivity",main=paste(title," missense variants in ",nbgene," genes",sep=""),asp=1,cex.main=2,cex.lab=1.8,cex.axis=1.5)
	lines(roc(data[,1],data[,74]),col=1,xaxs="i",yaxs="i",lty=1,lwd=1.5)
	roc1=roc(data[,1],data[,74])
	text(0.2,0.9,"ROC curve for top 8",cex=1.5,font=2)
	names=c("ClinPred","ADA","RF","CONDEL","SIFT","SIFT4G","PolyPhen-HDIV","PolyPhen-HVAR","LRT","MutationTaster",
	      "MutationAssessor","FATHMM","PROVEAN","VEST4","metaSVM","metaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL",
	      "fathmm-XF","eigen","eigen-PC","GenoCanyon","fitCons","GERP++NR","GERP++RS","PhyloP100way","PhyloP30way","PhyloP17way","phastCons100way","phastCons30way","phastCons17way","SiPhy29way")
	nb=c(16:55)
	res=cbind(names,names,names,names)
	res[,2]=nb
	for (i in 1:dim(res)[1]){
		print(i)
		res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
		if(is.element(i,c(5,6,9,12,13))){
			roc2=roc(data$ClinVar,(-1)*data[,as.numeric(res[i,2])])
			res[i,4]=roc.test(roc1, roc2)$p.value
		}
		else{
			roc2=roc(data$ClinVar,data[,as.numeric(res[i,2])])
			res[i,4]=roc.test(roc1, roc2)$p.value
		}
	}
	res1=res
	res2=res[sort(as.numeric(res[,3]),index.return=T,decreasing=T)$ix,]
	res2=res2[-which(res2[,2]==88 | res2[,2]==16),]

	n=7
	colo=c("darkgoldenrod","#FF0018","#FFA52C","#FFFF41","chartreuse3","blue","purple")
	for (k in 1:n){
		i=as.numeric(res2[k,2])
		lines(roc(data$ClinVar,data[,i]),col=colo[k],lty=1,lwd=1.5)
		#text(0.5,0.35-k*0.035,paste("AUC ",res2[k,1]," = ",round(auc(roc(data$ClinVar,data[,i])),digits=3),sep=""))
	}
	legend(0.35,0.87,legend=c("MutScore",res2[1:n,1]),
	     col=c(1,colo),lwd=2,lty=1,cex=1.2)

	n=20
	colo[1:7]=c("darkgoldenrod","#FF0018","#FFA52C","#FFFF41","chartreuse3","blue","purple")
	colo[8:n]="gray"
	par(fig = c(0.2,0.98,0.15,0.65), new = T) 
	xx<-barplot(c(auc(roc(data$ClinVar,data[,74])),as.numeric(res2[1:n,3])),names.arg=c("MutScore",res2[1:n,1]),col=c(1,colo),las=2,ylim=c(0.775,1),xpd=F,main="",cex.names=0.95,cex.axis=1.2,panel.first=c(abline(h=(0.8-0.775)/0.225,col='grey'),abline(h=(0.85-0.775)/0.225,col='grey'),abline(h=(0.9-0.775)/0.225,col='grey'),abline(h=(0.95-0.775)/0.225,col='grey')))
  #grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
	text(13,0.98,"AUC of top 20 predictors",font=2,adj=0.5,cex=1.5)
	label=c(auc(roc(data$ClinVar,data[,74])),as.numeric(res2[1:n,3]))
	label=round(label,digits=3)
	text(x = xx, y = label-0.002, label = label, pos = 3, cex = 0.5)
	box()
	dev.off()
	write.table(res1[,c(1,3,4)],file=paste(here,dataLoc,"AUC-delong-",name,".tsv",sep=""),quote=F,sep="\t")
	print(roc(data[,1],data[,74]))
}

plotROC2 <- function(data,name,title) {

	pdf(file = paste(here,plot,"",name,".pdf",sep=""),width=10,height=10)

	nbgene=length(unique(data[,6]))
	par(mfrow=c(1,1))
	par(mai=c(1,1,1,1))
	plot(1,1,xaxs="i",yaxs="i",col=0,xlim=c(1.02,-0.02),ylim=c(-0.02,1.02),xlab="Specificity",ylab="Sensitivity",main=paste(title," missense variants in ",nbgene," genes",sep=""),asp=1,cex.main=2,cex.lab=1.8,cex.axis=1.5)
	lines(roc(data[,1],data[,74]),col=1,xaxs="i",yaxs="i",lty=1,lwd=1.5)
	roc1=roc(data[,1],data[,74])
	text(0.2,0.9,"ROC curve for top 8",cex=1.5,font=2)
	names=c("ClinPred","ADA","RF","CONDEL","SIFT","SIFT4G","PolyPhen-HDIV","PolyPhen-HVAR","LRT","MutationTaster",
	      "MutationAssessor","FATHMM","PROVEAN","VEST4","metaSVM","metaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL",
	      "fathmm-XF","eigen","eigen-PC","GenoCanyon","fitCons","GERP++NR","GERP++RS","PhyloP100way","PhyloP30way","PhyloP17way","phastCons100way","phastCons30way","phastCons17way","SiPhy29way")
	nb=c(16:55)
	res=cbind(names,names,names,names)
	res[,2]=nb
	for (i in 1:dim(res)[1]){
		print(i)
		res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
		if(is.element(i,c(5,6,9,12,13))){
			roc2=roc(data$ClinVar,(-1)*data[,as.numeric(res[i,2])])
			res[i,4]=roc.test(roc1, roc2)$p.value
		}
		else{
			roc2=roc(data$ClinVar,data[,as.numeric(res[i,2])])
			res[i,4]=roc.test(roc1, roc2)$p.value
		}
	}
	res1=res
	res2=res[sort(as.numeric(res[,3]),index.return=T,decreasing=T)$ix,]
	res2=res2[-which(res2[,2]==88 | res2[,2]==16),]

	n=7
	colo=c("darkgoldenrod","#FF0018","#FFA52C","#FFFF41","chartreuse3","blue","purple")
	for (k in 1:n){
		i=as.numeric(res2[k,2])
		lines(roc(data$ClinVar,data[,i]),col=colo[k],lty=1,lwd=1.5)
		#text(0.5,0.35-k*0.035,paste("AUC ",res2[k,1]," = ",round(auc(roc(data$ClinVar,data[,i])),digits=3),sep=""))
	}
	legend(0.35,0.87,legend=c("MutScore",res2[1:n,1]),
	     col=c(1,colo),lwd=2,lty=1,cex=1.2)

	n=20
	colo[1:7]=c("darkgoldenrod","#FF0018","#FFA52C","#FFFF41","chartreuse3","blue","purple")
	colo[8:n]="gray"
	par(fig = c(0.2,0.98,0.15,0.65), new = T) 
	xx<-barplot(c(auc(roc(data$ClinVar,data[,74])),as.numeric(res2[1:n,3])),names.arg=c("MutScore",res2[1:n,1]),col=c(1,colo),las=2,ylim=c(0.5,0.8),xpd=F,main="",cex.names=0.95,cex.axis=1.2,panel.first=c(abline(h=(0.8-0.775)/0.225,col='grey'),abline(h=(0.85-0.775)/0.225,col='grey'),abline(h=(0.9-0.775)/0.225,col='grey'),abline(h=(0.95-0.775)/0.225,col='grey')))
  #grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted",lwd = par("lwd"), equilogs = TRUE)
	text(13,0.98,"AUC of top 20 predictors",font=2,adj=0.5,cex=1.5)
	label=c(auc(roc(data$ClinVar,data[,74])),as.numeric(res2[1:n,3]))
	label=round(label,digits=3)
	text(x = xx, y = label-0.002, label = label, pos = 3, cex = 0.5)
	box()
	dev.off()
	write.table(res1[,c(1,3,4)],file=paste(here,dataLoc,"AUC-delong-",name,".tsv",sep=""),quote=F,sep="\t")
	print(roc(data[,1],data[,74]))
}



### Training all

train2=train[which((train[,1]=="PLP")),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4

name="Training-all"
title="All ClinVar PLP and BLB"

plotROC(data,name,title)


### Training somatic

train2=train[which((train[,1]=="PLP" & (train[,62]==2 | train[,62]==3))),]
train3=train[which(train[,1]=="BLB" & is.element(train[,6],unique(train2[,6]))),]
train4=rbind(train2,train3)
data=train4

name="Training-somatic"
title="Somatic PLP and all BLB ClinVar"

plotROC(data,name,title)


### Training germline

train2=train[which((train[,1]=="PLP" & train[,62]==1)),]
train3=train[which(train[,1]=="BLB" & is.element(train[,6],unique(train2[,6]))),]
train4=rbind(train2,train3)
data=train4

name="Training-germline"
title="Germline PLP and all BLB ClinVar"

plotROC(data,name,title)


### Training de novo

train2=train[which((train[,1]=="PLP" & (train[,62]==32 | train[,62]==33))),]
train3=train[which(train[,1]=="BLB" & is.element(train[,6],unique(train2[,6]))),]
train4=rbind(train2,train3)
data=train4

name="Training-denovo"
title="De novo PLP and all BLB ClinVar"

plotROC(data,name,title)


### DoCM cancer

cos2<-read.table(file=paste(here,"00_scripts/databases-processing/DoCM.pos.tsv",sep=""),sep="\t",header=F)
list=cos2[,1]
train2=train[which((train[,1]=="PLP" & is.element(train[,75],list))),]
train3=train[which(train[,1]=="BLB" & is.element(train[,6],unique(train2[,6]))),]
train4=rbind(train2,train3)
data=train4

name="Training-DoCM"
title="DoCM cancer driver PLP and BLB training ClinVar"

plotROC(data,name,title)


### 0 stars

train2=train[which(train[,1]=="PLP" & train[,60]=="no_assertion_criteria_provided"),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4

name="Training-0star"
title="0 star PLP and all BLB ClinVar"

plotROC(data,name,title)


### 1 star

train2=train[which(train[,1]=="PLP" & train[,60]=="criteria_provided,_single_submitter"),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4

name="Training-1star"
title="1 star PLP and all BLB ClinVar"

plotROC(data,name,title)


### 2 stars

train2=train[which(train[,1]=="PLP" & train[,60]=="criteria_provided,_multiple_submitters,_no_conflicts"),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4

name="Training-2stars"
title="2 stars PLP and all BLB ClinVar"

plotROC(data,name,title)


### 3 stars

train2=train[which(train[,1]=="PLP" & train[,60]=="reviewed_by_expert_panel"),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4

name="Training-3stars"
title="3 stars PLP and all BLB ClinVar"

plotROC(data,name,title)


### without positional

train2=train[which((train[,1]=="PLP")),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4[which(train4[,72]=="LOW"),]

name="Training-WO-positional"
title="No positional, PLP and all BLB training ClinVar"

plotROC(data,name,title)


### with positional

train2=train[which((train[,1]=="PLP")),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4[which(train4[,72]==0),]

name="Training-WITH-positional"
title="With positional, PLP and all BLB training ClinVar"

plotROC(data,name,title)


### with positional

train2=train[which((train[,1]=="PLP")),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4[which(train4[,71]>=20),]

name="Training-WITH-20PLPs-more"
title="PLP and all BLB training ClinVar more"

plotROC(data,name,title)

### with positional

train2=train[which((train[,1]=="PLP")),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4[which(train4[,71]<20),]

name="Training-WITH-20PLPs-less"
title="PLP and all BLB training ClinVar less"

plotROC(data,name,title)

### with positional

train2=train[which((train[,1]=="PLP")),]
train3=train[which(train[,1]=="BLB"),]
train4=rbind(train2,train3)
data=train4[which(train4[,71]<20 & train4[,71]>0),]

name="Training-WITH-10-20PLPs"
title="PLP and all BLB training ClinVar less"

plotROC(data,name,title)




############ Figure 3 and SX ####

load(file = paste(here,dataLoc,"regions-review.RData",sep=""))

##### plot for PLP regions

# resi
# 1=isoform / 2=gene / 3=%PLP-regionS / 4=%gene-regionS / 5=score / 6=Fisher-estimate / 7=Fisher-p.value

resi=resi1

a5=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])>=0.66),2])
a4=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.66 & as.numeric(resi[,5])>=0.33),2])
a3=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.33),2])
a2=unique(resi[which(as.numeric(resi[,8])>=0.05),2])
a1=unique(train[which(train[,70]>0),6])

a4=setdiff(a4,a5)
a3=setdiff(a3,c(a4,a5))
a2=setdiff(a2,c(a3,a4,a5))
a1=setdiff(a1,c(a2,a3,a4,a5))

low=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])<0.33)
med=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.33 & as.numeric(resi[,5])<0.66)
high=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.66)

pdf(file = paste(here,plot,"clustering-PLP-pie-review.pdf",sep=""))
slices <- c(length(a1),length(a2),length(a3),length(a4),length(a5))
lbls <- c("No regions with >5 PLPs","> 1 region and \n not significant","Sparse clustering","Medium clustering","Dense clustering")
a=paste("n = ",slices,sep="")
tot=length(a1)+length(a2)+length(a3)+length(a4)+length(a5)
perc <- paste(round(slices/tot*100,digits=1),"%",sep="")
pie(slices, labels = paste(lbls,"\n",a," / ",perc,sep=""), main="PLP clustering categories for most clustered isoform per gene", col=c("lightgrey","darkgrey","red","red3","red4"), border="white",init.angle=90,cex=0.5)
dev.off()

pdf(file = paste(here,plot,"clustering-PLP-points-review.pdf",sep=""))
par(pty="s")
plot(as.numeric(resi[,3]),1-as.numeric(resi[,4]),col="grey",xaxs="i",yaxs="i",pch=16,cex=0.6,ylim=c(-0.02,1.02),asp=1,xlim=c(-0.02,1.02),xlab="Clustering precision",ylab="Clustering density",main="Clustering of PLP variants")
points(as.numeric(resi[low,3]),1-as.numeric(resi[low,4]),col="red",pch=16,cex=0.6)
points(as.numeric(resi[med,3]),1-as.numeric(resi[med,4]),col="red3",pch=16,cex=0.6)
points(as.numeric(resi[high,3]),1-as.numeric(resi[high,4]),col="red4",pch=16,cex=0.6)
lines(c(0.33,0.33),c(0.33,2))
lines(c(0.66,0.66),c(0.66,2))
lines(c(0.66,2),c(0.66,0.66))
lines(c(0.33,2),c(0.33,0.33))
legend(0.1,0.25,legend=c("High clustering","Medium clustering","Low clustering","Non significant"),col=c("red4","red3","red","darkgrey"),pch=19,lty=0,cex=0.8)
dev.off()


##### plot for PLP regions AD vs AR

load(file = paste(here,dataLoc,"regions-review.RData",sep=""))

dom<-read.table(file=paste(here,"00_scripts/databases-processing/gene-dom-pure.tsv",sep=""),header=FALSE,sep="\t")[,1]
rec<-read.table(file=paste(here,"00_scripts/databases-processing/gene-rec-pure.tsv",sep=""),header=FALSE,sep="\t")[,1]
som<-read.table(file=paste(here,"00_scripts/databases-processing/gene-som-pure.tsv",sep=""),header=FALSE,sep="\t")[,1]

dom1=which(is.element(resi1[,2],dom))
rec1=which(is.element(resi1[,2],rec))
som1=which(is.element(resi1[,2],som))

resi=resi1[dom1,]

write.table(unique(resi[,2]),file=paste(here,"00_scripts/databases-processing/gene-dom-used.tsv",sep=""),quote=F,sep="\t")

a5=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])>=0.66),2])
a4=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.66 & as.numeric(resi[,5])>=0.33),2])
a3=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.33),2])
a2=unique(resi[which(as.numeric(resi[,8])>=0.05),2])
a1=unique(train[which(train[,70]>0 & is.element(train[,6],dom)),6])

a4=setdiff(a4,a5)
a3=setdiff(a3,c(a4,a5))
a2=setdiff(a2,c(a3,a4,a5))
a1=setdiff(a1,c(a2,a3,a4,a5))

low=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])<0.33)
med=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.33 & as.numeric(resi[,5])<0.66)
high=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.66)

pdf(file = paste(here,plot,"clustering-PLP-pie-AD-review.pdf",sep=""))
slices <- c(length(a1),length(a2),length(a3),length(a4),length(a5))
lbls <- c("No regions with >5 PLPs","> 1 region and \n not significant","Sparse clustering","Medium clustering","Dense clustering")
a=paste("n = ",slices,sep="")
tot=length(a1)+length(a2)+length(a3)+length(a4)+length(a5)
perc <- paste(round(slices/tot*100,digits=1),"%",sep="")
pie(slices, labels = paste(lbls,"\n",a," / ",perc,sep=""), main="", col=c("lightgrey","darkgrey","red","red3","red4"), border="white",init.angle=90,cex=1)
dev.off()


resi=resi1[rec1,]

write.table(unique(resi[,2]),file=paste(here,"00_scripts/databases-processing/gene-rec-used.tsv",sep=""),quote=F,sep="\t")

a5=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])>=0.66),2])
a4=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.66 & as.numeric(resi[,5])>=0.33),2])
a3=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.33),2])
a2=unique(resi[which(as.numeric(resi[,8])>=0.05),2])
a1=unique(train[which(train[,70]>0 & is.element(train[,6],rec)),6])

a4=setdiff(a4,a5)
a3=setdiff(a3,c(a4,a5))
a2=setdiff(a2,c(a3,a4,a5))
a1=setdiff(a1,c(a2,a3,a4,a5))

low=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])<0.33)
med=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.33 & as.numeric(resi[,5])<0.66)
high=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.66)

pdf(file = paste(here,plot,"clustering-PLP-pie-AR-review.pdf",sep=""))
slices <- c(length(a1),length(a2),length(a3),length(a4),length(a5))
lbls <- c("No regions with >5 PLPs","> 1 region and \n not significant","Sparse clustering","Medium clustering","Dense clustering")
a=paste("n = ",slices,sep="")
tot=length(a1)+length(a2)+length(a3)+length(a4)+length(a5)
perc <- paste(round(slices/tot*100,digits=1),"%",sep="")
pie(slices, labels = paste(lbls,"\n",a," / ",perc,sep=""), main="", col=c("lightgrey","darkgrey","red","red3","red4"), border="white",init.angle=90,cex=1)
dev.off()

resi=resi1[som1,]

write.table(unique(resi[,2]),file=paste(here,"00_scripts/databases-processing/gene-som-used.tsv",sep=""),quote=F,sep="\t")

a5=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])>=0.66),2])
a4=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.66 & as.numeric(resi[,5])>=0.33),2])
a3=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.33),2])
a2=unique(resi[which(as.numeric(resi[,8])>=0.05),2])
a1=unique(train[which(train[,70]>0 & is.element(train[,6],som)),6])

a4=setdiff(a4,a5)
a3=setdiff(a3,c(a4,a5))
a2=setdiff(a2,c(a3,a4,a5))
a1=setdiff(a1,c(a2,a3,a4,a5))

low=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])<0.33)
med=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.33 & as.numeric(resi[,5])<0.66)
high=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.66)

pdf(file = paste(here,plot,"clustering-PLP-pie-SOM-review.pdf",sep=""))
slices <- c(length(a1),length(a2),length(a3),length(a4),length(a5))
lbls <- c("No regions with >5 PLPs","> 1 region and \n not significant","Sparse clustering","Medium clustering","Dense clustering")
a=paste("n = ",slices,sep="")
tot=length(a1)+length(a2)+length(a3)+length(a4)+length(a5)
perc <- paste(round(slices/tot*100,digits=1),"%",sep="")
pie(slices, labels = paste(lbls,"\n",a," / ",perc,sep=""), main="", col=c("lightgrey","darkgrey","red","red3","red4"), border="white",init.angle=90,cex=1)
dev.off()


# plot cluster BLB 

load(file = paste(here,dataLoc,"regions-review.RData",sep=""))

#load(file = paste(here,dataLoc,"regions.RData",sep=""))

resi=resi2

a5=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])>=0.66),2])
a4=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.66 & as.numeric(resi[,5])>=0.33),2])
a3=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.33),2])
a2=unique(resi[which(as.numeric(resi[,8])>=0.05 | is.na(resi[,7])),2])
a1=unique(train[which(train[,70]>0),6])

length(a1)
length(a2)
length(a3)
length(a4)
length(a5)

not=setdiff(unique(train[,6]),a1)

a4=setdiff(a4,c(a5,not))
a3=setdiff(a3,c(a4,a5,not))
a2=setdiff(a2,c(a3,a4,a5,not))
a1=setdiff(a1,c(a2,a3,a4,a5,not))

low=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])<0.33)
med=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.33 & as.numeric(resi[,5])<0.66)
high=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.66)

pdf(file = paste(here,plot,"clustering-BLB-pie-review.pdf",sep=""))
slices <- c(length(a1),length(a2),length(a3),length(a4),length(a5))
lbls <- c("No regions with >10 BLBs","> 1 region and \n not significant","Sparse clustering","Medium clustering","Dense clustering")
a=paste("n = ",slices,sep="")
tot=length(a1)+length(a2)+length(a3)+length(a4)+length(a5)
perc <- paste(round(slices/tot*100,digits=1),"%",sep="")
pie(slices, labels = paste(lbls,"\n",a," / ",perc,sep=""), main="", col=c("lightgrey","darkgrey","dodgerblue2","blue1","blue4"), border="white",init.angle=90,cex=1)
dev.off()

pdf(file = paste(here,plot,"clustering-BLB-points-review.pdf",sep=""))
par(pty="s")
plot(as.numeric(resi[,3]),1-as.numeric(resi[,4]),col="grey",xaxs="i",yaxs="i",pch=16,cex=0.6,ylim=c(-0.02,1.02),asp=1,xlim=c(-0.02,1.02),xlab="Clustering precision",ylab="Clustering density",main="Clustering of BLB variants")
points(as.numeric(resi[low,3]),1-as.numeric(resi[low,4]),col="dodgerblue2",pch=16,cex=0.6)
points(as.numeric(resi[med,3]),1-as.numeric(resi[med,4]),col="blue1",pch=16,cex=0.6)
points(as.numeric(resi[high,3]),1-as.numeric(resi[high,4]),col="blue4",pch=16,cex=0.6)
lines(c(0.33,0.33),c(0.33,2))
lines(c(0.66,0.66),c(0.66,2))
lines(c(0.66,2),c(0.66,0.66))
lines(c(0.33,2),c(0.33,0.33))
legend(0.1,0.25,legend=c("High clustering","Medium clustering","Low clustering","Non significant"),col=c("blue4","blue1","dodgerblue2","darkgrey"),pch=19,lty=0,cex=0.8)
dev.off()

##### plot for BLB regions AD vs AR

load(file = paste(here,dataLoc,"regions-review.RData",sep=""))

dom<-read.table(file=paste(here,"00_scripts/databases-processing/gene-dom-pure.tsv",sep=""),header=FALSE,sep="\t")[,1]
rec<-read.table(file=paste(here,"00_scripts/databases-processing/gene-rec-pure.tsv",sep=""),header=FALSE,sep="\t")[,1]
som<-read.table(file=paste(here,"00_scripts/databases-processing/gene-som-pure.tsv",sep=""),header=FALSE,sep="\t")[,1]

dom1=which(is.element(resi2[,2],dom))
rec1=which(is.element(resi2[,2],rec))
som1=which(is.element(resi2[,2],som))

resi=resi2[dom1,]

a5=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])>=0.66),2])
a4=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.66 & as.numeric(resi[,5])>=0.33),2])
a3=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.33),2])
a2=unique(resi[which(as.numeric(resi[,8])>=0.05 | is.na(resi[,7])),2])
a1=unique(train[which(train[,70]>0 & is.element(train[,6],dom)),6])

not=setdiff(unique(train[,6]),a1)

a4=setdiff(a4,c(a5,not))
a3=setdiff(a3,c(a4,a5,not))
a2=setdiff(a2,c(a3,a4,a5,not))
a1=setdiff(a1,c(a2,a3,a4,a5,not))

low=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])<0.33)
med=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.33 & as.numeric(resi[,5])<0.66)
high=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.66)

pdf(file = paste(here,plot,"clustering-BLB-pie-AD-review.pdf",sep=""))
slices <- c(length(a1),length(a2),length(a3),length(a4),length(a5))
lbls <- c("No regions with >10 BLBs","> 1 region and \n not significant","Sparse clustering","Medium clustering","Dense clustering")
a=paste("n = ",slices,sep="")
tot=length(a1)+length(a2)+length(a3)+length(a4)+length(a5)
perc <- paste(round(slices/tot*100,digits=1),"%",sep="")
pie(slices, labels = paste(lbls,"\n",a," / ",perc,sep=""), main="", col=c("lightgrey","darkgrey","dodgerblue2","blue1","blue4"), border="white",init.angle=90,cex=1)
dev.off()


resi=resi2[rec1,]

a5=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])>=0.66),2])
a4=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.66 & as.numeric(resi[,5])>=0.33),2])
a3=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.33),2])
a2=unique(resi[which(as.numeric(resi[,8])>=0.05 | is.na(resi[,7])),2])
a1=unique(train[which(train[,70]>0 & is.element(train[,6],rec)),6])

not=setdiff(unique(train[,6]),a1)

a4=setdiff(a4,c(a5,not))
a3=setdiff(a3,c(a4,a5,not))
a2=setdiff(a2,c(a3,a4,a5,not))
a1=setdiff(a1,c(a2,a3,a4,a5,not))

low=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])<0.33)
med=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.33 & as.numeric(resi[,5])<0.66)
high=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.66)

pdf(file = paste(here,plot,"clustering-BLB-pie-AR-review.pdf",sep=""))
slices <- c(length(a1),length(a2),length(a3),length(a4),length(a5))
lbls <- c("No regions with >10 BLBs","> 1 region and \n not significant","Sparse clustering","Medium clustering","Dense clustering")
a=paste("n = ",slices,sep="")
tot=length(a1)+length(a2)+length(a3)+length(a4)+length(a5)
perc <- paste(round(slices/tot*100,digits=1),"%",sep="")
pie(slices, labels = paste(lbls,"\n",a," / ",perc,sep=""), main="", col=c("lightgrey","darkgrey","dodgerblue2","blue1","blue4"), border="white",init.angle=90,cex=1)
dev.off()

resi=resi2[som1,]

a5=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])>=0.66),2])
a4=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.66 & as.numeric(resi[,5])>=0.33),2])
a3=unique(resi[which(as.numeric(resi[,8])<0.05 & as.numeric(resi[,5])<0.33),2])
a2=unique(resi[which(as.numeric(resi[,8])>=0.05 | is.na(resi[,7])),2])
a1=unique(train[which(train[,70]>0 & is.element(train[,6],som)),6])

not=setdiff(unique(train[,6]),a1)

a4=setdiff(a4,c(a5,not))
a3=setdiff(a3,c(a4,a5,not))
a2=setdiff(a2,c(a3,a4,a5,not))
a1=setdiff(a1,c(a2,a3,a4,a5,not))

low=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])<0.33)
med=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.33 & as.numeric(resi[,5])<0.66)
high=which(as.numeric(resi[,8])<(0.05) & as.numeric(resi[,5])>0.66)

pdf(file = paste(here,plot,"clustering-BLB-pie-SOM-review.pdf",sep=""))
slices <- c(length(a1),length(a2),length(a3),length(a4),length(a5))
lbls <- c("No regions with >10 BLBs","> 1 region and \n not significant","Sparse clustering","Medium clustering","Dense clustering")
a=paste("n = ",slices,sep="")
tot=length(a1)+length(a2)+length(a3)+length(a4)+length(a5)
perc <- paste(round(slices/tot*100,digits=1),"%",sep="")
pie(slices, labels = paste(lbls,"\n",a," / ",perc,sep=""), main="", col=c("lightgrey","darkgrey","dodgerblue2","blue1","blue4"), border="white",init.angle=90,cex=1)
dev.off()

######## Fig. S3 (startified by AF) ########

names=c("ClinPred","ADA","RF","CONDEL","SIFT","SIFT4G","PolyPhen-HDIV","PolyPhen-HVAR","LRT","MutationTaster",
        "MutationAssessor","FATHMM","PROVEAN","VEST4","metaSVM","metaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL",
        "fathmm-XF","eigen","eigen-PC","GenoCanyon","fitCons","GERP++NR","GERP++RS","PhyloP100way","PhyloP30way","PhyloP17way","phastCons100way","phastCons30way","phastCons17way","SiPhy29way","MutScore")
nb=c(16:55,74)
res=cbind(names,names,names)
res[,2]=nb

data=train
for (i in 1:dim(res)[1]){
  res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
}

sepa=c(0,3.5,4,4.5,5,7)*(-1)
for (j in 1:(length(sepa)-1)){
  res=cbind(res,names)
  data=train[which(train[,56]<sepa[j] & train[,56]>sepa[j+1]),]
  for (i in 1:dim(res)[1]){
    res[i,j+3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
  }
}

res[sort(as.numeric(res[,3]),index.return=T)$ix,]
res1=res[sort(as.numeric(res[,3]),index.return=T)$ix,]
write.table(res1[,c(1,3:8)],file=paste(here,dataLoc,"AUC-freq.tsv",sep=""),quote=F,sep="\t")

pdf(file = paste(here,plot,"Freq.pdf",sep=""),width=10,height=10)
plot(as.numeric(res1[41,3:8]),ylim=c(0.82,0.98),xlim=c(0.7,6.3),xaxt="n",pch=19,col=6,xlab="-log10(AF range of variants)",ylab="AUC",main="AUC for the training set stratified by AF")
lines(2:6,as.numeric(res1[41,4:8]),col=6)
axis(1, at=1:6, labels=c("All","<3.5","3.5-4.0","4.0-4.5","4.5-5.0",">5.0"))
lines(c(1.5,1.5),c(0.75,1))
for (i in 38:40){
  points(as.numeric(res1[i,3:8]),col=i-36,pch=19)
  lines(2:6,as.numeric(res1[i,4:8]),col=i-36)
}
legend(4.5,0.84,legend=c("MutScore","REVEL","VEST4","ClinPred"),col=c(6,4,3,2),pch=19,lty=1)
dev.off()

################### Fig. S6 VUS and conflicting ##############################################

# VUS
genes=unique(train[which(train[,1]=="PLP"),6])
VUS=VUSCON[which(VUSCON[,1]=="VUS" & is.element(VUSCON[,6],genes)),]
CON=VUSCON[which(VUSCON[,1]=="CON" & is.element(VUSCON[,6],genes)),]

pred = cbind(train[,74],train[,74])
pred2 = cbind(VUS[,74],VUS[,74])
pred3 = cbind(CON[,74],CON[,74])
pred4 = cbind(gnomad[,74],gnomad[,74])

# MutScore
limits=seq(0,1,0.001)
lim1=limits
lim2=limits
for (i in 1:length(limits)){
  lim1[i]=length(which(train[,1]=="PLP" & pred[,2]>limits[i]))/(length(which(train[,1]=="PLP" & pred[,2]>limits[i])) + length(which(train[,1]=="BLB" & pred[,2]>limits[i])))
  lim2[i]=length(which(train[,1]=="BLB" & pred[,2]<limits[i]))/(length(which(train[,1]=="PLP" & pred[,2]<limits[i])) + length(which(train[,1]=="BLB" & pred[,2]<limits[i])))

}
muthigh=min(limits[which(lim1>0.95)])-0.001
mutlow=max(limits[which(lim2>0.95)])+0.001

pdf(file = paste(here,plot,"VUSCON-MutScore.pdf",sep=""),width=4,height=5)
par(mfrow=c(4,1),mar=c(2.1, 4.1, 2.1, 2.1))
hist(pred[which(train[,1]=="PLP"),2],n=50,main="",xlab="",cex.main=2,cex.lab=1.2,col='red',ylim=c(0,12500),yaxt="n")
axis(2, at=c(0,3000,6000,9000,12000),labels=c(0,3000,6000,9000,12000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(pred[which(train[,1]=="PLP"),2]>muthigh))
a2=length(which(pred[which(train[,1]=="PLP"),2]<mutlow))
a3=length(which(pred[which(train[,1]=="PLP"),2]>=mutlow & pred[which(train[,1]=="PLP"),2]<=muthigh))
text(mutlow/2,12000,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,12000,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,12000,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)

hist(pred[which(train[,1]=="BLB"),2],n=50,main="",xlab="MutScore",cex.main=2,cex.lab=1.2,col='blue',ylim=c(0,4300),yaxt="n")
axis(2, at=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(pred[which(train[,1]=="BLB"),2]>muthigh))
a2=length(which(pred[which(train[,1]=="BLB"),2]<mutlow))
a3=length(which(pred[which(train[,1]=="BLB"),2]>=mutlow & pred[which(train[,1]=="BLB"),2]<=muthigh))
text(mutlow/2,4150,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,4150,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,4150,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)

x=seq(0,1,1/50)
x=x[2:length(x)]
col=1:50
col[1:50]='orange'
col[which(x<mutlow)]='deepskyblue'
col[which(x>muthigh)]='indianred2'

hist(pred2[,2],n=50,main="",xlab="MutScore",cex.main=2,cex.lab=1.2,col=col,ylim=c(0,16000),yaxt="n")
axis(2, at=c(0,5000,10000,15000),labels=c(0,5000,10000,15000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(pred2[,2]>muthigh))
a2=length(which(pred2[,2]<mutlow))
a3=length(which(pred2[,2]>=mutlow & pred2[,2]<=muthigh))
text(mutlow/2,15500,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,15500,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,15500,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)


hist(pred3[,2],n=50,main="",xlab="MutScore",cex.main=2,cex.lab=1.2,col=col,ylim=c(0,2100),yaxt="n")
axis(2, at=c(0,500,1000,1500,2000),labels=c(0,500,1000,1500,2000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(pred3[,2]>muthigh))
a2=length(which(pred3[,2]<mutlow))
a3=length(which(pred3[,2]>=mutlow & pred3[,2]<=muthigh))
text(mutlow/2,2000,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,2000,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,2000,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
dev.off()


pdf(file = paste(here,plot,"VUSCON-MutScore3.pdf",sep=""),width=10,height=5)
a1=length(which(pred2[,2]>muthigh))
a2=length(which(pred2[,2]<mutlow))
a3=length(which(pred2[,2]>=mutlow & pred2[,2]<=muthigh))
par(mfrow=c(1,2))
xx<-barplot(c(a2,a3,a1),names.arg=c("Likely benign","VUS","Likely pathogenic"),col=c('deepskyblue','orange','indianred2'),ylab="Number of VUS",main="VUS",ylim=c(0,145000),cex.main=2.5)
label=c(paste("n=",length(which(pred2[,2]<mutlow)),sep=""),paste("n=",dim(VUS)[1]-length(which(pred2[,2]<mutlow))-length(which(pred2[,2]>muthigh)),sep=""),paste("n=",length(which(pred2[,2]>muthigh)),sep=""))
text(x = xx, y = c(a2,a3,a1), label = label, pos = 3, cex = 0.8)
label=c(paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""))
text(x = xx, y = c(a2,a3,a1)+6000, label = label, pos = 3, cex = 0.8)

a1=length(which(pred3[,2]>muthigh))
a2=length(which(pred3[,2]<mutlow))
a3=length(which(pred3[,2]>=mutlow & pred3[,2]<=muthigh))
xx<-barplot(c(a2,a3,a1),names.arg=c("Likely benign","VUS","Likely pathogenic"),col=c('deepskyblue','orange','indianred2'),ylab="Number of CI variants",main="CI variants",ylim=c(0,14000),cex.main=2.5)
label=c(paste("n=",length(which(pred3[,2]<mutlow)),sep=""),paste("n=",dim(CON)[1]-length(which(pred3[,2]<mutlow))-length(which(pred3[,2]>muthigh)),sep=""),paste("n=",length(which(pred3[,2]>muthigh)),sep=""))
text(x = xx, y = c(a2,a3,a1), label = label, pos = 3, cex = 0.8)
label=c(paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""))
text(x = xx, y = c(a2,a3,a1)+600, label = label, pos = 3, cex = 0.8)
dev.off()

# VEST4
n=29
predictor="VEST4"
limits=seq(0,1,0.001)
lim1=limits
lim2=limits
for (i in 1:length(limits)){
  lim1[i]=length(which(train[,1]=="PLP" & train[,n]>limits[i]))/(length(which(train[,1]=="PLP" & train[,n]>limits[i])) + length(which(train[,1]=="BLB" & train[,n]>limits[i])))
  lim2[i]=length(which(train[,1]=="BLB" & train[,n]<limits[i]))/(length(which(train[,1]=="PLP" & train[,n]<limits[i])) + length(which(train[,1]=="BLB" & train[,n]<limits[i])))
  
}
muthigh=min(limits[which(lim1>0.95)])-0.001
mutlow=max(limits[which(lim2>0.9498)])+0.001

pdf(file = paste(here,plot,"VUSCON-VEST4.pdf",sep=""),width=4,height=5)
par(mfrow=c(4,1),mar=c(2.1, 4.1, 2.1, 2.1))
hist(train[which(train[,1]=="PLP"),n],n=50,main="",xlab="",cex.main=2,cex.lab=1.2,col='red',ylim=c(0,12500),yaxt="n")
axis(2, at=c(0,3000,6000,9000,12000),labels=c(0,3000,6000,9000,12000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(train[which(train[,1]=="PLP"),n]>muthigh))
a2=length(which(train[which(train[,1]=="PLP"),n]<mutlow))
a3=length(which(train[which(train[,1]=="PLP"),n]>=mutlow & train[which(train[,1]=="PLP"),n]<=muthigh))
text(mutlow/2,12000,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,12000,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,12000,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)


hist(train[which(train[,1]=="BLB"),n],n=50,main="",xlab=predictor,cex.main=2,cex.lab=1.2,col='blue',ylim=c(0,4300),yaxt="n")
axis(2, at=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(train[which(train[,1]=="BLB"),n]>muthigh))
a2=length(which(train[which(train[,1]=="BLB"),n]<mutlow))
a3=length(which(train[which(train[,1]=="BLB"),n]>=mutlow & train[which(train[,1]=="BLB"),n]<=muthigh))
text(mutlow/2,4150,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,4150,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,4150,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)


x=seq(0,1,1/50)
x=x[2:length(x)]
col=1:50
col[1:50]='orange'
col[which(x<mutlow)]='deepskyblue'
col[which(x>muthigh)]='indianred2'

hist(VUS[,n],n=50,main="",xlab=predictor,cex.main=2,cex.lab=1.2,col=col,ylim=c(0,16000),yaxt="n")
axis(2, at=c(0,5000,10000,15000),labels=c(0,5000,10000,15000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(VUS[,n]>muthigh))
a2=length(which(VUS[,n]<mutlow))
a3=length(which(VUS[,n]>=mutlow & VUS[,n]<=muthigh))
text(mutlow/2,15500,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,15500,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,15500,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)

hist(CON[,n],n=50,main="",xlab=predictor,cex.main=2,cex.lab=1.2,col=col,ylim=c(0,2100),yaxt="n")
axis(2, at=c(0,500,1000,1500,2000),labels=c(0,500,1000,1500,2000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(CON[,n]>muthigh))
a2=length(which(CON[,n]<mutlow))
a3=length(which(CON[,n]>=mutlow & CON[,n]<=muthigh))
text(mutlow/2,2000,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,2000,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,2000,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
dev.off()

pdf(file = paste(here,plot,"VUSCON-VEST43.pdf",sep=""),width=10,height=5)
a1=length(which(VUS[,n]>muthigh))
a2=length(which(VUS[,n]<mutlow))
a3=length(which(VUS[,n]>=mutlow & VUS[,n]<=muthigh))
par(mfrow=c(1,2))
xx<-barplot(c(a2,a3,a1),names.arg=c("Likely benign","VUS","Likely pathogenic"),col=c('deepskyblue','orange','indianred2'),ylab="Number of VUS",main="VUS",ylim=c(0,145000),cex.main=2.5)
label=c(paste("n=",length(which(VUS[,n]<mutlow)),sep=""),paste("n=",dim(VUS)[1]-length(which(VUS[,n]<mutlow))-length(which(VUS[,n]>muthigh)),sep=""),paste("n=",length(which(VUS[,n]>muthigh)),sep=""))
text(x = xx, y = c(a2,a3,a1), label = label, pos = 3, cex = 0.8)
label=c(paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""))
text(x = xx, y = c(a2,a3,a1)+6000, label = label, pos = 3, cex = 0.8)

a1=length(which(CON[,n]>muthigh))
a2=length(which(CON[,n]<mutlow))
a3=length(which(CON[,n]>=mutlow & CON[,n]<=muthigh))
xx<-barplot(c(a2,a3,a1),names.arg=c("Likely benign","VUS","Likely pathogenic"),col=c('deepskyblue','orange','indianred2'),ylab="Number of CI variants",main="CI variants",ylim=c(0,14000),cex.main=2.5)
label=c(paste("n=",length(which(CON[,n]<mutlow)),sep=""),paste("n=",dim(CON)[1]-length(which(CON[,n]<mutlow))-length(which(CON[,n]>muthigh)),sep=""),paste("n=",length(which(CON[,n]>muthigh)),sep=""))
text(x = xx, y = c(a2,a3,a1), label = label, pos = 3, cex = 0.8)
label=c(paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""))
text(x = xx, y = c(a2,a3,a1)+600, label = label, pos = 3, cex = 0.8)
dev.off()

# REVEL
n=33
predictor="REVEL"
limits=seq(0,1,0.001)
lim1=limits
lim2=limits
for (i in 1:length(limits)){
  lim1[i]=length(which(train[,1]=="PLP" & train[,n]>limits[i]))/(length(which(train[,1]=="PLP" & train[,n]>limits[i])) + length(which(train[,1]=="BLB" & train[,n]>limits[i])))
  lim2[i]=length(which(train[,1]=="BLB" & train[,n]<limits[i]))/(length(which(train[,1]=="PLP" & train[,n]<limits[i])) + length(which(train[,1]=="BLB" & train[,n]<limits[i])))
  
}
muthigh=min(limits[which(lim1>0.95)])-0.001
mutlow=max(limits[which(lim2>0.95)])+0.001

pdf(file = paste(here,plot,"VUSCON-REVEL.pdf",sep=""),width=4,height=5)
par(mfrow=c(4,1),mar=c(2.1, 4.1, 2.1, 2.1))
hist(train[which(train[,1]=="PLP"),n],n=50,main="",xlab=predictor,cex.main=2,cex.lab=1.2,col='red',ylim=c(0,12500),yaxt="n")
axis(2, at=c(0,3000,6000,9000,12000),labels=c(0,3000,6000,9000,12000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(train[which(train[,1]=="PLP"),n]>muthigh))
a2=length(which(train[which(train[,1]=="PLP"),n]<mutlow))
a3=length(which(train[which(train[,1]=="PLP"),n]>=mutlow & train[which(train[,1]=="PLP"),n]<=muthigh))
text(mutlow/2,12000,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,12000,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,12000,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)

hist(train[which(train[,1]=="BLB"),n],n=50,main="",xlab=predictor,cex.main=2,cex.lab=1.2,col='blue',ylim=c(0,4300),yaxt="n")
axis(2, at=c(0,1000,2000,3000,4000),labels=c(0,1000,2000,3000,4000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(train[which(train[,1]=="BLB"),n]>muthigh))
a2=length(which(train[which(train[,1]=="BLB"),n]<mutlow))
a3=length(which(train[which(train[,1]=="BLB"),n]>=mutlow & train[which(train[,1]=="BLB"),n]<=muthigh))
text(mutlow/2,4150,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,4150,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,4150,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)

x=seq(0,1,1/50)
x=x[2:length(x)]
col=1:50
col[1:50]='orange'
col[which(x<mutlow)]='deepskyblue'
col[which(x>muthigh)]='indianred2'

hist(VUS[,n],n=50,main="",xlab=predictor,cex.main=2,cex.lab=1.2,col=col,ylim=c(0,16000),yaxt="n")
axis(2, at=c(0,5000,10000,15000),labels=c(0,5000,10000,15000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(VUS[,n]>muthigh))
a2=length(which(VUS[,n]<mutlow))
a3=length(which(VUS[,n]>=mutlow & VUS[,n]<=muthigh))
text(mutlow/2,15500,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,15500,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,15500,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)

hist(CON[,n],n=50,main="",xlab=predictor,cex.main=2,cex.lab=1.2,col=col,ylim=c(0,2100),yaxt="n")
axis(2, at=c(0,500,1000,1500,2000),labels=c(0,500,1000,1500,2000), las=2,cex.axis=1,pos=0)
abline(v=mutlow,col=1)
abline(v=muthigh,col=1)
a1=length(which(CON[,n]>muthigh))
a2=length(which(CON[,n]<mutlow))
a3=length(which(CON[,n]>=mutlow & CON[,n]<=muthigh))
text(mutlow/2,2000,paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((muthigh-mutlow)/2+mutlow,2000,paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
text((1-muthigh)/2+muthigh,2000,paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""),cex=0.9)
dev.off()

pdf(file = paste(here,plot,"VUSCON-REVEL3.pdf",sep=""),width=10,height=5)
a1=length(which(VUS[,n]>muthigh))
a2=length(which(VUS[,n]<mutlow))
a3=length(which(VUS[,n]>=mutlow & VUS[,n]<=muthigh))
par(mfrow=c(1,2))
xx<-barplot(c(a2,a3,a1),names.arg=c("Likely benign","VUS","Likely pathogenic"),col=c('deepskyblue','orange','indianred2'),ylab="Number of VUS",main="VUS",ylim=c(0,145000),cex.main=2.5)
label=c(paste("n=",length(which(VUS[,n]<mutlow)),sep=""),paste("n=",dim(VUS)[1]-length(which(VUS[,n]<mutlow))-length(which(VUS[,n]>muthigh)),sep=""),paste("n=",length(which(VUS[,n]>muthigh)),sep=""))
text(x = xx, y = c(a2,a3,a1), label = label, pos = 3, cex = 0.8)
label=c(paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""))
text(x = xx, y = c(a2,a3,a1)+6000, label = label, pos = 3, cex = 0.8)

a1=length(which(CON[,n]>muthigh))
a2=length(which(CON[,n]<mutlow))
a3=length(which(CON[,n]>=mutlow & CON[,n]<=muthigh))
xx<-barplot(c(a2,a3,a1),names.arg=c("Likely benign","VUS","Likely pathogenic"),col=c('deepskyblue','orange','indianred2'),ylab="Number of CI variants",main="CI variants",ylim=c(0,14000),cex.main=2.5)
label=c(paste("n=",length(which(CON[,n]<mutlow)),sep=""),paste("n=",dim(CON)[1]-length(which(CON[,n]<mutlow))-length(which(CON[,n]>muthigh)),sep=""),paste("n=",length(which(CON[,n]>muthigh)),sep=""))
text(x = xx, y = c(a2,a3,a1), label = label, pos = 3, cex = 0.8)
label=c(paste(round(a2/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a3/(a1+a2+a3)*100,digits=1),"%",sep=""),paste(round(a1/(a1+a2+a3)*100,digits=1),"%",sep=""))
text(x = xx, y = c(a2,a3,a1)+600, label = label, pos = 3, cex = 0.8)
dev.off()


#### Preparation testing sets

# HGMD

hgmd<-read.table(file=paste(here,"00_scripts/databases-processing/hgmd_2020.2_hg19_missense_with_date_of_entry.tsv",sep=""),sep="\t",skip=1,header=F)

hgmd[,6]=str_split_fixed(as.character(hgmd[,6]), "-", 2)[,1]

# putting OOB score for dbnfsp
a=which(is.element(dbnfsp[,75],train[,75]))
list=dbnfsp[a,75]
res=1:length(a)
for (i in 1:length(a)){
  res[i]=train[which(train[,75]==list[i]),74]
}
dbnfsp[a,74]=res

hgmd[,1]=gsub("chr","",hgmd[,1])
hgmd[,8]=hgmd[,2]+1
hgmd[,9]=apply(hgmd[,c(1,2,3,4)],1,paste,collapse="-")
hgmd[,9]=gsub(" ","",hgmd[,9])
hgmd[,1]=hgmd[,6]
hgmd=hgmd[,-6]

realgno=dbnfsp[which(dbnfsp[,56]!="-6"),]

save(hgmd,realgno, file = paste(here,dataLoc,"hgmd.RData",sep=""))

load(file = paste(here,dataLoc,"hgmd.RData",sep=""))

save(hgmd, file = paste(here,dataLoc,"hgmd-only.RData",sep=""))

listbef=hgmd[,8]

# take DM, not ClinVar and >=2017
hgmd=hgmd[which(hgmd[,5]=="DM"),]
hgmd=hgmd[which(hgmd[,6]=="NULL"),]
hgmd=hgmd[which(hgmd[,1]>=2017),]

list=hgmd[,8]


list3=unique(c(train[,75],gnomad[,75]))
a1=dbnfsp[which(is.element(dbnfsp[,75],list) & is.element(dbnfsp[,75],list3)==F),] # in HGMD and not in gnomad-BLB or train

# remove same AA than training
# temp=train[which(train[,1]=="PLP"),c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAtrainPLP=apply(temp,1,paste,collapse="-",sep="")
# temp=a1[,c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAa1=apply(temp,1,paste,collapse="-",sep="")
# a1=a1[which(is.element(sameAAa1,sameAAtrainPLP)==F),] # in HGMD and not in gnomad-BLB or train

a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a))){a1=a1[-which(a==T),]}
list2=unique(a1[,6])
list3=unique(c(train[,75],VUSCON[,75],gnomad[,75],listbef))

b1=realgno[which(realgno[,56]>(-3.75) & is.element(realgno[,75],list3)==F & is.element(realgno[,6],list2)==T),]

# remove same AA than training
# temp=rbind(train[which(train[,1]=="BLB"),c(2,63)],gnomad[,c(2,63)])
# temp[,2]=as.character(temp[,2])
# sameAAtrainBLB=apply(temp,1,paste,collapse="-",sep="")
# temp=b1[,c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAb1=apply(temp,1,paste,collapse="-",sep="")
# b1=b1[which(is.element(sameAAb1,sameAAtrainBLB)==F),] # in HGMD and not in gnomad-BLB or train

b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b))){b1=b1[-which(b==T),]}

a2=a1[,75]
b2=b1[,75]

save(a1,b1,a2,b2, file=paste(here,dataLoc,"HGMD-testing-review.RData",sep=""))

# DoCM

load(file = paste(here,dataLoc,"hgmd.RData",sep=""))
cos2<-read.table(file=paste(here,"00_scripts/databases-processing/DoCM.pos.tsv",sep=""),sep="\t",header=F)
list=cos2[,1]

list3=unique(c(train[,75],gnomad[,75]))

a1=dbnfsp[which(is.element(dbnfsp[,75],list) & is.element(dbnfsp[,75],list3)==F),]

# remove same AA than training
# temp=train[which(train[,1]=="PLP"),c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAtrainPLP=apply(temp,1,paste,collapse="-",sep="")
# temp=a1[,c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAa1=apply(temp,1,paste,collapse="-",sep="")
# a1=a1[which(is.element(sameAAa1,sameAAtrainPLP)==F),] # in HGMD and not in gnomad-BLB or train


a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a))){a1=a1[-which(a==T),]}
list2=unique(a1[,6])
list3=unique(c(train[,75],gnomad[,75],hgmd[,8],VUSCON[,75]))
b1=realgno[which(realgno[,56]>(-3.75) & is.element(realgno[,75],list3)==F & is.element(realgno[,6],list2) & is.element(realgno[,75],list)==F),]

# remove same AA than training
# temp=rbind(train[which(train[,1]=="BLB"),c(2,63)],gnomad[,c(2,63)])
# temp[,2]=as.character(temp[,2])
# sameAAtrainBLB=apply(temp,1,paste,collapse="-",sep="")
# temp=b1[,c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAb1=apply(temp,1,paste,collapse="-",sep="")
# b1=b1[which(is.element(sameAAb1,sameAAtrainBLB)==F),] # in HGMD and not in gnomad-BLB or train

b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b))){b1=b1[-which(b==T),]}

a2=a1[,75]
b2=b1[,75]
save(a1,b1,a2,b2, file=paste(here,dataLoc,"DoCM-testing-review.RData",sep=""))

# ClinVar

cos2<-read.table(file=paste(here,"testing-set-1/clinvar-20210919-PLP.pos.tsv",sep=""),sep="\t",header=F)
list1=cos2[,1]

cos2<-read.table(file=paste(here,"testing-set-1/clinvar-20210919-BLB.pos.tsv",sep=""),sep="\t",header=F)
list2=cos2[,1]

list3=unique(c(train[,75],gnomad[,75]))

a1=dbnfsp[which(is.element(dbnfsp[,75],list1) & is.element(dbnfsp[,75],list3)==F),]

# remove same AA than training
# temp=train[which(train[,1]=="PLP"),c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAtrainPLP=apply(temp,1,paste,collapse="-",sep="")
# temp=a1[,c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAa1=apply(temp,1,paste,collapse="-",sep="")
# a1=a1[which(is.element(sameAAa1,sameAAtrainPLP)==F),] # in HGMD and not in gnomad-BLB or train

a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a))){a1=a1[-which(a==T),]}
listgene=unique(a1[,6])
#list3=unique(c(train[,75],gnomad[,75],hgmd[,8],VUSCON[,75]))
b1=dbnfsp[which(is.element(dbnfsp[,75],list2) & is.element(dbnfsp[,75],list3)==F & is.element(dbnfsp[,6],listgene)),]

# remove same AA than training
# temp=rbind(train[which(train[,1]=="BLB"),c(2,63)],gnomad[,c(2,63)])
# temp[,2]=as.character(temp[,2])
# sameAAtrainBLB=apply(temp,1,paste,collapse="-",sep="")
# temp=b1[,c(2,63)]
# temp[,2]=as.character(temp[,2])
# sameAAb1=apply(temp,1,paste,collapse="-",sep="")
# b1=b1[which(is.element(sameAAb1,sameAAtrainBLB)==F),] # in HGMD and not in gnomad-BLB or train

b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b))){b1=b1[-which(b==T),]}

a2=a1[,75]
b2=b1[,75]
save(a1,b1,a2,b2, file=paste(here,dataLoc,"ClinVar-testing-review.RData",sep=""))


# MVP de novo

################ de novo vs de novo

cos2<-read.table(file=paste(here,"MVP/autism.tsv",sep=""),sep="\t",header=F)
list1=cos2[,1]

cos2<-read.table(file=paste(here,"MVP/control.tsv",sep=""),sep="\t",header=F)

list2=cos2[,1]

list3=unique(c(train[,75],gnomad[,75]))

a1=dbnfsp[which(is.element(dbnfsp[,75],list1) & is.element(dbnfsp[,75],list3)==F),]

a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a))){a1=a1[-which(a==T),]}
listgene=unique(a1[,6])
#list3=unique(c(train[,75],gnomad[,75],hgmd[,8],VUSCON[,75]))
b1=dbnfsp[which(is.element(dbnfsp[,75],list2) & is.element(dbnfsp[,75],list3)==F & is.element(dbnfsp[,6],listgene)),]

b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b))){b1=b1[-which(b==T),]}

a2=a1[,75]
b2=b1[,75]
save(a1,b1,a2,b2, file=paste(here,dataLoc,"autism-testing-DN.RData",sep=""))

################ de novo2 vs de novo

cos2<-read.table(file=paste(here,"MVP/autism2.tsv",sep=""),sep="\t",header=F)
list1=cos2[,1]

cos2<-read.table(file=paste(here,"MVP/control.tsv",sep=""),sep="\t",header=F)

list2=cos2[,1]

list3=unique(c(train[,75],gnomad[,75]))

a1=dbnfsp[which(is.element(dbnfsp[,75],list1) & is.element(dbnfsp[,75],list3)==F),]

a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a))){a1=a1[-which(a==T),]}
listgene=unique(a1[,6])
#list3=unique(c(train[,75],gnomad[,75],hgmd[,8],VUSCON[,75]))
b1=dbnfsp[which(is.element(dbnfsp[,75],list2) & is.element(dbnfsp[,75],list3)==F & is.element(dbnfsp[,6],listgene)),]

b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b))){b1=b1[-which(b==T),]}

a2=a1[,75]
b2=b1[,75]
save(a1,b1,a2,b2, file=paste(here,dataLoc,"autism-testing-DN2.RData",sep=""))


##############
################ de novo vs clinvar

cos2<-read.table(file=paste(here,"MVP/autism.tsv",sep=""),sep="\t",header=F)
list1=cos2[,1]

cos2<-read.table(file=paste(here,"testing-set-1/clinvar-20210404-BLB.pos.tsv",sep=""),sep="\t",header=F)

list2=cos2[,1]

list3=unique(c(train[,75],gnomad[,75]))

a1=dbnfsp[which(is.element(dbnfsp[,75],list1) & is.element(dbnfsp[,75],list3)==F),]

a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a))){a1=a1[-which(a==T),]}
listgene=unique(a1[,6])
#list3=unique(c(train[,75],gnomad[,75],hgmd[,8],VUSCON[,75]))
b1=dbnfsp[which(is.element(dbnfsp[,75],list2) & is.element(dbnfsp[,75],list3)==F & is.element(dbnfsp[,6],listgene)),]

b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b))){b1=b1[-which(b==T),]}

a2=a1[,75]
b2=b1[,75]
save(a1,b1,a2,b2, file=paste(here,dataLoc,"autism-testing-DC.RData",sep=""))

##############
################ de novo2 vs clinvar

cos2<-read.table(file=paste(here,"MVP/autism2.tsv",sep=""),sep="\t",header=F)
list1=cos2[,1]

cos2<-read.table(file=paste(here,"testing-set-1/clinvar-20210404-BLB.pos.tsv",sep=""),sep="\t",header=F)

list2=cos2[,1]

list3=unique(c(train[,75],gnomad[,75]))

a1=dbnfsp[which(is.element(dbnfsp[,75],list1) & is.element(dbnfsp[,75],list3)==F),]

a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a))){a1=a1[-which(a==T),]}
listgene=unique(a1[,6])
#list3=unique(c(train[,75],gnomad[,75],hgmd[,8],VUSCON[,75]))
b1=dbnfsp[which(is.element(dbnfsp[,75],list2) & is.element(dbnfsp[,75],list3)==F & is.element(dbnfsp[,6],listgene)),]

b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b))){b1=b1[-which(b==T),]}

a2=a1[,75]
b2=b1[,75]
save(a1,b1,a2,b2, file=paste(here,dataLoc,"autism-testing-DC2.RData",sep=""))
##############
################ de novo3 vs clinvar

cos2<-read.table(file=paste(here,"MVP/autism3.tsv",sep=""),sep="\t",header=F)
list1=cos2[,1]

cos2<-read.table(file=paste(here,"testing-set-1/clinvar-20210404-BLB.pos.tsv",sep=""),sep="\t",header=F)

list2=cos2[,1]

list3=unique(c(train[,75],gnomad[,75]))

a1=dbnfsp[which(is.element(dbnfsp[,75],list1) & is.element(dbnfsp[,75],list3)==F),]

a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a))){a1=a1[-which(a==T),]}
listgene=unique(a1[,6])
#list3=unique(c(train[,75],gnomad[,75],hgmd[,8],VUSCON[,75]))
b1=dbnfsp[which(is.element(dbnfsp[,75],list2) & is.element(dbnfsp[,75],list3)==F & is.element(dbnfsp[,6],listgene)),]

b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b))){b1=b1[-which(b==T),]}

a2=a1[,75]
b2=b1[,75]
save(a1,b1,a2,b2, file=paste(here,dataLoc,"autism-testing-DC3.RData",sep=""))

##############


load(file=paste(here,dataLoc,"autism-testing-DN.RData",sep=""))
data=rbind(a1,b1)
colnames(data)[1]="ClinVar"
write.table(data[which(data[,1]=="PLP"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-PLP.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(data[which(data[,1]=="BLB"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-BLB.tsv",sep=""),quote=F,sep="\t",row.names=F)
name="autism-testing-DN"
title="Testing set 4: ASD de novo variants"
plotROC2(data,name,title)

load(file=paste(here,dataLoc,"autism-testing-DN2.RData",sep=""))
data=rbind(a1,b1)
colnames(data)[1]="ClinVar"
write.table(data[which(data[,1]=="PLP"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-PLP.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(data[which(data[,1]=="BLB"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-BLB.tsv",sep=""),quote=F,sep="\t",row.names=F)
name="autism-testing-DN2"
title="Testing set 4: ASD de novo variants"
plotROC2(data,name,title)

load(file=paste(here,dataLoc,"autism-testing-DC.RData",sep=""))
data=rbind(a1,b1)
colnames(data)[1]="ClinVar"
write.table(data[which(data[,1]=="PLP"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-PLP.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(data[which(data[,1]=="BLB"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-BLB.tsv",sep=""),quote=F,sep="\t",row.names=F)
name="autism-testing-DC"
title="Testing set 4: ASD de novo variants"
plotROC2(data,name,title)

load(file=paste(here,dataLoc,"autism-testing-DC2.RData",sep=""))
data=rbind(a1,b1)
colnames(data)[1]="ClinVar"
write.table(data[which(data[,1]=="PLP"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-PLP.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(data[which(data[,1]=="BLB"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-BLB.tsv",sep=""),quote=F,sep="\t",row.names=F)
name="autism-testing-DC2"
title="Testing set 4: ASD de novo variants"
plotROC2(data,name,title)

load(file=paste(here,dataLoc,"autism-testing-DC3.RData",sep=""))
data=rbind(a1,b1)
colnames(data)[1]="ClinVar"
write.table(data[which(data[,1]=="PLP"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-PLP.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(data[which(data[,1]=="BLB"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-ASD-BLB.tsv",sep=""),quote=F,sep="\t",row.names=F)
name="autism-testing-DC3"
title="Testing set 4: ASD de novo variants"
plotROC2(data,name,title)

#######


# exclusion list between testing sets

load(file=paste(here,dataLoc,"HGMD-testing-review.RData",sep=""))
HGMD=c(a2,b2)

load(file=paste(here,dataLoc,"DoCM-testing-review.RData",sep=""))
DOCM=c(a2,b2)

load(file=paste(here,dataLoc,"ClinVar-testing-review.RData",sep=""))
clinvar=c(a2,b2)

removeHGMD=unique(c(intersect(HGMD,DOCM),intersect(HGMD,clinvar)))
removeclinvar=unique(intersect(clinvar,DOCM))

save(removeclinvar,removeHGMD,file = paste(here,dataLoc,"remove-review.RData",sep=""))


########## HGMD testing set ##############################

load(file=paste(here,dataLoc,"HGMD-testing-review.RData",sep=""))

data=rbind(a1,b1)
load(file = paste(here,dataLoc,"remove-review.RData",sep=""))
data=data[-which(is.element(data[,75],removeHGMD)),]

write.table(data[which(data[,1]=="PLP"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-HGMD-PLP.tsv",quote=F,sep="\t",row.names=F)
write.table(data[which(data[,1]=="BLB"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-HGMD-BLB.tsv",quote=F,sep="\t",row.names=F)

colnames(data)[1]="ClinVar"
table(data[,1])

name="HGMD-testing-2017-review"
title="Testing set 2: HGMD disease"

plotROC(data,name,title)

#### Fig. S5B

load(file = paste(here,dataLoc,"hgmd.RData",sep=""))

years=2000:2020

results<-matrix(nrow=42,ncol=24)
results2=1:24

dbnfspa=dbnfsp[which(is.element(dbnfsp[,75],listbef)),]
realgnoa=realgno[which(realgno[,56]>(-3.75)),]

k=3
for (j in years){
  print(j)
  load(file = paste(here,dataLoc,"hgmd-only.RData",sep=""))
  listbef=hgmd[,8]
  hgmd=hgmd[which(hgmd[,5]=="DM"),]
  hgmd=hgmd[which(hgmd[,6]=="NULL"),]
  hgmd=hgmd[which(hgmd[,1]==j),]
  list=hgmd[,8]
  list3=unique(c(train[,75],gnomad[,75]))
  a1=dbnfspa[which(is.element(dbnfspa[,75],list) & is.element(dbnfspa[,75],list3)==F),]
  a1[,1]="PLP"
  a=duplicated(a1[,75])
  if(length(which(a)==T)>0){a1=a1[-which(a==T),]}
  list2=unique(a1[,6])
  list3=unique(c(train[,75],VUSCON[,75],gnomad[,75]))
  b1=realgnoa[which(is.element(realgnoa[,75],listbef)==F & is.element(realgnoa[,75],list3)==F & is.element(realgnoa[,6],list2)==T),]
  b1[,1]="BLB"
  b=duplicated(b1[,75])
  if(length(which(b)==T)>0){b1=b1[-which(b==T),]}
  data=rbind(a1,b1)
  colnames(data)[1]="ClinVar"
  names=c("AA prob","Position","ClinPred","ADA","RF","CONDEL","SIFT","SIFT4G","PolyPhen-HDIV","PolyPhen-HVAR","LRT","MutationTaster",
          "MutationAssessor","FATHMM","PROVEAN","VEST4","metaSVM","metaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL",
          "fathmm-XF","eigen","eigen-PC","GenoCanyon","fitCons","GERP++NR","GERP++RS","PhyloP100way","PhyloP30way","PhyloP17way","phastCons100way","phastCons30way","phastCons17way","SiPhy29way")
  nb=c(69,70,16:55)
  res=cbind(names,names,names)
  res[,2]=nb
  for (i in 1:dim(res)[1]){
    res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
  }
  results[,k]=res[,3]
  results2[k]=auc(roc(data$ClinVar,data[,74]))
  k=k+1
}
results[,1]=names
results2[1]="MutScore"
results=rbind(results,results2)

write.table(results,file=paste(here,dataLoc,"AUC-HGMD-years.tsv",sep=""),quote=F,sep="\t")

pdf(file = paste(here,plot,"HGMD-years.pdf",sep=""),height=10,width=10)
par(mfrow=c(1,1))
plot(years,results[43,3:23],ylim=c(0.8,1),xaxt="n",xlim=c(2000,2020),main="AUC of top-4 predictors for HGMD variants published per year",ylab="AUC",xlab="Year of publication for HGMD pathogenic variants",cex.main=1.8,pch=19,col=6,cex.axis=1.5,cex.lab=2)
axis(1, at=2000:2020,labels=2000:2020, las=2,cex.axis=1.5)
abline(v=2000:2020, col="gray", lty=3)
abline(v=2013.5,col=3,lwd=4)
abline(v=2015.5,col=4,lwd=4)
abline(v=2016.5,col=2,lwd=4)
lines(years,results[43,3:23],col=6)
n=20
points(years,results[n,3:23],pch=19,col=4)
lines(years,results[n,3:23],col=4)
n=22
points(years,results[n,3:23],pch=19,col=3)
lines(years,results[n,3:23],col=3)
n=16
points(years,results[n,3:23],pch=19,col=2)
lines(years,results[n,3:23],col=2)
legend(2003,0.85,legend=c("MutScore","REVEL","MVP","VEST4"),col=c(6,4,3,2),pch=19,lty=1,bg='white',cex=1.5)
dev.off()

pdf(file = paste(here,plot,"HGMD-years2.pdf",sep=""),height=10,width=10)
par(mfrow=c(1,1))

unbiased=3:23
for (i in 3:23){
  unbiased[i]=mean(as.numeric(results[c(7,11,15,36),i]))
}

plot(years,results[43,3:23],ylim=c(0.75,1),xlim=c(2000,2020),main="AUC of top-4 predictors for HGMD variants published per year",ylab="AUC",xlab="Year of publication for HGMD pathogenic variants",cex.main=1.8,pch=19,col=6)
lines(years,results[43,3:23],col=6)
n=20
points(years,results[n,3:23],pch=19,col=4)
lines(years,results[n,3:23],col=4)
n=22
points(years,results[n,3:23],pch=19,col=3)
lines(years,results[n,3:23],col=3)
n=16
points(years,results[n,3:23],pch=19,col=2)
lines(years,results[n,3:23],col=2)

points(years,unbiased[3:23],pch=19,col=5)
lines(years,unbiased[3:23],col=5)

legend(2000,0.79,legend=c("MutScore","REVEL","MVP","VEST4","Unbiased average (SIFT, LRT, PROVEAN and PhyloP100way)"),col=c(6,4,3,2,5),pch=19,lty=1)
abline(v=2016.5)
dev.off()

#### Fig. S5A

load(file = paste(here,dataLoc,"hgmd.RData",sep=""))
hgmd=hgmd[which(hgmd[,5]=="DM"),]
hgmd=hgmd[which(hgmd[,6]=="NULL"),]
hgmd=hgmd[which(hgmd[,1]<2017),]
list=hgmd[,8]
list3=unique(c(train[,75],gnomad[,75]))
a1=dbnfsp[which(is.element(dbnfsp[,75],list) & is.element(dbnfsp[,75],list3)==F),]
a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a)==T)>0){a1=a1[-which(a==T),]}
list2=unique(a1[,6])
list3=unique(c(train[,75],VUSCON[,75],gnomad[,75]))
b1=realgno[which(realgno[,56]>(-3.75) & is.element(realgno[,75],listbef)==F & is.element(realgno[,75],list3)==F & is.element(realgno[,6],list2)==T),]
b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b)==T)>0){b1=b1[-which(b==T),]}
data=rbind(a1,b1)
colnames(data)[1]="ClinVar"
names=c("ClinPred","ADA","RF","CONDEL","SIFT","SIFT4G","PolyPhen-HDIV","PolyPhen-HVAR","LRT","MutationTaster",
        "MutationAssessor","FATHMM","PROVEAN","VEST4","metaSVM","metaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL",
        "fathmm-XF","eigen","eigen-PC","GenoCanyon","fitCons","GERP++NR","GERP++RS","PhyloP100way","PhyloP30way","PhyloP17way","phastCons100way","phastCons30way","phastCons17way","SiPhy29way")
nb=c(16:55)
res=cbind(names,names,names)
res[,2]=nb
for (i in 1:dim(res)[1]){
  res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
}
res1=res
res1=rbind(c("MutScore","X",auc(roc(data$ClinVar,data[,74]))),res1)
write.table(res1[,c(1,3)],file=paste(here,dataLoc,"AUC-testingHGMD-before2017.tsv",sep=""),quote=F,sep="\t")

load(file = paste(here,dataLoc,"hgmd.RData",sep=""))
hgmd=hgmd[which(hgmd[,5]=="DM"),]
hgmd=hgmd[which(hgmd[,6]=="NULL"),]
hgmd=hgmd[which(hgmd[,1]>=2017),]
list=hgmd[,8]
list3=unique(c(train[,75],VUSCON[,75],gnomad[,75]))
a1=dbnfsp[which(is.element(dbnfsp[,75],list) & is.element(dbnfsp[,75],list3)==F),]
a1[,1]="PLP"
a=duplicated(a1[,75])
if(length(which(a)==T)>0){a1=a1[-which(a==T),]}
list2=unique(a1[,6])
b1=realgno[which(realgno[,56]>(-3.75) & is.element(realgno[,75],listbef)==F & is.element(realgno[,75],list3)==F & is.element(realgno[,6],list2)==T),]
b1[,1]="BLB"
b=duplicated(b1[,75])
if(length(which(b)==T)>0){b1=b1[-which(b==T),]}
data=rbind(a1,b1)
colnames(data)[1]="ClinVar"
names=c("ClinPred","ADA","RF","CONDEL","SIFT","SIFT4G","PolyPhen-HDIV","PolyPhen-HVAR","LRT","MutationTaster",
        "MutationAssessor","FATHMM","PROVEAN","VEST4","metaSVM","metaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL",
        "fathmm-XF","eigen","eigen-PC","GenoCanyon","fitCons","GERP++NR","GERP++RS","PhyloP100way","PhyloP30way","PhyloP17way","phastCons100way","phastCons30way","phastCons17way","SiPhy29way")
nb=c(16:55)
res=cbind(names,names,names)
res[,2]=nb
for (i in 1:dim(res)[1]){
  res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
}
res2=res
res2=rbind(c("MutScore","X",auc(roc(data$ClinVar,data[,74]))),res2)

pdf(file = paste(here,plot,"HGMD-beforeafter.pdf",sep=""),height=10,width=10)
par(mfrow=c(1,1))
plot(res1[,3],res2[,3],xlim=c(0.7,0.96),ylim=c(0.7,0.96),xlab="AUC HGMD dataset before 2017",ylab="AUC HGMD testing set in 2017 and after",cex.main=1.8,cex.axis=1.5,cex.lab=2)
text(x = as.numeric(res1[,3]), y = as.numeric(res2[,3])+0.003, label = res1[,1], cex = 1)
abline(a=0,b=1)
dev.off()


########## DoCM testing set ##############

load(file=paste(here,dataLoc,"DoCM-testing-review.RData",sep=""))

data=rbind(a1,b1)
colnames(data)[1]="ClinVar"

write.table(data[which(data[,1]=="PLP"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-DoCM-PLP.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(data[which(data[,1]=="BLB"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-DoCM-BLB.tsv",sep=""),quote=F,sep="\t",row.names=F)

name="DoCM-testing-review"
title="Testing set 3: DoCM somatic cancer"

plotROC(data,name,title)

table(data[,1])
length(unique(data[,6]))
write.table(res1[,c(1,3)],file=paste(here,dataLoc,"AUC-testing3-docm.tsv",sep=""),quote=F,sep="\t")


########## ClinVar testing set ###########################################################################################################

load(file=paste(here,dataLoc,"ClinVar-testing-review.RData",sep=""))

data=rbind(a1,b1)
load(file = paste(here,dataLoc,"remove-review.RData",sep=""))
if(length(removeclinvar)>0){
	data=data[-which(is.element(data[,75],removeclinvar)),]
}
colnames(data)[1]="ClinVar"

write.table(data[which(data[,1]=="PLP"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-clinvar-PLP-review.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(data[which(data[,1]=="BLB"),c(8,9,11,12)],file=paste(here,dataLoc,"sets-clinvar-BLB-review.tsv",sep=""),quote=F,sep="\t",row.names=F)

name="ClinVar-testing-review"
title="Testing set 1: Recent ClinVar"

plotROC(data,name,title)

table(data[,1])
length(unique(data[,6]))
write.table(res1[,c(1,3)],file=paste(here,dataLoc,"AUC-testing1-ClinVar.tsv",sep=""),quote=F,sep="\t")
  








test<-read.table(file=paste(here,"data-20210919/","testing.tsv",sep=""),header=FALSE,sep="\t")  # OK

load(file=paste(here,dataLoc,"ClinVar-testing-review.RData",sep=""))

data=rbind(a1,b1)
load(file = paste(here,dataLoc,"remove-review.RData",sep=""))
if(length(removeclinvar)>0){
	data=data[-which(is.element(data[,75],removeclinvar)),]
}
colnames(data)[1]="ClinVar"

dataALL=data

### Training de novo

l1=test[which(test[,1]=="PLP" & (test[,4]==32 | test[,4]==33)),2]
p1=dataALL[which(is.element(dataALL[,75],l1)),]
l2=test[which(test[,1]=="BLB"),2]
p2=dataALL[which(is.element(dataALL[,75],l2) & is.element(dataALL[,6],unique(p1[,6]))),]
data=rbind(p1,p2)

name="ClinVar-testing-denovo-review"
title="De novo PLP and all BLB ClinVar"

plotROC(data,name,title)


### Training germline

l1=test[which(test[,1]=="PLP" & (test[,4]==1)),2]
p1=dataALL[which(is.element(dataALL[,75],l1)),]
l2=test[which(test[,1]=="BLB"),2]
p2=dataALL[which(is.element(dataALL[,75],l2) & is.element(dataALL[,6],unique(p1[,6]))),]
data=rbind(p1,p2)

name="ClinVar-testing-germline-review"
title="Germline PLP and all BLB ClinVar"

plotROC(data,name,title)

### 0 stars

l1=test[which(test[,1]=="PLP" & (test[,3]=="no_assertion_criteria_provided")),2]
p1=dataALL[which(is.element(dataALL[,75],l1)),]
l2=test[which(test[,1]=="BLB"),2]
p2=dataALL[which(is.element(dataALL[,75],l2)),]
data=rbind(p1,p2)
name="ClinVar-testing-0star-review"
title="0 star PLP and all BLB ClinVar"
plotROC(data,name,title)

### 1 stars

l1=test[which(test[,1]=="PLP" & (test[,3]=="criteria_provided,_single_submitter")),2]
p1=dataALL[which(is.element(dataALL[,75],l1)),]
l2=test[which(test[,1]=="BLB"),2]
p2=dataALL[which(is.element(dataALL[,75],l2)),]
data=rbind(p1,p2)
name="ClinVar-testing-1star-review"
title="1 star PLP and all BLB ClinVar"
plotROC(data,name,title)

### 2 stars

l1=test[which(test[,1]=="PLP" & (test[,3]=="criteria_provided,_multiple_submitters,_no_conflicts")),2]
p1=dataALL[which(is.element(dataALL[,75],l1)),]
l2=test[which(test[,1]=="BLB"),2]
p2=dataALL[which(is.element(dataALL[,75],l2)),]
data=rbind(p1,p2)
name="ClinVar-testing-2star-review"
title="2 stars PLP and all BLB ClinVar"
plotROC(data,name,title)

### 3 stars

l1=test[which(test[,1]=="PLP" & (test[,3]=="reviewed_by_expert_panel")),2]
p1=dataALL[which(is.element(dataALL[,75],l1)),]
l2=test[which(test[,1]=="BLB"),2]
p2=dataALL[which(is.element(dataALL[,75],l2)),]
data=rbind(p1,p2)
name="ClinVar-testing-3star-review"
title="3 stars PLP and all BLB ClinVar"
plotROC(data,name,title)





















### scores for all ###########

# train,VUSCON,gnomad,dbnfsp
load(file = paste(here,dataLoc,"data10.RData",sep=""))

wrong=c("FPGT-TNNI3K","UGT1A5","UGT1A8","UGT1A6","UGT1A7","UGT1A3","UGT1A9","UGT1A4","UGT1A10","ABHD14A-ACY1","SMN2","KAAG1","CNPY3-GNMT","PGBD3","INS-IGF2","NDUFC2-KCTD14","FXYD6-FXYD2","BIVM-ERCC5","ST20-MTHFS","CORO7-PAM16","CORO7-PAM16","KCNE1B","CBSL","U2AF1L5","CRYAA2","SIK1B","LOC102724788","OPN1MW3","OPN1MW2")
VUSCON=VUSCON[-which(is.element(VUSCON[,6],wrong)),]
gnomad=gnomad[-which(is.element(gnomad[,6],wrong)),]

list1=unique(train[,75])
list2=unique(VUSCON[,75])
list3=unique(gnomad[,75])
list4=unique(dbnfsp[,75])

sel2=setdiff(list2,list1)
sel3=setdiff(list3,c(sel2,list1))
sel4=setdiff(list4,c(sel3,sel2,list1))

VUSCON=VUSCON[which(is.element(VUSCON[,75],sel2)),]
gnomad=gnomad[which(is.element(gnomad[,75],sel3)),]
dbnfsp=dbnfsp[which(is.element(dbnfsp[,75],sel4)),]

DF <- as.matrix(data.frame(table(VUSCON[,75])))
b=which(as.numeric(DF[,2])>1)
VUSCONtemp=VUSCON[which(is.element(VUSCON[,75],DF[b,1])),]
VUSCON=VUSCON[-which(is.element(VUSCON[,75],DF[b,1])),]
new=c()
for (i in 1:length(b)){
  a=which(VUSCONtemp[,75]==DF[b[i],1])
  dat=VUSCONtemp[a[1],]
  dat[,70]=mean(as.numeric(VUSCONtemp[a,70]))
  dat[,74]=mean(as.numeric(VUSCONtemp[a,74]))
  new=rbind(new,dat)
}
VUSCON=rbind(VUSCON,new)

DF <- as.matrix(data.frame(table(gnomad[,75])))
b=which(as.numeric(DF[,2])>1)
gnomadtemp=gnomad[which(is.element(gnomad[,75],DF[b,1])),]
gnomad=gnomad[-which(is.element(gnomad[,75],DF[b,1])),]
new=c()
for (i in 1:length(b)){
  a=which(gnomadtemp[,75]==DF[b[i],1])
  dat=gnomadtemp[a[1],]
  dat[,70]=mean(as.numeric(gnomadtemp[a,70]))
  dat[,74]=mean(as.numeric(gnomadtemp[a,74]))
  new=rbind(new,dat)
}
gnomad=rbind(gnomad,new)

DF <- as.matrix(data.frame(table(dbnfsp[,75])))
b=which(as.numeric(DF[,2])>1)
dbnfsptemp=dbnfsp[which(is.element(dbnfsp[,75],DF[b,1])),]
dbnfsp=dbnfsp[-which(is.element(dbnfsp[,75],DF[b,1])),]
new=c()
for (i in 1:length(b)){
  print(i)
  a=which(dbnfsptemp[,75]==DF[b[i],1])
  dat=dbnfsptemp[a[1],]
  dat[,70]=mean(as.numeric(dbnfsptemp[a,70]))
  dat[,74]=mean(as.numeric(dbnfsptemp[a,74]))
  new=rbind(new,dat)
}
colnames(dbnfsp)[1]="V1"
dbnfsp=rbind(dbnfsp,new)

colnames(train)[1]="V1"
data=rbind(train,VUSCON,gnomad,dbnfsp)

save(data,file = paste(here,dataLoc,"data11.RData",sep=""))

load(file = paste(here,dataLoc,"data11.RData",sep=""))

mutscore=data[,c(8:12,69,70,74)]

save(mutscore,file = paste(here,dataLoc,"data11-mutscore.RData",sep=""))

load(file = paste(here,dataLoc,"data11-mutscore.RData",sep=""))

mutscore[,6]=round(mutscore[,6],digit=3)
mutscore[,7]=round(mutscore[,7],digit=3)
mutscore[,8]=round(mutscore[,8],digit=3)
colnames(mutscore)=c("Chr","Pos","Pos2","Ref","Alt","AA-score","positional-score","MutScore")
mutscore2=mutscore[,c(1,2,4,5,8)]
mutscore2=mutscore2[order(mutscore2[,1],mutscore2[,2],mutscore2[,3],mutscore2[,4]),]
rownames(mutscore2) <- NULL
save(mutscore2,file = paste(here,dataLoc,"data12-mutscore.RData",sep=""))
write.table(mutscore2,file=paste(here,dataLoc,"mutscore-hg19.tsv",sep=""),quote=F,sep="\t",row.names=F,col.names=F)

load(file = paste(here,dataLoc,"data12-mutscore.RData",sep=""))

# creating database for ANNOVAR

mutscore[,3]=mutscore[,2]
colnames(mutscore)=c("#Chr","Start","End","Ref","Alt","AA-score","positional-score","MutScore")
mutscore=mutscore[order(mutscore[,1],mutscore[,2]),]
write.table(mutscore,file=paste(here,"00_scripts/clinvar_annotation/humandb/hg19_mutscore.txt",sep=""),quote=F,sep="\t",row.names=F)

# write different sets of variants

trainBLB=train[which(train[,1]=="BLB"),c(8,9,11,12)]
trainPLP=train[which(train[,1]=="PLP"),c(8,9,11,12,60,62)]

trainPLP[which(trainPLP[,5]=="no_assertion_criteria_provided"),5]="0-star"
trainPLP[which(trainPLP[,5]=="criteria_provided,_single_submitter"),5]="1-star"
trainPLP[which(trainPLP[,5]=="criteria_provided,_multiple_submitters,_no_conflicts"),5]="2-star"
trainPLP[which(trainPLP[,5]=="reviewed_by_expert_panel"),5]="3-star"
trainPLP[which(trainPLP[,5]=="practice_guideline"),5]="4-star"

trainPLP[which(trainPLP[,6]==1),6]="germline"
trainPLP[which(trainPLP[,6]==2),6]="somatic"
trainPLP[which(trainPLP[,6]==32 | trainPLP[,6]==33),6]="de-novo"
trainPLP[which(trainPLP[,6]!="germline" & trainPLP[,6]!="somatic" & trainPLP[,6]!="de-novo" ),6]="other"

gnomadBLB=gnomad[,c(8,9,11,12)]

list=unique(gnomad[,73])
genes=unique(train[which(train[,1]=="PLP"),6])
VUS=VUSCON[which(VUSCON[,1]=="VUS" & is.element(VUSCON[,6],genes)),]
VUS=VUS[-which(is.element(VUS[,73],list)),]
CON=VUSCON[which(VUSCON[,1]=="CON" & is.element(VUSCON[,6],genes)),]
CON=CON[-which(is.element(CON[,73],list)),]
VUSBLB=VUS[,c(8,9,11,12)]
CONBLB=CON[,c(8,9,11,12)]

write.table(trainBLB,file=paste(here,dataLoc,"sets-training-BLB.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(trainPLP,file=paste(here,dataLoc,"sets-training-PLP.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(gnomadBLB,file=paste(here,dataLoc,"sets-gnomad-BLB.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(VUSBLB,file=paste(here,dataLoc,"sets-VUS.tsv",sep=""),quote=F,sep="\t",row.names=F)
write.table(CONBLB,file=paste(here,dataLoc,"sets-CON.tsv",sep=""),quote=F,sep="\t",row.names=F)


#### colagen 

load(file = paste(here,dataLoc,"data11.RData",sep=""))

#COL3A1
iso="NM_000090.3"
l1=154
l2=1221

cola=data[which(data[,2]==iso),]
cola=cola[,c(1,2,3,4,14,74)]

cola[,2]=substring(cola[,5], 1, 1)
cola[,3]=substring(cola[,4], 4, 100)
cola[,3]=substring(cola[,3], 1, nchar(cola[,3])-1)

cola[,3]=as.numeric(cola[,3])

a=1:2000
x1=a*3
x2=a*3-1
x3=a*3-2

Gin=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2)
length(intersect(cola[Gin,3],x1))
length(intersect(cola[Gin,3],x2))
length(intersect(cola[Gin,3],x3))

Gin=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2 & is.element(cola[,3],x1))
Ginno=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2 & is.element(cola[,3],x1)==F)
Gout=which(cola[,2]=="G" & (cola[,3]<l1 | cola[,3]>l2))

Xin=which(cola[,2]!="G" & cola[,3]>=l1 & cola[,3]<=l2)
Xout=which(cola[,2]!="G" & (cola[,3]<l1 | cola[,3]>l2))

dim(cola)
length(Gin)+length(Ginno)+length(Gout)+length(Xin)+length(Xout)

#boxplot(cola[Gin,6],cola[Ginno,6],cola[Gout,6],cola[Xin,6],cola[Xout,6])

par(mfrow=c(5,1))
hist(cola[Gin,6],breaks=seq(0,1,0.025))
hist(cola[Ginno,6],breaks=seq(0,1,0.025))
hist(cola[Gout,6],breaks=seq(0,1,0.025))
hist(cola[Xin,6],breaks=seq(0,1,0.025))
hist(cola[Xout,6],breaks=seq(0,1,0.025))

pdf(file=paste(here,plot,"COL3A1.pdf",sep=""))
par(mfrow=c(3,1))
hist(cola[Gin,6],breaks=seq(0,1,0.01),main="COL3A1: Glycine changes in the triple helix domain (one each 3aa)",xlab="MutScore")
hist(cola[Ginno,6],breaks=seq(0,1,0.01),main="COL3A1: Glycine changes in the triple helix domain (not one of each 3)",xlab="MutScore")
hist(cola[Xin,6],breaks=seq(0,1,0.01),main="COL3A1: Other amino acid changes in the triple helix domain",xlab="MutScore")
dev.off()

### COL1A1

iso="NM_000088.3"
l1=162
l2=1218

cola=data[which(data[,2]==iso),]
cola=cola[,c(1,2,3,4,14,74)]

cola[,2]=substring(cola[,5], 1, 1)
cola[,3]=substring(cola[,4], 4, 100)
cola[,3]=substring(cola[,3], 1, nchar(cola[,3])-1)

cola[,3]=as.numeric(cola[,3])

a=1:2000
x1=a*3
x2=a*3-1
x3=a*3-2

Gin=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2)
length(intersect(cola[Gin,3],x1))
length(intersect(cola[Gin,3],x2))
length(intersect(cola[Gin,3],x3))

Gin=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2 & is.element(cola[,3],x2))
Ginno=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2 & is.element(cola[,3],x2)==F)
Gout=which(cola[,2]=="G" & (cola[,3]<l1 | cola[,3]>l2))

Xin=which(cola[,2]!="G" & cola[,3]>=l1 & cola[,3]<=l2)
Xout=which(cola[,2]!="G" & (cola[,3]<l1 | cola[,3]>l2))

dim(cola)
length(Gin)+length(Ginno)+length(Gout)+length(Xin)+length(Xout)

#boxplot(cola[Gin,6],cola[Ginno,6],cola[Gout,6],cola[Xin,6],cola[Xout,6])

par(mfrow=c(5,1))
hist(cola[Gin,6],breaks=seq(0,1,0.025))
hist(cola[Ginno,6],breaks=seq(0,1,0.025))
hist(cola[Gout,6],breaks=seq(0,1,0.025))
hist(cola[Xin,6],breaks=seq(0,1,0.025))
hist(cola[Xout,6],breaks=seq(0,1,0.025))

pdf(file=paste(here,plot,"COL1A1.pdf",sep=""))
par(mfrow=c(3,1))
hist(cola[Gin,6],breaks=seq(0,1,0.01),main="COL1A1: Glycine changes in the triple helix domain (one each 3aa)",xlab="MutScore")
hist(cola[Ginno,6],breaks=seq(0,1,0.01),main="COL1A1: Glycine changes in the triple helix domain (in between triplets)",xlab="MutScore")
hist(cola[Xin,6],breaks=seq(0,1,0.01),main="COL1A1: Other amino acid changes in the triple helix domain",xlab="MutScore")
dev.off()

### COL1A2

iso="NM_000089.3"
l1=80
l2=1119

cola=data[which(data[,2]==iso),]
cola=cola[,c(1,2,3,4,14,74)]

cola[,2]=substring(cola[,5], 1, 1)
cola[,3]=substring(cola[,4], 4, 100)
cola[,3]=substring(cola[,3], 1, nchar(cola[,3])-1)

cola[,3]=as.numeric(cola[,3])

a=1:2000
x1=a*3
x2=a*3-1
x3=a*3-2

Gin=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2)
length(intersect(cola[Gin,3],x1))
length(intersect(cola[Gin,3],x2))
length(intersect(cola[Gin,3],x3))

Gin=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2 & is.element(cola[,3],x3))
Ginno=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2 & is.element(cola[,3],x3)==F)
Gout=which(cola[,2]=="G" & (cola[,3]<l1 | cola[,3]>l2))

Xin=which(cola[,2]!="G" & cola[,3]>=l1 & cola[,3]<=l2)
Xout=which(cola[,2]!="G" & (cola[,3]<l1 | cola[,3]>l2))

dim(cola)
length(Gin)+length(Ginno)+length(Gout)+length(Xin)+length(Xout)

#boxplot(cola[Gin,6],cola[Ginno,6],cola[Gout,6],cola[Xin,6],cola[Xout,6])

par(mfrow=c(5,1))
hist(cola[Gin,6],breaks=seq(0,1,0.025))
hist(cola[Ginno,6],breaks=seq(0,1,0.025))
hist(cola[Gout,6],breaks=seq(0,1,0.025))
hist(cola[Xin,6],breaks=seq(0,1,0.025))
hist(cola[Xout,6],breaks=seq(0,1,0.025))

pdf(file=paste(here,plot,"COL1A2.pdf",sep=""))
par(mfrow=c(3,1))
hist(cola[Gin,6],breaks=seq(0,1,0.01),main="COL1A2: Glycine changes in the triple helix domain (one each 3aa)",xlab="MutScore")
hist(cola[Ginno,6],breaks=seq(0,1,0.01),main="COL1A2: Glycine changes in the triple helix domain (in between triplets)",xlab="MutScore")
hist(cola[Xin,6],breaks=seq(0,1,0.01),main="COL1A2: Other amino acid changes in the triple helix domain",xlab="MutScore")
dev.off()

### COL2A1

iso="NM_001844.5"
l1=80
l2=1119

cola=data[which(data[,2]==iso),]
cola=cola[,c(1,2,3,4,14,74)]

cola[,2]=substring(cola[,5], 1, 1)
cola[,3]=substring(cola[,4], 4, 100)
cola[,3]=substring(cola[,3], 1, nchar(cola[,3])-1)

cola[,3]=as.numeric(cola[,3])

a=1:2000
x1=a*3
x2=a*3-1
x3=a*3-2

Gin=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2)
length(intersect(cola[Gin,3],x1))
length(intersect(cola[Gin,3],x2))
length(intersect(cola[Gin,3],x3))

Gin=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2 & is.element(cola[,3],x1))
Ginno=which(cola[,2]=="G" & cola[,3]>=l1 & cola[,3]<=l2 & is.element(cola[,3],x1)==F)
Gout=which(cola[,2]=="G" & (cola[,3]<l1 | cola[,3]>l2))

Xin=which(cola[,2]!="G" & cola[,3]>=l1 & cola[,3]<=l2)
Xout=which(cola[,2]!="G" & (cola[,3]<l1 | cola[,3]>l2))

dim(cola)
length(Gin)+length(Ginno)+length(Gout)+length(Xin)+length(Xout)

#boxplot(cola[Gin,6],cola[Ginno,6],cola[Gout,6],cola[Xin,6],cola[Xout,6])

par(mfrow=c(5,1))
hist(cola[Gin,6],breaks=seq(0,1,0.025))
hist(cola[Ginno,6],breaks=seq(0,1,0.025))
hist(cola[Gout,6],breaks=seq(0,1,0.025))
hist(cola[Xin,6],breaks=seq(0,1,0.025))
hist(cola[Xout,6],breaks=seq(0,1,0.025))

pdf(file=paste(here,plot,"COL2A1.pdf",sep=""))
par(mfrow=c(3,1))
hist(cola[Gin,6],breaks=seq(0,1,0.01),main="COL2A1: Glycine changes in the triple helix domain (one each 3aa)",xlab="MutScore")
hist(cola[Ginno,6],breaks=seq(0,1,0.01),main="COL2A1: Glycine changes in the triple helix domain (in between triplets)",xlab="MutScore")
hist(cola[Xin,6],breaks=seq(0,1,0.01),main="COL2A1: Other amino acid changes in the triple helix domain",xlab="MutScore")
dev.off()


