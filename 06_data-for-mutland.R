
library(MASS)
library(data.table)

### TO BE CHANGED ACCORDINGLY
here="/home/mquinodo/SYNO/WES/Clinvar2/"
version="20201121"


dataLoc=paste("data-",version,"/",sep="")
plot=paste("plots-",version,"/",sep="")

load(file = paste(here,dataLoc,"data2.RData",sep=""))

# shortening of names in column 1
data[which(data[,1]==paste("clinvar-",version,"-PLP",sep="")),1]="PLP"
data[which(data[,1]==paste("clinvar-",version,"-BLB",sep="")),1]="BLB"
data[which(data[,1]==paste("clinvar-",version,"-VUS",sep="")),1]="VUS"
data[which(data[,1]==paste("clinvar-",version,"-CON",sep="")),1]="CON"

PLPgene=unique(data[which(data[,1]=="PLP"),6])
#data=data[which(is.element(data[,6],PLPgene)),]
#data2=data2[which(is.element(data2[,6],PLPgene)),]

db=which(data[,1]=="dbNFSP4.0")
data=data[-db,]

genelist=unique(data[,6])
data[,64]=0
data[,65]=0
for (i in 1:length(genelist)){
  print(i)
  PLP=which(data[,6]==genelist[i] & data[,1]=="PLP")
  ALL=which(data[,6]==genelist[i])
  maxPLP=max(data[PLP,56])
  nb=length(unique(data[PLP,8:12])[,1])
  data[ALL,64]=nb
  data[ALL,65]=maxPLP
  if(nb>9){
    data[ALL[which(data[ALL,1]=="gnomAD-all" & data[ALL,56]>maxPLP)],1]="OK"
  }
}
data[which(data[,1]=="OK"),57]=(-7)
data[which(data[,1]=="OK"),1]="BLB"
data=data[-which(data[,1]=="gnomAD-all"),]

data=data[,c(1,2,4:9,11,12,48,49,52,56,57,58,60,63,65)]
colnames(data)=c("cat","iso","prot","type","gene","exon","chr","begin","ref","alt","GERP","phylo","phast","gno1","gno2","pext","stars","pos","maxaf")

data[,c(6:13,16)]=0

save(data, file=paste(here,"Shiny/data-all.RData",sep=""))

#####

data2<-read.table(file=paste(here,dataLoc,"ALL-",version,"-dbNFSP4.0.tsv",sep=""),header=FALSE,sep="\t")  # OK


data2=data2[which(data2[,5]=="nonsynonymous SNV"),]
# transforming protein notation to position and store in column 63
data2[,66]=0
temp=gsub("f","",gsub(".{1}$","",substring(gsub("p.","",data2[,4]),2)))
temp2=temp
temp3=temp
for (i in 1:length(temp)){
  temp2[i]=strsplit(temp[i],"_")[[1]][1]
  temp3[i]=strsplit(temp2[i],"del")[[1]][1]
}
data2[,66]=as.numeric(temp3)

data2[,7]=substring(data2[,7], 5)

# transforming protein notation to position and store in column 63
data2[,67]=0
temp=substring(gsub("c.","",data2[,3]),1,nchar(data2[,3])-5)
data2[,67]=as.numeric(temp)

other=data2[,c(2,7,48,49,52,58,66,62,3,67,4,8,9,11,12,56,60,61)]
colnames(other)=c("iso","exon","GERP","phylo","phast","pext","pos","mutscore","c.","c-pos","p.","chr","begin","ref","alt","gnomAD","aa-score","posi-score")
rownames(other) <- c()

other[,4]=round(as.numeric(other[,4]),digits=3)
other[,5]=round(as.numeric(other[,5]),digits=3)
other[,6]=round(as.numeric(other[,6]),digits=3)
other[,8]=round(as.numeric(other[,8]),digits=3)
other[,2]=as.numeric(other[,2])
other[,3]=as.numeric(other[,3])
other[,16]=as.numeric(other[,16])
other[,17]=as.numeric(other[,17])
other[,18]=as.numeric(other[,18])

other2=other[,c(1,9,10,11,7,8)]
other3=other2[-which(is.na(other2[,6])),]

other3=other3[which(other3[,5]!=1),]

save(other3, file=paste(here,"Shiny/other.RData",sep=""))

other2=unique(other[,1:8])
scores=other2

keys <- c("iso","exon","pos")

X <- as.data.table(scores)
X2=X[,lapply(.SD,mean),keys]

scores2=X2[,c(1,2,4,5,6,7,8,3)]

scores2[,3]=round(scores2[,3],2)
scores2[,4]=round(scores2[,4],2)
scores2[,5]=round(scores2[,5],2)
scores2[,6]=round(scores2[,6],2)
scores2[,7]=round(scores2[,7],2)

scores=scores2

save(scores, file=paste(here,"Shiny/scores-all.RData",sep=""))


################# clean regions for mutland

load(file = paste(here,dataLoc,"regions.RData",sep=""))

write.table(resi1,file=paste(here,dataLoc,"cluster-PLP.tsv",sep=""),quote=F,sep="\t")
write.table(resi2,file=paste(here,dataLoc,"cluster-BLB.tsv",sep=""),quote=F,sep="\t")

p=unique(resi1[which(as.numeric(resi1[,7])<0.05),1])
b=unique(resi2[which(as.numeric(resi2[,7])<0.05),1])

regionsp=regionsp[which(is.element(regionsp[,2],p)),]
regionsb=regionsb[which(is.element(regionsb[,2],b)),]

# regionsp
# 1=gene / 2=isoform / 3=begin-region / 4=end-region / 5=chr / 6=begin-pos / 7=end-pos / 8=nb-PLP / 9=nb-BLB / 10=%-PLP-inside / 11=%-region-gene / 12=nb-PLP-tot / 13=size-gene / 14=size-region
# resi
# 1=isoform / 2=gene / 3=%PLP-regionS / 4=%gene-regionS / 5=score / 6=Fisher-estimate / 7=Fisher-p.value

regionsp1=regionsp[,c(1,2,3,4,8,9)]
regionsb1=regionsb[,c(1,2,3,4,8,9)]

resip1=resi1[,c(2,1,3,4,5,6,7)]
resip1[,3]=round(as.numeric(resip1[,3]),digits=2)
resip1[,4]=round(as.numeric(resip1[,4]),digits=2)
resip1[,5]=round(as.numeric(resip1[,5]),digits=2)
resip1[,6]=round(as.numeric(resip1[,6]),digits=2)
resip1[,7]=signif(as.numeric(resip1[,7]),digits=3)

resib1=resi2[,c(2,1,3,4,5,6,7)]
resib1[,3]=round(as.numeric(resib1[,3]),digits=2)
resib1[,4]=round(as.numeric(resib1[,4]),digits=2)
resib1[,5]=round(as.numeric(resib1[,5]),digits=2)
resib1[,6]=round(as.numeric(resib1[,6]),digits=2)
resib1[,7]=signif(as.numeric(resib1[,7]),digits=3)

resi1=resi1[,c(2,1,3,4,5,6,7,8)]
resi1[,3]=round(as.numeric(resi1[,3]),digits=2)
resi1[,4]=round(as.numeric(resi1[,4]),digits=2)
resi1[,5]=round(as.numeric(resi1[,5]),digits=2)
resi1[,6]=round(as.numeric(resi1[,6]),digits=2)
resi1[,7]=signif(as.numeric(resi1[,7]),digits=3)
resi1[,8]=signif(as.numeric(resi1[,8]),digits=3)

resi2=resi2[,c(2,1,3,4,5,6,7,8)]
resi2[,3]=round(as.numeric(resi2[,3]),digits=2)
resi2[,4]=round(as.numeric(resi2[,4]),digits=2)
resi2[,5]=round(as.numeric(resi2[,5]),digits=2)
resi2[,6]=round(as.numeric(resi2[,6]),digits=2)
resi2[,7]=signif(as.numeric(resi2[,7]),digits=3)
resi2[,8]=signif(as.numeric(resi2[,8]),digits=3)

save(regionsp,regionsb,resip1,resib1,regionsp1,regionsb1,resi1,resi2,file=paste(here,"Shiny/regions.RData",sep=""))

#### hg38

mutscore<-read.table(file=paste(here,dataLoc,"mutscore-hg38.tsv",sep=""),header=FALSE,sep="\t")  # OK

colnames(mutscore)=c("Chr","Pos","Ref","Alt","MutScore")
mutscore=mutscore[order(mutscore[,1],mutscore[,2],mutscore[,3],mutscore[,4]),]
rownames(mutscore) <- NULL
save(mutscore,file = paste(here,dataLoc,"data12-mutscore-hg38.RData",sep=""))


