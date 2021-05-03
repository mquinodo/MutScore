
here="/home/mquinodo/SYNO/WES/Clinvar2/"

refseq<-read.table(paste(here,"Shiny/refse.txt",sep=""),header=TRUE,sep="\t",fill=TRUE)
uniprot<-read.table(paste(here,"Shiny/uniprot.txt",sep=""),header=TRUE,sep="\t",fill=TRUE)
conv=read.table(file=paste(here,"Shiny/ENST-NM.pro.txt",sep=""),skip=1,header=F)
conv[,5]=""
conv=conv[-which(is.na(conv[,3]) & is.na(conv[,4])),]

save(refseq,uniprot,conv,file=paste(here,"Shiny/preprocessed.RData",sep=""))


load(file = paste(here,"Shiny/data-all.RData",sep=""))
load(file = paste(here,"Shiny/scores-all.RData",sep=""))

maxpos=unique(rbind(data[,c(2,18)],scores[,c(1,8)]))

isos=unique(maxpos[,1])

maxiso=cbind(isos,isos)
for (i in 1:length(isos)){
	print(i)
	maxiso[i,2]=max(maxpos[which(maxpos[,1]==isos[i]),2],na.rm=T)
}
maxiso[,2]=as.numeric(maxiso[,2])

save(maxiso,file=paste(here,"Shiny/maxiso-all.RData",sep=""))

#################

load(file=paste(here,"Shiny/preprocessed.RData",sep=""))
load(file =paste(here,"Shiny/regions.RData",sep=""))
load(file=paste(here,"Shiny/data-all.RData",sep=""))
load(file = paste(here,"Shiny/maxiso.RData",sep=""))
load(file=paste(here,"Shiny/scores-all.RData",sep=""))

load(file=paste(here,"Shiny/data-all.RData",sep=""))

load(file=paste(here,"Shiny/other.RData",sep=""))

#genes=sort(as.character(unique(data[which(data[,2]!="." & data[,1]=="PLP"),5])))
genes=sort(as.character(unique(data[which(data[,2]!="."),5])))

dataallgenes=data
data1=dataallgenes[which(dataallgenes[,4]=="nonsynonymous SNV"),c(5,2,1,18,3,15)]
data1[which(data1[,6]=="-7"),3]="gnomAD"
data1=data1[,1:5]

scores[,2]=as.numeric(unlist(scores[,2]))
scores[,3]=as.numeric(unlist(scores[,3]))
scores[,4]=as.numeric(unlist(scores[,4]))
scores[,5]=as.numeric(unlist(scores[,5]))
scores[,6]=as.numeric(unlist(scores[,6]))
scores[,7]=as.numeric(unlist(scores[,7]))


iso=unique(c(other3[,1],scores[,1],dataallgenes[,2],data1[,2]))
iso=cbind(iso,iso)
iso[,2]=sample(1:500, dim(iso)[1], replace=T)

scoresT=scores
data1T=data1
dataallgenesT=dataallgenes

for(i in 1:500){
	print(i)
  sel=iso[which(iso[,2]==i),1]
  mut1=other3[which(is.element(other3[,1],sel)),]
  scores=scoresT[which(is.element(scoresT[,1],sel)),]
  data1=data1T[which(is.element(data1T[,2],sel)),]
  dataallgenes=dataallgenesT[which(is.element(dataallgenesT[,2],sel)),]
  save(mut1,scores,data1,dataallgenes, file=paste(here,"Shiny/cut/data-",i,".RData",sep=""))
}

iso[which(iso[,1]=="NM_022455.4"),]
# 111

isosel=iso
save(isosel, file=paste(here,"Shiny/cut/isosel.RData",sep=""))

geneinfo=dataallgenesT[,c(2,5)]

save(regionsb,regionsp,maxiso,conv,refseq,uniprot,genes,regionsp1,regionsb1,resip1,resib1,resi1,resi2,geneinfo,file=paste(here,"Shiny/data-mutland-small-all.RData",sep=""))

#save(scores,regionsb,regionsp,maxiso,conv,refseq,uniprot,genes,dataallgenes,data1,regionsp1,regionsb1,resip1,resib1,resi1,resi2,file=paste(here,"Shiny/data-mutland2.RData",sep=""))

save(genes,maxiso, file=paste(here,"Shiny/init.RData",sep=""))

