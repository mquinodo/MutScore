
library(gridExtra)
library(MASS)
library(pROC)
library(dplyr)
library(randomForest)
require(caTools)
library(rpart)
library(rpart.plot)

### TO BE CHANGED ACCORDINGLY
here="/home/mquinodo/SYNO/WES/Clinvar2/"
dataLoc="data-20201121/"
plot="plots-20201121/"
###

load(file=paste(here,dataLoc,"data4.RData",sep=""))

# transform protein change notation to amino acid change (p.R154G -> RG)
temp=gsub("p.","",data[,4])
temp2=gsub('[[:digit:]]+', '', temp)
data[,14]=temp2

data[,67:74]=0

wrong=c("FPGT-TNNI3K","UGT1A5","UGT1A8","UGT1A6","UGT1A7","UGT1A3","UGT1A9","UGT1A4","UGT1A10","ABHD14A-ACY1","SMN2","KAAG1","CNPY3-GNMT",
  "PGBD3","INS-IGF2","NDUFC2-KCTD14","FXYD6-FXYD2","BIVM-ERCC5","ST20-MTHFS","CORO7-PAM16","CORO7-PAM16","KCNE1B","CBSL","U2AF1L5","CRYAA2",
  "SIK1B","LOC102724788","OPN1MW3","OPN1MW2")
data=data[-which(is.element(data[,6],wrong)),]

# separating training and passenger variants
train=data[which(data[,1]=="PLP" | data[,1]=="BLB"),]
train[,1]=as.factor(train[,1])
train[,73]=apply(train[,8:12],1,paste,collapse=":")
gno=train[which(train[,57]=="-7"),]
train=train[which(train[,57]!="-7"),]

posi=unique(train[,73])
n=length(posi)/10
ss=c(rep(1,n),rep(2,n),rep(3,n),rep(4,n),rep(5,n),rep(6,n),rep(7,n),rep(8,n),rep(9,n),rep(10,n))
ss=sample(ss)

save(train,posi,ss,gno, file=paste(here,dataLoc,"data5-testing.RData",sep=""))

load(file=paste(here,dataLoc,"data5-testing.RData",sep=""))

save(train,posi,ss, file=paste(here,dataLoc,"data5-testing-only.RData",sep=""))

out=c("BLB",2,3,4,5,6)

for (party in 1:10){
  print(party)
  load(file=paste(here,dataLoc,"data5-testing-only.RData",sep=""))
  a=which(ss==party)

  passenger=train[which(is.element(train[,73],posi[a])),]

  temp=passenger[,c(1,8,9,11,12,13)]
  temp[,6]=party

  out=rbind(out,unique(temp))
}
out=out[-1,]

write.table(out,file=paste(here,dataLoc,"sets-cross-validation.tsv",sep=""),quote=F,sep="\t",row.names=F)

names=c("ClinPred","ADA","RF","CONDEL","SIFT","SIFT4G","PolyPhen-HDIV","PolyPhen-HVAR","LRT","MutationTaster","MutationAssessor","FATHMM",
  "PROVEAN","VEST4","metaSVM","metaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL","fathmm-XF",
  "eigen","eigen-PC","GenoCanyon","fitCons","GERP++NR","GERP++RS","PhyloP100way","PhyloP30way","PhyloP17way","phastCons100way","phastCons30way",
  "phastCons17way","SiPhy29way","MutScore")
nb=c(16:55,74)
restrain=cbind(names,nb)
restest=restrain
resboth=restrain

for (party in 1:10){
  load(file=paste(here,dataLoc,"data5-testing.RData",sep=""))
  a=which(ss==party)

  passenger=train[which(is.element(train[,73],posi[a])),]
  train=train[which(is.element(train[,73],posi[a])==F),]

  train=rbind(train,gno)

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

  # annotating aa change score for all variants
  for (i in 1:dim(freq)[1]){
    print(i)
    train[which(train[,14]==freq[i,1]),69]=freq[i,4]
    passenger[which(passenger[,14]==freq[i,1]),69]=freq[i,4]
  }
  train[,69]=as.numeric(as.character(train[,69]))
  passenger[,69]=as.numeric(as.character(passenger[,69]))

  ####### positional score
  train[,6]=as.character(train[,6])
  train[,1]=as.factor(as.character(train[,1]))
  passenger[,6]=as.character(passenger[,6])
  passenger[,1]=as.factor(as.character(passenger[,1]))

  PLPgene=unique(train[which(train[,1]=="PLP"),6])
  passenger2=passenger

  genelist=unique(c(as.character(train[,2]),as.character(passenger2[,2])))

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

  i="LOW"
  seltrain=which(train[,72]==i)
  train[seltrain,70]=0
  selpass=which(passenger2[,72]==i)
  passenger2[selpass,70]=0

  passenger=passenger2
  rm(passenger2)

  gnomad=train[which(train[,57]==(-7)),]
  train=train[-which(train[,57]==(-7)),]

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

  passenger1=passenger

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
  passenger=passenger1[which(taken==1),]

  passenger[,57]="-8"

  data=rbind(train,passenger)

  # putting median where missing values

  a=duplicated(data[,73])

  ### PREDICTORS
  #SIFT and others
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

  passenger=data[which(data[,57]=="-8"),]
  train=data[which(data[,57]!="-8"),]


  table(train[,1])
  genes=unique(train[which(train[,1]=="PLP"),6])
  train=train[which(is.element(train[,6],genes)),]
  table(train[,1])

  #######
  train[which(is.na(train[,70])),70]=0
  passenger[which(is.na(passenger[,70])),70]=0

  colnames(train)[1]="ClinVar"
  colnames(passenger)[1]="ClinVar"
  train[,1]=as.character(train[,1])
  train[,1]=as.factor(train[,1])

  rf <- randomForest(
    ClinVar ~ V20 + V24 + V21 + V28 + V48 + V49 + V50 + V51 + V52 + V53 + V54 + V55 + V17 + V18 + V58 + V59 + V69 + V70,
    data=train, ntree = 1000, importance = TRUE
  )

  rf.roc<-roc(train$ClinVar,rf$votes[,2])
  auc(rf.roc)

  pred2 = predict(rf, newdata=passenger,"prob")
  rf2.roc<-roc(passenger$ClinVar,pred2[,2])
  auc(rf2.roc)


  train[,74]=rf$votes[,2]
  train2=train[which((train[,1]=="PLP")),]
  train3=train[which(train[,1]=="BLB"),]
  train4=rbind(train2,train3)
  data=train4
  datatrain=data
  names=c("ClinPred","ADA","RF","CONDEL","SIFT","SIFT4G","PolyPhen-HDIV","PolyPhen-HVAR","LRT","MutationTaster",
          "MutationAssessor","FATHMM","PROVEAN","VEST4","metaSVM","metaLR","M-CAP","REVEL","MutPred","MVP","MPC","PrimateAI","DEOGEN2","CADD","DANN","fathmm-MKL",
          "fathmm-XF","eigen","eigen-PC","GenoCanyon","fitCons","GERP++NR","GERP++RS","PhyloP100way","PhyloP30way","PhyloP17way","phastCons100way","phastCons30way","phastCons17way","SiPhy29way")
  nb=c(16:55)
  res=cbind(names,names,names)
  res[,2]=nb
  for (i in 1:dim(res)[1]){
    res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
  }
  restrain=cbind(restrain,c(res[,3],auc(roc(data$ClinVar,data[,74]))))

  passenger[,74]=pred2[,2]
  train2=passenger[which((passenger[,1]=="PLP")),]
  train3=passenger[which(passenger[,1]=="BLB"),]
  train4=rbind(train2,train3)
  data=train4
  datatest=data
  res=cbind(names,names,names)
  res[,2]=nb
  for (i in 1:dim(res)[1]){
    res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
  }
  restest=cbind(restest,c(res[,3],auc(roc(data$ClinVar,data[,74]))))

  data=rbind(datatrain,datatest)
  res=cbind(names,names,names)
  res[,2]=nb
  for (i in 1:dim(res)[1]){
    res[i,3]=round(auc(roc(data$ClinVar,data[,as.numeric(res[i,2])])),digits=3)
  }
  resboth=cbind(resboth,c(res[,3],auc(roc(data$ClinVar,data[,74]))))

  write.matrix(restrain,file=paste(here,dataLoc,"restrain.txt",sep=""))
  write.matrix(restest,file=paste(here,dataLoc,"restest.txt",sep=""))
  write.matrix(resboth,file=paste(here,dataLoc,"resboth.txt",sep=""))

}

save(restrain,restest,resboth, file=paste(here,dataLoc,"res-tests.RData",sep=""))

load(file=paste(here,dataLoc,"res-tests.RData",sep=""))

tr=as.numeric(restrain[41,3:12])
te=as.numeric(restest[41,3:12])

pdf(file=paste(here,plot,"Cross-validation-boxplot.pdf",sep=""),width=6,height=6)
boxplot(tr,te,ylim=c(0.9,1),names=c("Training","Validation"),main="Boxplot of AUCs for 10-fold cross-validation on the training set")
points(rep(1,10),tr)
points(rep(2,10),te)
lines(c(1,1),c(0.96,0.99))
lines(c(2,2),c(0.96,0.99))
lines(c(1,2),c(0.99,0.99))
text(1.5,0.993,paste("n.s. (p=",round(t.test(tr,te)$p.value,digits=3),")",sep=""))
dev.off()


