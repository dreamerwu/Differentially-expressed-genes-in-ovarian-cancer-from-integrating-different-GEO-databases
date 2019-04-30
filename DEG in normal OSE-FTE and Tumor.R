###########################################################################
#To calculate differentially expressed genes using multi-batches data from GEO datasets
#install packages
source("https://bioconductor.org/biocLite.R")
biocLite("genefilter")
library(genefilter)
install.packages("Rcpp")
library("Rcpp")
biocLite("sva")
library(sva)
install.packages("pamr")
library(pamr)
library(limma)



#Step1: read CEL data, normalize and calculate expression
library("affy")
raw_data=ReadAffy(celfile.path="~/../../")
normalized_data=rma(raw_data)
write.csv(normalized,"~/../../../../output.csv")    #to find the sample name and order
#normalized_data2=read.delim("~/./../../output.csv",head=T,sep="\t")
#sva_data=as.numeric(as.character(normalized_data2))

matrix_normalized_data=as.matrix(normalized_data)
matrix_normalized_data=as.matrix(matrix_normalized_data)
t(matrix_normalized_data)


#Step2:remove batch effect
library("sva")
pheno=read.delim("~/../../")   #pheno is a list using mode() to test
mod=model.matrix(~as.factoor(group),data=pheno)
mod0=model.matrix(~1,data=pheno)
batch=pheno$batch
modcombat=model.matrix(~1,data=pheno)
adjusted_normalized_data=ComBat(dat=matrix_normalized_data,batch=batch,mod=modcombat,par.prior=TRUE,prior.plots=FALSE)
write.csv(adjusted_normalized_data,"~/.././../xxx.csv")
pValuesComBat=f.pvalue(adjusted_normalized_data,mod,mod0)
qValuesComBat=p.adjust(pValuesComBat,method="BH")
F_test_result=cbind(pValuesComBat,qValuesComBat)
write.csv(F_test_result,"~/../../../..")


#Step3:calculate differentially expressed genes
#input=read.csv("~/../../../../xxx.csv")
#matrix=data.matrix(input)
#adjusted_normalized_data=matrix[,2:1041]



list=read.delim("~/../../../",head=T,sep="\t")
mod_limma=model.matrix(~ 0+factor(list))
colnames(mod_limma)=c("xxx","xxxx")
fit=lmFit(adjusted_normalized_data,mod_limma)
contrast.matrix=makeContrasts(xxx-xxxx,levels=mod_limma)
fit2=contrasts.fit(fit,contrast.matrix)
fit22=eBayes(fit2)
result=toptable(fit22,number=60000,adjust.method="BH",p.value=1,lfc=0)
write.csv(result,"~/../../../..")


















