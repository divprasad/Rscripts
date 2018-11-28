
#library(foreign)
#library(flowCore)
#library(Rtsne)
library(treeWAS)
#library(bugwas)

setwd('/linuxhome/tmp/divyae/CWI/VariationGraphs')
#setwd('/ufs/divyae/GWAS/VariationGraphs')
P2<-as.matrix(read.table(file='bugwas_input.unique_rows.binary',header=TRUE,row.names = 1))
Pid2<-as.matrix(read.table(file='ID2.new', header=TRUE))
#Pid<-as.matrix(read.table(file='bugwas_input.id_phenotype', header=TRUE))
#temp<-read.table(file='temp')
#temp1<-t(temp)
#Pid2t<-t(Pid2)

#Pid3<-as.matrix(read.table(file='1new.ID2.new'))


#Pid3t<-t(Pid3)
#write.table(Pid3t, file = "2new.ID2.tsv", append = FALSE, quote = FALSE, sep = "\t")

snps<-t(P2)
df<-read.table(file='ID2.new', header=TRUE)

dfT<-t(df)
phen<-as.vector(unlist(dfT))
names(phen)<-rownames(dfT)
out<-treeWAS(snps,phen)
save(out, file = "treeWAS_try.RData")
