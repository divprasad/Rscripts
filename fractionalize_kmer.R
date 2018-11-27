
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("one argument must be supplied ", call.=FALSE)
} 

#for viral contigs

#setwd('/home/divyae/Kmeroutput/Vir/comb/trim/')
for (i in 1:args[1])
  { 
  X<-read.table(file=paste(i,"_merTrim_vir.tsv", sep=""),sep = '\t', header = FALSE)
  #Y<-X
  #Y$contig_id <- NULL
  #smR1<-rowSums(Y)
  #cols=ncol(X)
  smR=rowSums(X)
  Z<-X/smR
  write.table(Z, file=paste(i,"_merTfrac_vir.tsv", sep=""), sep='\t', row.names=FALSE, col.name=FALSE, quote =FALSE)
  
  }


#for bacterial c
#setwd('/home/divyae/Kmeroutput/Bac/Bac1/trim/')
for (i in 1:args[1])
{ 
  X<-read.table(file=paste(i,"_merTrim_bac.tsv", sep=""),sep = '\t', header = FALSE)
  #Y<-X
  #Y$contig_id <- NULL
  #smR1<-rowSums(Y)
  #cols=ncol(X)
  smR=rowSums(X)
  Z<-X/smR
  write.table(Z, file=paste(i,"_merTfrac_bac.tsv", sep=""), sep='\t', row.names=FALSE, col.name=FALSE, quote =FALSE)
  
}





#Y$Afrac<-X$A/smR
#Y$Cfrac<-X$C/smR

#Z1<-X$contig_id
#Z1$Afrac<-X[,1]/smR
#Z1$Cfrac<-X[,2]/smR
#Z1$A <- NULL
#Z1$C <- NULL
