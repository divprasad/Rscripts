

#install.packages("Rtsne")

# load packages
library(foreign)
library(flowCore)
library(Rtsne)
library(treeWAS)

setwd('/linuxhome/tmp/divyae/CWI/VariationGraphs')

P<-as.matrix(read.table(file='bugwas_input.unique_rows.binary',header=TRUE))
#Pid<-as.matrix(read.table(file='bugwas_input.id_phenotype', header=TRUE))

N<-P[-c(1),]
N1<-N[,-c(1)]


Pid2<-as.matrix(read.table(file='ID2.new', header=TRUE))

Nid<-t(Pid2)



out<-treeWAS(N1,Nid)




#####################






Ntop<-data.matrix(N1)
Nflip<-t(Ntop)

nsub <- 10000
set.seed(123)  # set random seed
data <- data[sample(1:nrow(N2), nsub), ]


Ntop.pca <- prcomp(Ntop, center = TRUE, scale. = TRUE) 
plot(Ntop.pca, type = "l")
PCA.top<-summary(Ntop.pca)

Nflip.pca <- prcomp(Nflip, center = TRUE, scale. = FALSE) 
plot(Nflip.pca, type = "l")
PCA.flip<-summary(Nflip.pca)




N.svd<-svd(Ntop)
#SVD<-N.svd$d
X<-N.svd$d
U<-N.svd$u
V<-N.svd$v
write.table(X, file='X.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(U, file='U.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(V, file='V.tsv', quote=FALSE, sep='\t', col.names = NA)

screeplot(N.pca)


#g = ggbiplot(pca, choices = 1:2, obs.scale = 1, var.scale = 1, groups = dataBC[,13], ellipse = T, circle = T)
#g+theme(legend.position='top') + scale_colour_manual(values=c('green','red'))