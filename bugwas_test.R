
library(bugwas)
library(ape)
library(phangorn)

gen <- read.table("/hosts/linuxhome/wildtype1/divyae2/BUGWAS/example/gen.txt", header=TRUE,	sep="\t")
pheno <- read.table("/hosts/linuxhome/wildtype1/divyae2/BUGWAS/example/pheno.txt", header=TRUE,	sep="\t")
phylo <- read.table("/hosts/linuxhome/wildtype1/divyae2/BUGWAS/example/tree.txt", header=TRUE,	sep="\t")
prefix <- "test_bugwas"
gem.path <- "./gemma"

#source('~/divyae2/BUGWAS/bugwas/R/linLocGEMMA.R')
data <- lin_loc(gen = gen, pheno = pheno, phylo = phylo, prefix = prefix, gem.path = gem.path)
data <- linLocGEMMA(gen = gen, pheno = pheno, phylo = phylo, prefix = prefix, gem.path = gem.path)
