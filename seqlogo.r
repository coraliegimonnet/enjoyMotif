#!/bin/Rscript

library(seqLogo)

args <- commandArgs(TRUE)

ma <- read.table(args[1],header = F, sep = "\t",dec = ".",na.strings = "NA",fill = TRUE)

x <- ncol(ma)
matrix = ma[,-x] #suppr la dernière colonne qui est vide (défaut script perl)
p <- makePWM(matrix,alphabet = "DNA")

png(args[2])
seqLogo(p,ic.scale=FALSE)
dev.off()