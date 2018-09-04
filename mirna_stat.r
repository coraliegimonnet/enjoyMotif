#!/bin/Rscript

args <- commandArgs(TRUE)

mirna <- read.table(args[1], header = T, sep =",")

Random=mirna$random
Match=mirna$match

row.names(mirna) = mirna[,1]
mirnas = mirna[,-1]

result = chisq.test(mirnas)

sink(args[2])
print(result)
sink()