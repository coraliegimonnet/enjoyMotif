#!/bin/Rscript

args <- commandArgs(TRUE)

data <- read.table(args[1], header = T, sep =",")

Random=data$random
Match=data$match

row.names(data) = data[,1]
data1 = data[,-1]

result = kruskal.test(data1)

sink(args[2])
print(args[1])
print(result)
sink()