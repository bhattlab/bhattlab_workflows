# Takes arguments - 1: Ribo bam file, 2: output RData file
library(riboSeqR)

args = commandArgs(trailingOnly=TRUE)

ribofiles <-  (args[1])
#Load in asssembly and annotate ORFs
riboDat <- readRibodata(ribofiles, replicates = 1)
save(riboDat, file=args[2])
