# Takes arguments - 1: fasta file, 2: output RData file
library(riboSeqR)

args = commandArgs(trailingOnly=TRUE)

chlamyFasta <- (args[1])
#Load in asssembly and annotate ORFs
fastaCDS <- findCDS(fastaFile = chlamyFasta,startCodon = c("ATG"),stopCodon = c("TAG", "TAA", "TGA"))
save(fastaCDS, file=args[2])
