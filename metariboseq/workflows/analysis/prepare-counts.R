#Takes arguments - 1: fasta .RData file, 2: Ribo bam .RData file, 3: output counts .RData file
library(riboSeqR)
library(tools)

args = commandArgs(trailingOnly=TRUE)

chlamyFastaState <- (args[1])
load(chlamyFastaState) # fastaCDA from prepare-fastaCDA.R

ribofilesState <-  (args[2])
load(ribofilesState) # riboDat from prepare-ribodat.R

#Count reads across ORFs, keeping track of what frame they fall into
fCs <- frameCounting(riboDat, fastaCDS, length = c(20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38))
fS <- readingFrame(rC = fCs, lengths = c(20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38))

save(fCs, fS, file=args[3])
