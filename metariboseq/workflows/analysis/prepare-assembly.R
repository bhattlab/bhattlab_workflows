# Takes arguments - 1: fasta file
library(riboSeqR)

chlamyFasta <- (args[1])
chlamyFastaState <- paste(basename(chlamyFasta), ".RData", sep = '')
fastaCDS <- findCDS(fastaFile = chlamyFasta,startCodon = c("ATG"),stopCodon = c("TAG", "TAA", "TGA"))

# Load in asssembly and annotate ORFs
save(fastaCDA, file=chlamyFastaState)
