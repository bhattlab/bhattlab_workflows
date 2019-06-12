# post processing of binning results in a more intelligent way
# than that awful shell script we were using
options(stringsAsFactors = F)
# use the snakemake inputs iteratively
# prokka
prokka.files <- snakemake@input[['prokka']]
quast.files <- snakemake@input[['quast']]
checkm.files <- snakemake@input[['checkm']]
trna.files <- snakemake@input[['trna']]
rrna.files <- snakemake@input[['rrna']]
classify.files <- snakemake@input[['classify']]
coverage.files <- snakemake@input[['coverage']]
bins <- snakemake@params[['bins']]
sample.name <- snakemake@params[['sample']]

names(prokka.files) <- bins
names(quast.files) <- bins
names(trna.files) <- bins
names(rrna.files) <- bins
names(coverage.files) <- bins

# process prokka
# but skip the unbinned bin
bin.gene.count <- c()
for (b in bins){
  if (b %in% c('bin.unbinned', 'bin.unbinned.contigs')){
    bin.gene.count[b] <- NA
  } else {
    f <- prokka.files[b]
    fl <- readLines(f)
    gene.count <- length(grep('CDS', fl, value = T))
    bin.gene.count[b] <- gene.count
  }
}

# process quast
quast.df <- data.frame()
for (b in bins){
  f <- quast.files[b]
  df <- read.table(f, sep='\t', comment.char = '', quote = '')
  if(nrow(quast.df)==0){
    quast.df <- df
  } else {
    quast.df <- cbind(quast.df, df[,2])
  }
}
quast.df <- as.data.frame(t(quast.df))
colnames(quast.df) <- quast.df[1,]
quast.df <- quast.df[2:nrow(quast.df), ]
rownames(quast.df) <- quast.df$Assembly

# process checkm
checkm.df <- read.table(checkm.files[1], sep='\t', fill=T, comment.char = '', header=T)
rownames(checkm.df) <- checkm.df$Bin.Id

# process trna
bin.trna.count <- c()
for (b in bins){
  f <- trna.files[b]
  fl <- readLines(f)
  total.string <- grep('Total', fl, value = T)
  total.number <- as.numeric(tail(strsplit(total.string, split=' ')[[1]],1))
  bin.trna.count[b] <- total.number
}

# process rrna
bin.16S.count <- c()
bin.23S.count <- c()
bin.5S.count <- c()
for (b in bins){
  f <- rrna.files[b]
  fl <- readLines(f)
  count.16S <- length(grep('16S', fl, value = T))
  count.23S <- length(grep('23S', fl, value = T))
  count.5S <- length(grep('5S', fl, value = T))
  bin.16S.count[b] <- count.16S
  bin.23S.count[b] <- count.23S
  bin.5S.count[b] <- count.5S
}

# process classify
classify.df <- read.table(classify.files[1], sep='\t', header=T)
classify.df$Bin <- gsub('.fa', '', classify.df$Bin)
rownames(classify.df) <- classify.df$Bin

# process coverage
coverage.list <- list()
for (b in bins){
  f <- coverage.files[b]
  df <- read.table(f, sep='\t')
  coverage.list[[b]] <- as.character(df[1,])
  
}
coverage.df <- do.call(rbind, coverage.list)
colnames(coverage.df) <- c('Sample', 'Bin', 'Coverage')

# merge them all together
out.df <- data.frame(Bin=bins,
                     Sample=sample.name, 
                     Genes=bin.gene.count[bins])
out.df <- cbind(out.df, quast.df[bins, 2:ncol(quast.df)])
out.df <- cbind(out.df, checkm.df[bins, 2:ncol(checkm.df)])
out.df$tRNA <- bin.trna.count[bins]
out.df$rna.16S <- bin.16S.count[bins]
out.df$rna.23S <- bin.23S.count[bins]
out.df$rna.5S <- bin.5S.count[bins]
out.df <- cbind(out.df, classify.df[bins, 2:ncol(classify.df)])
out.df$Coverage <- coverage.df[bins, 'Coverage']
out.df[is.na(out.df)] <- 0
out.df <- out.df[order(out.df$Bin),]

simple.columns <- c("Sample", "Bin", "lca_species", "lca_level", 
                    "lca_fraction", "best_species", "best_level", 
                    "best_fraction", "Size.Mb", "Coverage", "Completeness", 
                    "Contamination", "Strain.heterogeneity", "# contigs (>= 0 bp)", 
                    "Largest contig", "N50", "N75")
out.df.simple <- out.df[,simple.columns]

write.table(out.df, snakemake@output[["full"]], sep = "\t", quote = F, row.names = F)
write.table(out.df.simple, snakemake@output[["simple"]], sep = "\t", quote = F, row.names = F)
