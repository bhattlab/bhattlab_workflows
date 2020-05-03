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
  df <- read.table(f, sep='\t', comment.char = '', quote = '', skip=1, row.names=1)
  if(nrow(quast.df)==0){
    quast.df <- df
  } else {
    quast.df <- cbind(quast.df, df[,1])
  }
}
quast.df <- as.data.frame(t(quast.df))
rownames(quast.df) <- bins

# process checkm
checkm.df <- read.table(checkm.files[1], sep='\t', fill=T, comment.char = '', header=T, quote = '')
rownames(checkm.df) <- checkm.df$Bin.Id

# process trna
bin.trna.count <- c()
for (b in bins){
  f <- trna.files[b]
  fl <- readLines(f)
  total.string <- grep('Total', fl, value = T)
  if (length(total.string) >0){
    total.number <- as.numeric(tail(strsplit(total.string, split=' ')[[1]],1))
  } else {
    total.number <- 0
  }
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
classify.df <- read.table(classify.files[1], sep='\t', header=T, quote = '')
classify.df$Bin <- gsub('.fa', '', classify.df$Bin)
rownames(classify.df) <- classify.df$Bin

# process coverage
coverage.list <- list()
for (b in bins){
  f <- coverage.files[b]
  df <- read.table(f, sep='\t', quote = '')
  coverage.list[[b]] <- as.character(df[1,])
  
}
coverage.df <- do.call(rbind, coverage.list)
colnames(coverage.df) <- c('Sample', 'Bin', 'Coverage')

# merge them all together
out.df <- data.frame(Bin=bins,
                     Sample=sample.name, 
                     Genes=bin.gene.count[bins])
out.df <- cbind(out.df, quast.df[bins, ])
out.df <- cbind(out.df, checkm.df[bins, 2:ncol(checkm.df)])
out.df$tRNA <- bin.trna.count[bins]
out.df$rna.16S <- bin.16S.count[bins]
out.df$rna.23S <- bin.23S.count[bins]
out.df$rna.5S <- bin.5S.count[bins]
out.df <- cbind(out.df, classify.df[bins, 2:ncol(classify.df)])
out.df$Coverage <- coverage.df[bins, 'Coverage']
out.df[is.na(out.df)] <- 0
out.df <- out.df[order(out.df$Bin),]

# set fields as numeric for calculating quality
out.df$Completeness <- as.numeric(out.df$Completeness)
out.df$Contamination <- as.numeric(out.df$Contamination)
out.df$N50 <- as.numeric(out.df$N50)
out.df[,'# contigs (>= 0 bp)'] <- as.numeric(out.df[,'# contigs (>= 0 bp)'])
out.df$Coverage <- round(as.numeric(out.df$Coverage), 2)
out.df$tRNA <- round(as.numeric(out.df$tRNA), 2)
out.df$rna.16S <- round(as.numeric(out.df$rna.16S), 2)
out.df$rna.23S <- round(as.numeric(out.df$rna.23S), 2)
out.df$rna.5S <- round(as.numeric(out.df$rna.5S), 2)
# implement some quality calls and standards in this 
out.df$bin.quality.numeric <- round(out.df$Completeness - (5* out.df$Contamination),2) 
# different quality thresholds
out.df$low.quality <- out.df$Completeness < 50 & out.df$Contamination <10
out.df$med.quality <- out.df$Completeness >= 50 & out.df$Contamination <10
out.df$high.quality.nayfach <- out.df$Completeness >= 90 & 
                               out.df$Contamination <= 5 & 
                               out.df$N50 >=10000 & 
                               out.df[,'# contigs (>= 0 bp)'] <= 500 &
                               out.df$Coverage >=5
out.df$high.quality.bowers <- out.df$Completeness >= 90 & 
                              out.df$Contamination <= 5 & 
                              out.df$tRNA >= 18 & 
                              out.df$rna.16S >= 1 & 
                              out.df$rna.23S >= 1 & 
                              out.df$rna.5S >= 1
# ensure low quality is on evrything with higher quality
out.df[out.df$med.quality, 'low.quality'] <- TRUE
quality.map <- c('0) really bad', '1) low quality', '2) medium quality', '3) high quality Nayfach', '4) high quality Bowers')
out.df$bin.quality.call <- quality.map[rowSums(out.df[, c('low.quality', 'med.quality', 'high.quality.nayfach', 'high.quality.bowers')])+1]

new.col.order <- c("Sample", "Bin", "bin.quality.numeric", "bin.quality.call", "Completeness", 
                    "Contamination", "Strain.heterogeneity", "lca_species", "lca_level", 
                    "lca_fraction", "best_species", "best_level", 
                    "best_fraction", "N50", "Size.Mb", "Coverage", "# contigs (>= 0 bp)", 
                    "Largest contig", "Genes", "tRNA", "rna.16S", "rna.23S", "rna.5S",
                    "low.quality", "med.quality", "high.quality.nayfach", "high.quality.bowers",
                    "# contigs (>= 10000 bp)","# contigs (>= 50000 bp)","# contigs (>= 100000 bp)","# contigs (>= 250000 bp)","# contigs (>= 500000 bp)","# contigs (>= 1000000 bp)","# contigs (>= 2000000 bp)","# contigs (>= 3000000 bp)",
                    "Total length (>= 0 bp)","Total length (>= 10000 bp)","Total length (>= 50000 bp)","Total length (>= 100000 bp)","Total length (>= 250000 bp)","Total length (>= 500000 bp)","Total length (>= 1000000 bp)","Total length (>= 2000000 bp)","Total length (>= 3000000 bp)",
                    "# contigs","Total length",
                    "N75","L50","L75","# N's per 100 kbp","Marker.lineage","X..genomes","X..markers","X..marker.sets","X0","X1","X2","X3","X4","X5."
                    )
out.df <- out.df[, new.col.order]
simple.columns <- new.col.order[1:23]
out.df.simple <- out.df[,simple.columns]

write.table(out.df, snakemake@output[["full"]], sep = "\t", quote = F, row.names = F, col.names=T)
write.table(out.df.simple, snakemake@output[["simple"]], sep = "\t", quote = F, row.names = F, col.names=T)
