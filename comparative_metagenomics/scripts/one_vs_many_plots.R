# plot results of the one vs many pipeline
# takes in compiled nucmer report file
library(ggplot2)
# library(reshape2)
options(stringsAsFactors = F)

# testing args
# nucmer.file <- '~/scg/bmt_alldata_preprocessing/compare_one_vs_many/p4018_ecoli_complete_genomes/report_stats.txt'
# outdir <- '~/scg/bmt_alldata_preprocessing/compare_one_vs_many/4018_ecoli_complete_genomes'

nucmer.file <- snakemake@input[[1]]
names <- snakemake@params[['names']]
outdir <- snakemake@params[['outdir']]
# out.basename <- strsplit(basename(compare.file), split='\\.')[[1]][1]
# write.table(sample.names, '~/krak_names.txt', sep='\t', quote=F, row.names = F, col.names = F)
# sample.names <- c(read.table('~/scg/krak_names.txt')[,1])

nucmer.df <- read.table(nucmer.file, sep='\t', quote='', header=F, colClasses = c('character', 'numeric','numeric'), comment.char = '')
colnames(nucmer.df) <- c('comp.name', 'avg.identity', 'alignment.length')
rownames(nucmer.df) <- nucmer.df$comp.name

# get names of samples from comp name
nucmer.df$sample.1 <- sapply(nucmer.df$comp.name, function(x) strsplit(x, split='__')[[1]][1])
nucmer.df$sample.2 <- sapply(nucmer.df$comp.name, function(x) strsplit(x, split='__')[[1]][2])

# filter to greater than zero
nucmer.df <- nucmer.df[nucmer.df$alignment.length > 0 & nucmer.df$avg.identity >0, ]

# filter for plotting
nucmer.df.plot <- nucmer.df[nucmer.df$alignment.length >= max(nucmer.df$alignment.length) * 0.95, ]
pdf(file.path(outdir, 'pident_scatter.pdf'))
g <- ggplot(nucmer.df.plot, aes(x=alignment.length, y=avg.identity)) + 
      geom_point(alpha=0.5) + 
      labs(title=paste(names[1], 'reference genome comparison')) 
ggExtra::ggMarginal(g, type = "histogram")
dev.off()

# top comparisons
nucmer.df <- nucmer.df[order(nucmer.df$avg.identity, decreasing = T), ]
# save to file
write.table(nucmer.df, file.path(outdir, 'nucmer_results.txt'), sep='\t', row.names = F, col.names = T, quote=F)