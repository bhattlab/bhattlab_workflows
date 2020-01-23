library(ggplot2)
library(reshape2)
library(RColorBrewer)

options(stringsAsFactors = F)
# inf <- '~/scg_lab/bootcamp/zeiser/preprocessing/01_processing/readcounts.tsv'
# outf <- '~/scg_lab/bootcamp/zeiser/preprocessing/01_processing/readcounts.pdf'
inf <- snakemake@input[[1]]
outf <- snakemake@output[[1]]

readcounts <- read.table(inf, sep='\t', quote='', header=T)
readcounts$raw_frac <- 1.0
readcounts.melt.count <- melt(readcounts[,c('Sample', 'raw_reads', 
                                            'dedup_reads', 'trimmed_reads', 
                                            'host_removed_reads')], id.vars = 'Sample') 
readcounts.melt.frac <- melt(readcounts[,c('Sample', 'raw_frac',
                                            'dedup_frac', 'trimmed_frac',
                                            'host_removed_frac')], id.vars = 'Sample') 


g.count <- ggplot(readcounts.melt.count, aes(x=Sample, y=value, fill=variable)) + 
    geom_bar(stat='identity' , position='dodge') + 
    theme_bw() + 
    scale_fill_brewer(palette = 'Set2') +
    labs(title='Readcounts at each processing step', 
         y = 'Read count') + 
    guides(fill=guide_legend(title="Processing level")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
# g.count

g.frac <- ggplot(readcounts.melt.frac, aes(x=Sample, y=value, fill=variable)) + 
    geom_bar(stat='identity' , position='dodge') + 
    theme_bw() + 
    scale_fill_brewer(palette = 'Set2') +
    labs(title='Readcounts at each processing step', 
         y = 'Fraction of Raw ') + 
    guides(fill=guide_legend(title="Processing level")) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
# g.frac

# set plot width to be a function of the number of samples
nsamp <- nrow(readcounts)
plot.width <- 2 + (nsamp/2)
plot.height <- 6

pdf(outf, height = plot.height, width=plot.width)
print(g.count)
print(g.frac)
dev.off()