# plot results of the snippy one vs many pipeline
# takes in compiled nucmer report file
library(ggplot2)
library(rafalib)
options(stringsAsFactors = F)

## testing args
# basedir <- '~/scg/bmt_alldata_preprocessing/snippy_compare_one_vs_many/p4018_ecoli_complete_genomes/01_snippy/'
# comparison.names <- list.dirs(basedir, full.names = F,recursive = F)
# outdir <- '~/scg/bmt_alldata_preprocessing/snippy_compare_one_vs_many/p4018_ecoli_complete_genomes'
# sample.name <- 'p4018'


comparison.names <- snakemake@params[['comparison_names']]
outdir <- snakemake@params[['outdir']]
basedir <- snakemake@params[['basedir']]
sample.name <- snakemake@params[['sample_name']]

reference.names <- sapply(comparison.names, function(x) strsplit(x, split="__")[[1]][2])

# for each comparison, read in the snp files
snp.count.list <- list()
for (comp in comparison.names){
  # print(comp)
  ref <- reference.names[comp]
  # print(ref)
  snp.count.f <- file.path(basedir, comp, "snps_filtered_count_by_type.txt")
  depth.count.f <- file.path(basedir, comp, "chromosome_min_depth_positions.txt")
  snp.count <- read.table(snp.count.f, sep='\t', col.names = c('chr', 'type', 'count'))
  depth.count <- read.table(depth.count.f, sep='\t', col.names = c('chr', 'total.bases', 'min.depth.bases'))
  # check to make sure we have lines in these inputs
  if ((nrow(snp.count) > 0) & (nrow(depth.count) > 0)){
      snp.count$count <- as.numeric(snp.count$count)
      depth.count$total.bases <- as.numeric(depth.count$total.bases)
      depth.count$min.depth.bases <- as.numeric(depth.count$min.depth.bases)
      # probably filter this to the chromosome, which is the chr wit the most bases
      keep.chr <- depth.count$chr[order(depth.count$total.bases, decreasing = T)[1]]
      depth.count <- depth.count[depth.count$chr==keep.chr, ]
      snp.count <- snp.count[snp.count$chr==keep.chr, ]
      # remove any NA rows
      snp.count <- snp.count[!(apply(snp.count, 1, function(x) any(is.na(x)))),]
      # fraction depth
      depth.count$frac.min.depth <- depth.count$min.depth.bases / depth.count$total.bases 
      # add total count
      snp.count <- rbind(snp.count, c(snp.count$chromosome[1], 'TOTAL', sum(as.numeric(snp.count$count), na.rm=T), ref))
      snp.count$count <- as.numeric(snp.count$count)
      # limit to certain types of variants and stuff the rest into other
      keep.variants <- c('TOTAL', 'SNP', 'MNP', 'INDEL', 'OTHER')
      add.to.other <- sum(snp.count[!(snp.count$type %in% keep.variants), 'count'])
      snp.count[snp.count$type=='OTHER', 'count'] <- snp.count[snp.count$type=='OTHER', 'count'] + add.to.other
      snp.count <- snp.count[snp.count$type %in% keep.variants,]
      # per covered kb
      snp.count$count.per.covered.kb <- sapply(snp.count$count, function(x) as.numeric(x)/depth.count$min.depth.bases[1]*1000)
      # just add depth info here
      snp.count$total.chr.bases <- depth.count$total.bases
      snp.count$total.covered.bases <- depth.count$min.depth.bases
      snp.count$frac.covered <- depth.count$frac.min.depth
      # add reference name column
      snp.count$reference <- ref
      snp.count$sample.name <- sample.name
      snp.count.list[[ref]] <- snp.count
    }
}

snp.count.df <- do.call(rbind, snp.count.list)
snp.count.df[,'sample.name'] <- sample.name
snp.count.df$count <- as.numeric(snp.count.df$count)
# reorganize some columns
col.order <- c("sample.name","reference","chr","type","count","count.per.covered.kb",
    "total.chr.bases","total.covered.bases","frac.covered")
snp.count.df <- snp.count.df[, col.order]
# print(head(snp.count.df))
# remove NA counts - where do these come from?
write.table(snp.count.df, file.path(outdir, 'snp_count_df.txt'), sep='\t', quote=F)
snp.count.df <- snp.count.df[!(is.na(snp.count.df$count)),]
snp.count.df <- snp.count.df[!(is.na(snp.count.df$count.per.covered.kb)),]
write.table(snp.count.df, file.path(outdir, 'snp_count_df_nna.txt'), sep='\t', quote=F)
# total only
snp.count.df$count <- as.numeric(snp.count.df$count)
snp.count.df$count.per.covered.kb <- as.numeric(snp.count.df$count.per.covered.kb)
snp.count.total.df <- snp.count.df[snp.count.df$type=='TOTAL', ]
write.table(snp.count.total.df, file.path(outdir, 'snp_count_df_total.txt'), sep='\t', quote=F)


# tables to save
top.total.df <- snp.count.total.df[order(snp.count.total.df$count.per.covered.kb), ][1:10,]
top.total.references <- top.total.df$reference
top.all.df <- do.call(rbind, lapply(top.total.references, function(x) snp.count.df[snp.count.df$reference == x, ]))
total.df <- snp.count.total.df[order(snp.count.total.df$count.per.covered.kb), ]
total.references <- total.df$reference
all.df <- do.call(rbind, lapply(total.references, function(x) snp.count.df[snp.count.df$reference == x, ]))
# output files
outf.1 <- file.path(outdir, 'top10_reference_SNP_total.txt')
outf.2 <- file.path(outdir, 'top10_reference_SNP_alltypes.txt')
outf.3 <- file.path(outdir, 'all_reference_SNP_total.txt')
outf.4 <- file.path(outdir, 'all_reference_SNP_alltypes.txt')
# save tables
write.table(top.total.df, outf.1, sep='\t', quote=F, row.names = F, col.names = T)
write.table(top.all.df, outf.2, sep='\t', quote=F, row.names = F, col.names = T)
write.table(total.df, outf.3, sep='\t', quote=F, row.names = F, col.names = T)
write.table(all.df, outf.4, sep='\t', quote=F, row.names = F, col.names = T)

out.plot <- file.path(outdir, 'variant_distribution_plots.pdf')
pdf(out.plot, width = 9, height=6)
ggplot(snp.count.df, aes(x=frac.covered, y=count)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_wrap(~type, scales='free') + 
  labs(title='Variants compared to reference genomes, by type', x = 'Fraction of reference covered >= 10x', y='Count of variants')

ggplot(snp.count.df, aes(x=frac.covered, y=count.per.covered.kb)) +
  geom_point(alpha=0.5) +
  theme_bw() +
  facet_wrap(~type, scales='free') + 
  labs(title='Variants compared to reference genomes, by type', x = 'Fraction of reference covered >= 10x', y='Variants per covered kb')

mypar(2, 1)
hist(snp.count.total.df$count, breaks=50, main='Histogram of total variant counts', xlab='Total variants', col='steelblue')
hist(snp.count.total.df$count.per.covered.kb, breaks=50, xlab='Total variants per covered kb', main='', col='steelblue')
dev.off()



