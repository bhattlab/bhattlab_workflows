# script to cluster and make a heatmap of sourmash jaccard 
library(gplots)
library(viridis)
options(stringsAsFactors = F)

# testing args
# compare.file <- '~/scg/HIV_microbiome/allsamples_preprocessed/sourmash/04_sourmash_compare/compare_k21.csv'
# outdir <- '~/scg/HIV_microbiome/allsamples_preprocessed/sourmash/04_sourmash_compare/'

compare.file <- snakemake@input[[1]]
outdir <- snakemake@params[['outdir']]
out.basename <- strsplit(basename(compare.file), split='\\.')[[1]][1]

compare.mat <- read.table(compare.file, sep=',', header=F)
sample.names <- gsub('_concat.fq.abundtrim','', basename(as.character(compare.mat[1,])))
compare.mat <- data.matrix(compare.mat[2:nrow(compare.mat),])
colnames(compare.mat) <- sample.names
rownames(compare.mat) <- sample.names
nsamp <- nrow(compare.mat)
# set diag to NA
diag(compare.mat) <- NA
# cut at 99.9% quantile
compare.thresh <- quantile(compare.mat, c(0.9999), na.rm = T)
compare.mat[compare.mat > compare.thresh] <- compare.thresh

hclust.methods <- c('ward.D2','single','complete')
hclust.outfiles <- sapply(hclust.methods, function(x) 
    file.path(outdir, paste(out.basename, '_heatmap_', x, '.pdf', sep='')))
names(hclust.outfiles) <- hclust.methods

# file size depends on number of samples
pdf.size <- 10 + (nsamp * 0.05)
for (m in hclust.methods){
    pdf(hclust.outfiles[m], height = pdf.size, width=pdf.size)
    heatmap.2(compare.mat, col = viridis(32), trace='none', na.color = 'grey50', margins = c(10,10),
              main=paste('Sourmash ', out.basename, sep=''),
              cexRow = 0.5, cexCol = 0.5,
              distfun = function(x) as.dist(1-x),
              hclustfun = function(x) hclust(x, method=m))
    dev.off()
}
    
