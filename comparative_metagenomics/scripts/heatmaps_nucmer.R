# script to cluster and make a heatmap nucmer
library(gplots)
library(viridis)
# library(reshape2)
options(stringsAsFactors = F)

# testing args
# nucmer.file <- '~/scg_lab/transmit_crass/compare_krak_reads/crassphage/04_nucmer_pairwise/report_stats.txt'
# outdir <- '~/scg_lab/transmit_crass/compare_krak_reads/crassphage/04_nucmer_pairwise/plots'

nucmer.file <- snakemake@input[[1]]
sample.names <- snakemake@params[['sample_names']]
outdir <- snakemake@params[['outdir']]
min.identity.plot <- as.numeric(snakemake@params[['min_identity_plot']])
# out.basename <- strsplit(basename(compare.file), split='\\.')[[1]][1]
# write.table(sample.names, '~/krak_names.txt', sep='\t', quote=F, row.names = F, col.names = F)
# sample.names <- c(read.table('~/scg/krak_names.txt')[,1])

nucmer.df <- read.table(nucmer.file, sep='\t', quote='', header=F, colClasses = c('character', 'numeric','numeric'))
colnames(nucmer.df) <- c('comp.name', 'avg.identity', 'alignemnt.length')
rownames(nucmer.df) <- nucmer.df$comp.name

# get names of samples from comp name
nucmer.df$sample.1 <- sapply(nucmer.df$comp.name, function(x) strsplit(x, split='__')[[1]][1])
nucmer.df$sample.2 <- sapply(nucmer.df$comp.name, function(x) strsplit(x, split='__')[[1]][2])

# First filter for everything - at least 1kb
min.align.length <- 1000
nucmer.df <- nucmer.df[nucmer.df$alignemnt.length >= min.align.length, ]

# filter to alignments at least this long
# without a reference genome, a way to do this is the alignment 
# must be 75% of the 90% quantile of all alignments
min.align.length <- quantile(nucmer.df$alignemnt.length, c(0.90)) * 0.75
nucmer.df.filt <- nucmer.df[nucmer.df$alignemnt.length >= min.align.length, ]

all.sample.names <- unique(c(nucmer.df$sample.1, nucmer.df$sample.2))
nsamp <- length(all.sample.names)
all.sample.names.filt <- unique(c(nucmer.df.filt$sample.1, nucmer.df.filt$sample.2))
nsamp.filt <- length(all.sample.names.filt)

# fill this in with a for loop oh god WHY
nucmer.mat <- matrix(0, nrow=nsamp, ncol=nsamp, dimnames=list(all.sample.names, all.sample.names))
length.mat <- matrix(0, nrow=nsamp, ncol=nsamp, dimnames=list(all.sample.names, all.sample.names))
for (i in 1:nrow(nucmer.df)){
    x <- nucmer.df[i, ]
    nucmer.mat[x[1,'sample.1'], x[1,'sample.2']] <- x[1,'avg.identity']
    nucmer.mat[x[1,'sample.2'], x[1,'sample.1']] <- x[1,'avg.identity']
    length.mat[x[1,'sample.1'], x[1,'sample.2']] <- x[1,'alignemnt.length']
    length.mat[x[1,'sample.2'], x[1,'sample.1']] <- x[1,'alignemnt.length']
}
# filtered version
nucmer.mat.filt <- matrix(0, nrow=nsamp.filt, ncol=nsamp.filt, dimnames=list(all.sample.names.filt, all.sample.names.filt))
length.mat.filt <- matrix(0, nrow=nsamp.filt, ncol=nsamp.filt, dimnames=list(all.sample.names.filt, all.sample.names.filt))
for (i in 1:nrow(nucmer.df.filt)){
    x <- nucmer.df.filt[i, ]
    nucmer.mat.filt[x[1,'sample.1'], x[1,'sample.2']] <- x[1,'avg.identity']
    nucmer.mat.filt[x[1,'sample.2'], x[1,'sample.1']] <- x[1,'avg.identity']
    length.mat.filt[x[1,'sample.1'], x[1,'sample.2']] <- x[1,'alignemnt.length']
    length.mat.filt[x[1,'sample.2'], x[1,'sample.1']] <- x[1,'alignemnt.length']
}

# set diagonal to NA for plotting
# except in the case of very small matrices
if(nrow(nucmer.mat) >2){ diag(nucmer.mat) <- NA
} else { diag(nucmer.mat) <- 100 }
if(nrow(nucmer.mat.filt) >2){ diag(nucmer.mat.filt) <- NA 
} else {diag(nucmer.mat.filt) <- 100 }
    
# set cells that are 0 to NA for easier visualization
nucmer.mat.plot <- nucmer.mat
nucmer.mat.plot[nucmer.mat.plot==0] <- NA
nucmer.mat.filt.plot <- nucmer.mat.filt
nucmer.mat.filt.plot[nucmer.mat.filt.plot==0] <- NA

# if limiting to a certain percentage, set them to NA 
if (min.identity.plot >0) {
    nucmer.mat.plot[nucmer.mat.plot < min.identity.plot] <- NA    
    nucmer.mat.filt.plot[nucmer.mat.filt.plot < min.identity.plot] <- NA    
}

# unfiltered version
if (all(nucmer.mat[!is.na(nucmer.mat)] == 0)){
    pdf(file.path(outdir, 'nucmer_identity_heatmap_complete.pdf'))
    plot(1,1, col='white')
    text(1,1, "Not enough nonzero samples for heatmap")
} else {
    hclust.methods <- c('ward.D2','single','complete')
    hclust.outfiles <- sapply(hclust.methods, function(x) 
        file.path(outdir, paste('nucmer_identity_heatmap_unfiltered_', x, '.pdf', sep='')))
    names(hclust.outfiles) <- hclust.methods

    # file size depends on number of samples
    pdf.size <- 10 + (nsamp * 0.05)
    for (m in hclust.methods){
        pdf(hclust.outfiles[m], height = pdf.size, width=pdf.size)
        heatmap.2(nucmer.mat.plot, col = viridis(128), trace='none', 
            na.color = 'grey50', margins = c(10,10),
            main='nucmer assembled genome average identity',
            cexRow = 0.5, cexCol = 0.5,
            distfun = function(x) as.dist(100-nucmer.mat),
            key.xlab = '% identity', key.title = '',
            hclustfun = function(x) hclust(x, method=m))
        dev.off()
    }
}

# filtered version
if (all(nucmer.mat.filt[!is.na(nucmer.mat.filt)] == 0)){
    pdf(file.path(outdir, 'nucmer_identity_heatmap_complete.pdf'))
    plot(1,1, col='white')
    text(1,1, "Not enough nonzero samples for heatmap")
} else {
    hclust.methods <- c('ward.D2','single','complete')
    hclust.outfiles <- sapply(hclust.methods, function(x) 
        file.path(outdir, paste('nucmer_identity_heatmap_', x, '.pdf', sep='')))
    names(hclust.outfiles) <- hclust.methods

    # file size depends on number of samples
    pdf.size <- 10 + (nsamp * 0.05)
    for (m in hclust.methods){
        pdf(hclust.outfiles[m], height = pdf.size, width=pdf.size)
        heatmap.2(nucmer.mat.filt.plot, col = viridis(128), trace='none', 
            na.color = 'grey50', margins = c(10,10),
            main='nucmer assembled genome average identity',
            cexRow = 0.5, cexCol = 0.5,
            distfun = function(x) as.dist(100-nucmer.mat.filt),
            key.xlab = '% identity', key.title = '',
            hclustfun = function(x) hclust(x, method=m))
        dev.off()
    }
}

# note version, using unfiltered 
if (any(nucmer.mat[!is.na(nucmer.mat)] > 0)){
    hclust.methods <- c('ward.D2','single','complete')
    hclust.outfiles.note1 <- sapply(hclust.methods, function(x) 
        file.path(outdir, paste('nucmer_identity_heatmap_unfiltered_note_pident_', x, '.pdf', sep='')))
    hclust.outfiles.note2 <- sapply(hclust.methods, function(x) 
        file.path(outdir, paste('nucmer_identity_heatmap_unfiltered_note_length_', x, '.pdf', sep='')))
    names(hclust.outfiles.note1) <- hclust.methods
    names(hclust.outfiles.note2) <- hclust.methods
    # file size depends on number of samples
    pdf.size <- 10 + (nsamp * 0.05)
    for (m in hclust.methods){
        # note version with % identity
        pdf(hclust.outfiles.note1[m], height = pdf.size, width=pdf.size)
        heatmap.2(nucmer.mat.plot, col = viridis(128), trace='none', 
                  na.color = 'grey50', margins = c(10,10),
                  main='nucmer assembled genome average identity',
                  cexRow = 0.5, cexCol = 0.5,
                  distfun = function(x) as.dist(100-nucmer.mat),
                  hclustfun = function(x) hclust(x, method=m),
                  key.xlab = '% identity', key.title = '',
                  cellnote = round(nucmer.mat, 2), notecex = 0.4, notecol = 'magenta')
        dev.off()
        # note version with length of alignment
        pdf(hclust.outfiles.note2[m], height = pdf.size, width=pdf.size)
        heatmap.2(nucmer.mat.plot, col = viridis(128), trace='none', 
                  na.color = 'grey50', margins = c(10,10),
                  main='nucmer assembled genome average identity',
                  cexRow = 0.5, cexCol = 0.5,
                  distfun = function(x) as.dist(100-nucmer.mat),
                  hclustfun = function(x) hclust(x, method=m),
                  key.xlab = '% identity', key.title = '',
                  cellnote = round(length.mat/1000), notecex = 0.4, notecol = 'magenta')
        dev.off()
    }
}

# note version, using filtered 
if (any(nucmer.mat[!is.na(nucmer.mat)] > 0)){
  hclust.methods <- c('ward.D2','single','complete')
  hclust.outfiles.note1 <- sapply(hclust.methods, function(x) 
    file.path(outdir, paste('nucmer_identity_heatmap_note_pident_', x, '.pdf', sep='')))
  hclust.outfiles.note2 <- sapply(hclust.methods, function(x) 
    file.path(outdir, paste('nucmer_identity_heatmap_note_length_', x, '.pdf', sep='')))
  names(hclust.outfiles.note1) <- hclust.methods
  names(hclust.outfiles.note2) <- hclust.methods
  # file size depends on number of samples
  pdf.size <- 10 + (nsamp * 0.05)
  for (m in hclust.methods){
    # note version with % identity
    pdf(hclust.outfiles.note1[m], height = pdf.size, width=pdf.size)
    heatmap.2(nucmer.mat.filt.plot, col = viridis(128), trace='none', 
              na.color = 'grey50', margins = c(10,10),
              main='nucmer assembled genome average identity',
              cexRow = 0.5, cexCol = 0.5,
              distfun = function(x) as.dist(100-nucmer.mat.filt),
              hclustfun = function(x) hclust(x, method=m),
              key.xlab = '% identity', key.title = '',
              cellnote = round(nucmer.mat.filt, 2), notecex = 0.4, notecol = 'magenta')
    dev.off()
    # note version with length of alignment
    pdf(hclust.outfiles.note2[m], height = pdf.size, width=pdf.size)
    heatmap.2(nucmer.mat.filt.plot, col = viridis(128), trace='none', 
              na.color = 'grey50', margins = c(10,10),
              main='nucmer assembled genome average identity',
              cexRow = 0.5, cexCol = 0.5,
              distfun = function(x) as.dist(100-nucmer.mat.filt),
              hclustfun = function(x) hclust(x, method=m),
              key.xlab = '% identity', key.title = '',
              cellnote = round(length.mat.filt/1000), notecex = 0.4, notecol = 'magenta')
    dev.off()
  }
}

# remove pesky rplots.pdf
if(file.exists('Rplots.pdf')){
    file.remove('Rplots.pdf')
}