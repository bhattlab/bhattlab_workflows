# Heatmaps for instrain results
# in the model of my nucmer heatmaps

# takes in a instrain genome_wide file from many samples
library(gplots)
library(viridis)
options(stringsAsFactors = F)

# testing args
if (F){
# clus <- '22_2'
# instrain.base <- '~/scg/tx88/drep_alignment_comparison/instrain_compare_all/'
# gw.f <- paste(clus, 'genomeWide_compare.tsv', sep='_')
# instrain.f <-file.path(instrain.base, clus, "output", gw.f)
ins.f <- "~/scg/tx88/drep_alignment_comparison/instrain_compare_all/22_2/output/22_2_genomeWide_compare.tsv"
outdir <- '~/scg/tx88/drep_alignment_comparison/instrain_compare_all/22_2/heatmaps'
cluster.name <- 'Lactobacillus rhamnosus'
min.identity.plot <- 0
# relatively lax cutoff of  20% of the genome compared
min.frac.compared <- 0.2
}

# snakemake args
ins.f <- snakemake@input[[1]]
outdir <- snakemake@params[['outdir']]
min.identity.plot <- as.numeric(snakemake@params[['min_identity_plot']])
min.frac.compared <- as.numeric(snakemake@params[['min_frac_compared']])
cluster.name <- snakemake@params[['cluster_name']]

# function to read an instrain report into a manageable
# dataframe
parse_instrain <- function(ins.f, min.frac.compared) {
    # read instrain df
    ins.orig <- read.table(ins.f, sep='\t', header=T, quote='')
    # first remove rows with no coverage
    ins <- ins.orig[!is.na(ins.orig$coverage_overlap),]
    # remove low coverage
    ins <- ins[ins$percent_compared > min.frac.compared, ]
    # order by ani
    ins <- ins[order(ins$popANI, decreasing = T),]
    # fix names by removing .sorted and .bam
    ins$name1 <- gsub('\\.bam', '', ins$name1)
    ins$name2 <- gsub('\\.bam', '', ins$name2)
    ins$name1 <- gsub('\\.sorted', '', ins$name1)
    ins$name2 <- gsub('\\.sorted', '', ins$name2)
    return(ins)
}


# make a heatmap from the instrain data, since thier pipeline doesn't really work
# convert instrain mappig to a matrix of pairwise
ins_to_matrix_list <- function(ins.raw, ani.thresh=0, use.col='popANI'){
    ins <- ins.raw[ins.raw[, use.col] > ani.thresh, ]
    all.sample.names <- unique(c(ins$name1, ins$name2))
    nsamp <- length(all.sample.names)
    # fill this in with a for loop oh god WHY
    ani.mat <- matrix(0, nrow=nsamp, ncol=nsamp, dimnames=list(all.sample.names, all.sample.names))
    length.mat <- matrix(0, nrow=nsamp, ncol=nsamp, dimnames=list(all.sample.names, all.sample.names))
    for (i in 1:nrow(ins)){
        x <- ins[i, ]
        ani.mat[x[1,'name1'], x[1,'name2']] <- x[1, use.col]
        ani.mat[x[1,'name2'], x[1,'name1']] <- x[1, use.col]
        length.mat[x[1,'name1'], x[1,'name2']] <- x[1, 'percent_compared']
        length.mat[x[1,'name2'], x[1,'name1']] <- x[1, 'percent_compared']
    }
    # set diagonal to NA for plotting
    # except in the case of very small matrices
    if(nrow(ani.mat) >2){ diag(ani.mat) <- NA
    } else { diag(ani.mat) <- 1 }
    if(nrow(length.mat) >2){ diag(length.mat) <- NA
    } else { diag(length.mat) <- 1 }

    ani.mat <- ani.mat * 100
    return(list(ani.mat, length.mat))
}

# check for a run without enough samples
ins.test <- read.table(ins.f, nrows=1, sep='\t')
if (ins.test[1,1] == "NOT ENOUGH FILTERD SAMPLES PRESENT"){
    print(snakemake@output[[1]][1])
    pdf(snakemake@output[[1]][1])
    plot(1,1)
    dev.off()
    quit()
}

# load the instrain dataframe
ins <- parse_instrain(ins.f, min.frac.compared)
if (nrow(ins)<2){
    comp <- "popANI"
    outdir.comp <- file.path(outdir, comp)
    if (!(dir.exists(outdir.comp))){dir.create(outdir.comp, recursive = T)}
    pdf(file.path(outdir.comp, paste0(comp, '_heatmap_unfiltered_complete.pdf')))
    plot(1,1, col='white')
    text(1,1, "Not enough nonzero samples for heatmap")
    quit(status=0)
}

# do all of this for popANI and conANI
for (comp in c('popANI', 'conANI')){
    outdir.comp <- file.path(outdir, comp)
    if (!(dir.exists(outdir.comp))){dir.create(outdir.comp, recursive = T)}

    ins.ml <- ins_to_matrix_list(ins, use.col=comp)
    ani.mat <- ins.ml[[1]]
    length.mat <- ins.ml[[2]]
    all.sample.names <- unique(c(ins$name1, ins$name2))
    nsamp <- length(all.sample.names)

    # set cells that are 0 to NA for easier visualization
    ani.mat.plot <- ani.mat
    ani.mat.plot[ani.mat.plot==0] <- NA
    if (min.identity.plot >0) {
        ani.mat.plot[ani.mat.plot < min.identity.plot] <- NA
    }
    # if all values are identical, can't plot the heatmap. 
    # sub one of them
    if(all(ani.mat.plot[!(is.na(ani.mat.plot))] == 100)) {
      ani.mat.plot[1,2] <- 99.9999999
    }

    # unfiltered version
    if (all(ani.mat[!is.na(ani.mat)] == 0)){
        pdf(file.path(outdir.comp, paste0(comp, '_heatmap_complete.pdf')))
        plot(1,1, col='white')
        text(1,1, "Not enough nonzero samples for heatmap")
    } else {
        hclust.methods <- c('ward.D2','single','complete')
        hclust.outfiles <- sapply(hclust.methods, function(x)
            file.path(outdir.comp, paste(comp, '_heatmap_unfiltered_', x, '.pdf', sep='')))
        names(hclust.outfiles) <- hclust.methods
        # layout params
        lmat=rbind( c(4, 3), c(2,1 ) )
        lwid=c(0.2, 0.85)
        lhei=c(0.15, 0.85)

        # file size depends on number of samples
        pdf.size <- 7 + (nsamp * 0.05)
        for (m in hclust.methods){
            pdf(hclust.outfiles[m], height = pdf.size, width=pdf.size)
            heatmap.2(ani.mat.plot, col = viridis(128), trace='none',
                      na.color = 'grey50', margins = c(10,10),
                      main=paste(cluster.name, comp),
                      cexRow = 0.5, cexCol = 0.5,
                      distfun = function(x) as.dist(100-ani.mat),
                      key.xlab = comp, key.title = '',
                      hclustfun = function(x) hclust(x, method=m),
                      lmat=lmat, lhei=lhei, lwid=lwid, key.par = list(mar=c(4,2.5,1,0.5)))
            dev.off()
        }

        # note version, using unfiltered
        if (any(ani.mat[!is.na(ani.mat)] > 0)){
            hclust.methods <- c('ward.D2','single','complete')
            hclust.outfiles.note1 <- sapply(hclust.methods, function(x)
                file.path(outdir.comp, paste(comp, '_heatmap_unfiltered_note_pident_', x, '.pdf', sep='')))
            hclust.outfiles.note2 <- sapply(hclust.methods, function(x)
                file.path(outdir.comp, paste(comp, '_heatmap_unfiltered_note_length_', x, '.pdf', sep='')))
            names(hclust.outfiles.note1) <- hclust.methods
            names(hclust.outfiles.note2) <- hclust.methods
            # file size depends on number of samples
            pdf.size <- 7 + (nsamp * 0.05)
            for (m in hclust.methods){
                # note version with % identity
                pdf(hclust.outfiles.note1[m], height = pdf.size, width=pdf.size)
                heatmap.2(ani.mat.plot, col = viridis(128), trace='none',
                          na.color = 'grey50', margins = c(10,10),
                          main=paste(cluster.name, comp),
                          cexRow = 0.5, cexCol = 0.5,
                          distfun = function(x) as.dist(100-ani.mat),
                          hclustfun = function(x) hclust(x, method=m),
                          key.xlab = comp, key.title = '',
                          lmat=lmat, lhei=lhei, lwid=lwid, key.par = list(mar=c(4,2.5,1,0.5)),
                          cellnote = round(ani.mat, 3), notecex = 0.4, notecol = 'magenta')
                          # cellnote = formatC(100-ani.mat, digits = 3, format='e'), notecex = 0.4, notecol = 'magenta')
                dev.off()
                # note version with length of alignment
                pdf(hclust.outfiles.note2[m], height = pdf.size, width=pdf.size)
                heatmap.2(ani.mat.plot, col = viridis(128), trace='none',
                          na.color = 'grey50', margins = c(10,10),
                          main=paste(cluster.name, comp),
                          cexRow = 0.5, cexCol = 0.5,
                          distfun = function(x) as.dist(100-ani.mat),
                          hclustfun = function(x) hclust(x, method=m),
                          key.xlab = comp, key.title = '',
                          lmat=lmat, lhei=lhei, lwid=lwid, key.par = list(mar=c(4,2.5,1,0.5)),
                          cellnote = round(length.mat, 2), notecex = 0.4, notecol = 'magenta')
                dev.off()
            }
        }
    }
}

# remove pesky rplots.pdf
if(file.exists('Rplots.pdf')){
    file.remove('Rplots.pdf')
}




# filtered versions not implemented yet
if (F){
# filtered version
if (all(ins.filt[!is.na(ins.filt)] == 0)){
    pdf(file.path(outdir.comp, 'ani_heatmap_complete.pdf'))
    plot(1,1, col='white')
    text(1,1, "Not enough nonzero samples for heatmap")
} else {
    hclust.methods <- c('ward.D2','single','complete')
    hclust.outfiles <- sapply(hclust.methods, function(x)
        file.path(outdir.comp, paste('ani_heatmap_', x, '.pdf', sep='')))
    names(hclust.outfiles) <- hclust.methods

    # file size depends on number of samples
    pdf.size <- 10 + (nsamp * 0.05)
    for (m in hclust.methods){
        pdf(hclust.outfiles[m], height = pdf.size, width=pdf.size)
        heatmap.2(ins.filt.plot, col = viridis(128), trace='none',
                  na.color = 'grey50', margins = c(10,10),
                  main='nucmer assembled genome average identity',
                  cexRow = 0.5, cexCol = 0.5,
                  distfun = function(x) as.dist(100-ins.filt),
                  key.xlab = '% identity', key.title = '',
                  hclustfun = function(x) hclust(x, method=m))
        dev.off()
    }
}


# note version, using filtered
if (any(ani.mat[!is.na(ani.mat)] > 0)){
    hclust.methods <- c('ward.D2','single','complete')
    hclust.outfiles.note1 <- sapply(hclust.methods, function(x)
        file.path(outdir.comp, paste('ani_heatmap_note_pident_', x, '.pdf', sep='')))
    hclust.outfiles.note2 <- sapply(hclust.methods, function(x)
        file.path(outdir.comp, paste('ani_heatmap_note_length_', x, '.pdf', sep='')))
    names(hclust.outfiles.note1) <- hclust.methods
    names(hclust.outfiles.note2) <- hclust.methods
    # file size depends on number of samples
    pdf.size <- 10 + (nsamp * 0.05)
    for (m in hclust.methods){
        # note version with % identity
        pdf(hclust.outfiles.note1[m], height = pdf.size, width=pdf.size)
        heatmap.2(ani.mat.filt.plot, col = viridis(128), trace='none',
                  na.color = 'grey50', margins = c(10,10),
                  main='nucmer assembled genome average identity',
                  cexRow = 0.5, cexCol = 0.5,
                  distfun = function(x) as.dist(100-ani.mat.filt),
                  hclustfun = function(x) hclust(x, method=m),
                  key.xlab = '% identity', key.title = '',
                  cellnote = round(ani.mat.filt, 2), notecex = 0.4, notecol = 'magenta')
        dev.off()
        # note version with length of alignment
        pdf(hclust.outfiles.note2[m], height = pdf.size, width=pdf.size)
        heatmap.2(ani.mat.filt.plot, col = viridis(128), trace='none',
                  na.color = 'grey50', margins = c(10,10),
                  main='nucmer assembled genome average identity',
                  cexRow = 0.5, cexCol = 0.5,
                  distfun = function(x) as.dist(100-ani.mat.filt),
                  hclustfun = function(x) hclust(x, method=m),
                  key.xlab = '% identity', key.title = '',
                  cellnote = round(length.mat.filt/1000), notecex = 0.4, notecol = 'magenta')
        dev.off()
    }
}
}
