# annotate the ARG table with information on bin, organism, etc
options(stringsAsFactors = F)

# short read data
bin.folder <- '~/scg/bmt_microbiome_data/binning_das_tool/'
rgi.folder <- '~/scg/bmt_microbiome_data/arg_detection/rgi/'
outfolder <- '~/scg/bmt_microbiome_data/arg_detection/rgi_improved/'
type <- 'SR'

# also do for 10x
# bin.folder <- '~/scg/tx88/binning_das_tool/'
# rgi.folder <- '~/scg/tx88/arg_detection/rgi/'
# outfolder <- '~/scg/tx88/arg_detection/rgi_improved/'
# type <- '10x'

sample.list <- sapply(list.files(rgi.folder, pattern='txt'), function(x) strsplit(x, split='\\.')[[1]][1])
# sample.list <- c('p5160_2016-07-19_rep2', 'p6161_2016-10-05', 'p6929_2018-07-16', 'p7042_2016-12-27')

for (s in sample.list){
    print(s)
    bin.final.f <- file.path(bin.folder, s, 'final', paste(s, '.tsv', sep=''))
    contig2bin.f <- file.path(bin.folder, s, "DAS_tool/_DASTool_scaffolds2bin.txt")
    rgi.f <- file.path(rgi.folder, paste(s, '.rgi.txt', sep=''))

    rgi <- read.table(rgi.f, sep='\t', quote='', header=T, fill=T, comment.char = '')
    remove.columns <- c('Predicted_DNA','Predicted_Protein','CARD_Protein_Sequence')
    rgi <- rgi[,!(colnames(rgi) %in% remove.columns)]
    # remove last underscore from contig
    rgi$Contig.match <- sapply(rgi$Contig, function(x) {
        a <- strsplit(x, split='_')[[1]]
        paste.min <- min(4, length(a)-1)
        paste(a[1:paste.min], collapse='_')
    })

    bin.final <- read.table(bin.final.f, sep='\t', quote='', header=T, fill=T, comment.char = '')

    if (file.exists(contig2bin.f)){
        contig2bin.df <- read.table(contig2bin.f, sep='\t', quote='', header=F, colClasses = 'character')
        contig2bin <- contig2bin.df[,2]
        if (type=='SR'){
            names(contig2bin) <- sapply(contig2bin.df[,1], function(x) {
                a <- strsplit(x, split='_')[[1]]
                paste.min <- min(4, length(a)-1)
                paste(a[1:paste.min], collapse='_')
            })
        } else {
            names(contig2bin) <- contig2bin.df[,1]
        }
    } else {
        print(contig2bin.f)
        contig2bin <- rep('unbinned', times=length(unique(rgi$Contig.match)))
        names(contig2bin) <- unique(rgi$Contig.match)
    }


    if (!all(rgi$Contig.match %in% names(contig2bin))){
        unmatching <- rgi$Contig.match[!(rgi$Contig.match %in% names(contig2bin))]
        c2bnew <- rep('unbinned', times=length(unmatching))
        names(c2bnew) <- unmatching
        contig2bin <- c(contig2bin, c2bnew)
    }

    rgi$bin <- contig2bin[rgi$Contig.match]
    rgi$organism <- sapply(rgi$bin, function(x) {
        if (x %in% bin.final$Bin){
            bin.final[bin.final$Bin==x,'lca_species' ]
        } else {
            NA
        }
    })
    rgi$bin.quality.call <- sapply(rgi$bin, function(x) {
        if (x %in% bin.final$Bin){
            bin.final[bin.final$Bin==x,"bin.quality.call" ]
        } else {
            NA
        }
    })
    rgi$sample <- paste(s, type, sep='_')
    # reorder columns
    col.order <- c(27, 23:26, 3:22, 1,2)
    rgi.new <- rgi[, col.order]
    rgi.new <- rgi.new[order(rgi.new$bin, rgi.new$Cut_Off),]
    outf <- file.path(outfolder, paste(s, '.rgi.txt', sep=''))

    write.table(rgi.new, outf, sep='\t', quote=F, row.names = F, col.names = T)
}

# combine files into a single large table
# head -n 1 rgi_improved/bas_bb_spri.rgi.txt > rgi_improved_allsamples.txt
# for i in rgi_improved/*.rgi.txt; do echo $i; tail -n +2 $i >> rgi_improved_allsamples.txt; done
# for i in ~/bmt_microbiome_data/arg_detection/rgi_improved/*.rgi.txt; do echo $i; tail -n +2 $i >> rgi_improved_allsamples.txt; done
