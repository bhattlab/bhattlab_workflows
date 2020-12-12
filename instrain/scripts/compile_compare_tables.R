# quick script to run after instrain compare to compile the resutls from all clusters
# and calculate some informative stats on them all
# takes in a directory and a two column tsv with sample patient metadata
options(stringsAsFactors = F)
# OPTIONS
# minimum fraction of the genome compared between samples
# I leave this at a relatively lax 30% here becuase you can always filter more
min.frac.compared <- 0.3
### TESTING ARGS
if (F){
    sample.groups.f <- '~/bhatt_local/tx88/processing_temp/basksg/sample_groups.tsv'
    instrain.compare.filtered.f <- '~/bhatt_local/tx88/processing_temp/basksg/instrain_compare_filtered/'
}

# snakemake args
outdir <- snakemake@params[['outdir']]
sample.groups.f <- snakemake@params[['sample_groups']]
# source the instrain tools for parsing
instrain.tools.f <- file.path(snakemake@scriptdir, 'instrain_tools.R')
source(instrain.tools.f)
# output the tables in this directory
outf1 <- file.path(outdir, 'instrain_compare_compiled.tsv')
outf2 <- file.path(outdir, 'instrain_compare_compiled_diffpt.tsv')

# get a list of clusters from the folders within
cluster.folders <- list.dirs(outdir, full.names = F,recursive = F)
# result files within each
res.files <- sapply(cluster.folders, function(x) file.path(outdir, x, 'output', paste0(x, '_genomeWide_compare.tsv')))
#all(file.exists(res.files))

# load them all
res.list <- lapply(res.files, function(x) parse_instrain(x, min.frac.compared))
print('all loaded')
# add cluster annotation
for (c in names(res.list)){
    if(nrow(res.list[[c]]) >0){
        res.list[[c]]$cluster <- c
    }
}
res.df <- do.call(rbind, res.list)

# if a sample groups annotation is given, we can add that to the dataframe
# otherwise the pipeline will guess based on the sample name
if (file.exists(sample.groups.f)){
    sample.groups <- read.table(sample.groups.f, sep='\t', quote='', header=F)
    colnames(sample.groups) <- c('sample', 'group')
    rownames(sample.groups) <- sample.groups$sample
    res.df$pt1 <- sample.groups[res.df$name1, "group"]
    res.df$pt2 <- sample.groups[res.df$name2, "group"]
    res.df$same.pt <- res.df$pt1 == res.df$pt2
}

# sort and filter to separate pts
res.df <- res.df[order(res.df$popANI, decreasing = T), ]
res.df.diffpt <- res.df[!res.df$same.pt, ]

# write to the designated output files
write.table(res.df, outf1, sep='\t', quote=F, row.names = F, col.names = T)
write.table(res.df.diffpt, outf2, sep='\t', quote=F, row.names = F, col.names = T)
