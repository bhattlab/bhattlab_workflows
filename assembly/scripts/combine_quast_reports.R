# merge quast results from many files and pick a common set of columns
options(stringsAsFactors = F)

read_quast_report <- function(f, name, keep.columns=c('Assembly', 'contigs')){
  df <- tryCatch({
    read.table(f, sep='\t', quote='', header=T, comment.char = '')  
  }, error = function(e){
    return(data.frame(matrix(nrow=1, ncol=length(keep.columns),dimnames = list('', keep.columns))))
  })
  # print(head(df))
  keep.columns.mod <- keep.columns[keep.columns %in% colnames(df)]  
  keep.columns.mod <- keep.columns.mod[1:2]
  df <- df[, keep.columns.mod]
  colnames(df)[2] <- name
  return(df)
}

sample.names <- snakemake@params[['sample_names']]
assembly.dir <- snakemake@params[['assembly_dir']]

result.files <- sapply(sample.names, function(x) file.path(assembly.dir, x, 'quast/report.tsv'))

if(!all(file.exists(result.files))) {exit('Cant find all result files!')}

df.list <- lapply(sample.names, function(x) read_quast_report(result.files[x], x, keep.columns=c("Assembly", "contigs", paste0(x, ".contigs"), gsub('-', '_', paste0(x, ".contigs")))))

merge.colname <- 'Assembly'
merge.temp <- suppressWarnings(Reduce(function(x,y) merge(x, y, all=TRUE, by=merge.colname, sort=F), df.list))
merge.temp <- as.data.frame(t(merge.temp))
colnames(merge.temp) <- merge.temp[1,]
merge.temp$sample <- rownames(merge.temp)
# print(head(merge.temp))

column.order <- c("sample", "N50", "N75", "# contigs", "Largest contig", "Total length", "# contigs (>= 0 bp)", "# contigs (>= 1000 bp)", "# contigs (>= 5000 bp)", "# contigs (>= 10000 bp)", "# contigs (>= 25000 bp)", "# contigs (>= 50000 bp)", "Total length (>= 0 bp)", "Total length (>= 1000 bp)", "Total length (>= 5000 bp)", "Total length (>= 10000 bp)", "Total length (>= 25000 bp)", "Total length (>= 50000 bp)", "Reference length", "GC (%)", "Reference GC (%)", "L50", "L75", "# total reads", "# left", "# right", "Mapped (%)", "Reference mapped (%)", "Properly paired (%)", "Reference properly paired (%)", "Avg. coverage depth", "Reference avg. coverage depth", "Coverage >= 1x (%)", "Reference coverage >= 1x (%)", "# misassemblies", "# misassembled contigs", "Misassembled contigs length", "# local misassemblies", "# scaffold gap ext. mis.", "# scaffold gap loc. mis.", "# structural variations", "# unaligned mis. contigs", "# unaligned contigs", "Unaligned length", "Genome fraction (%)", "Duplication ratio", "# N's per 100 kbp", "# mismatches per 100 kbp", "# indels per 100 kbp", "Largest alignment", "Total aligned length", "NA50", "NGA50", "LA50", "NA75", "LA75")
column.order <- column.order[column.order %in% colnames(merge.temp)]
merge.temp <- merge.temp[, column.order]
merge.temp <- merge.temp[c('Assembly', sort(merge.temp[2:nrow(merge.temp), 'sample'])), ]
write.table(merge.temp, snakemake@output[[1]], sep = '\t', quote=F, row.names = F, col.names = F)
