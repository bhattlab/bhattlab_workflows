# get the instrain results from this cluster
parse_instrain <- function(ins.f, min.frac.compared) {
    # read instrain df
    ins.orig <- read.table(ins.f, sep='\t', header=T, quote='',
                           colClasses = c(rep('character', 3),
                                          rep('numeric', 7)))
    if(nrow(ins.orig) >0){
        # first remove rows with no coverage
        ins <- ins.orig[!is.na(ins.orig$coverage_overlap),]
        # remove low coverage
        ins <- ins[ins$percent_compared >= min.frac.compared, ]
        # order by ani
        ins <- ins[order(ins$popANI, decreasing = T),]
        # fix names by removing .sorted and .bam
        ins$name1 <- gsub('\\.bam', '', ins$name1)
        ins$name2 <- gsub('\\.bam', '', ins$name2)
        ins$name1 <- gsub('\\.sorted', '', ins$name1)
        ins$name2 <- gsub('\\.sorted', '', ins$name2)

        ins$pt1 <- sapply(ins$name1, function(x) strsplit(x, split='_')[[1]][1])
        ins$pt2 <- sapply(ins$name2, function(x) strsplit(x, split='_')[[1]][1])
        ins$same.pt <- ins$pt1==ins$pt2

        return(ins)
    } else {
        return(data.frame())
    }
}

# convert instrain mappig to a matrix of pairwise
ins_to_matrix_list <- function(ins.raw, ani.thresh=0, use.col='popANI'){
    if(ani.thresh > 1){
        warning("ani.thresh should be a fraction of 1. Dividing your value by 100")
        ani.thresh <- ani.thresh / 100
    }
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

# annotate an instrain df with the matrix of swapping
# read from a file into a matrix
annotate_with_swap_mat <- function(ins, swap.mat.f, swap.thresh=0.3){
    swap.mat <- read.table(swap.mat.f, sep='\t', quote='', header=T)
    rownames(swap.mat) <- make.names(rownames(swap.mat))

    # add this swap mat to the compare df
    a <- apply(ins, 1, function(x) {
        swap.mat[make.names(x['name1']), make.names(x['name2'])]
    })

    if(is.null(a)) {
        ins$swap.frac <- NA
    } else {
        a[sapply(a, is.null)] <- NA
        ins$swap.frac <- unlist(a)
    }

    ins$swap.frac[is.na(ins$swap.frac)] <- 0
    # check the top swappers, it's not many but they are significant
    # ins[order(ins$swap.frac, decreasing = T), ][1:10,]

    # remove cases where different patients and there was > 30% swapping
    ins$swapping.suspected <- FALSE
    ins[ins$swap.frac > swap.thresh & !ins$same.pt, 'swapping.suspected'] <- TRUE
    # could also just set all these values to NA
    # ins[ins$swap.frac > swap.thresh & !ins$same.pt, "popANI"] <- NA
    # ins[ins$swap.frac > swap.thresh & !ins$same.pt, "conANI"] <- NA
    return(ins)

}

# check for barcode swapping in a matrix
check.swap <- function(swap.mat, name){
    report <- swap.mat[make.names(name),,drop=F]
    report <- unlist(c(report[which(!is.na(report))]))
    report <- sort(report, decreasing = T)
    return(report)
}
