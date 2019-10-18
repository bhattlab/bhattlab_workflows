# functions for reading and processing snp data

# read a snp calling file from my pipeline
# and applies some corrections
# returns only fixed sites
read_snp_file_fixed <- function(snp.f){
  # tab file has all the variant  information we need
  snp.df.samp <- read.table(snp.f, sep='\t', quote='', header=T)
  # ensure at least one row
  if(nrow(snp.df.samp)==0){
    return(NA)
  }
  # any duplicated rows?
  snp.df.samp$POS_ALT <- paste(snp.df.samp$POS, snp.df.samp$ALT, sep='_')
  if (any(duplicated(snp.df.samp$POS_ALT))){
    dup.pos <- snp.df.samp$POS_ALT[duplicated(snp.df.samp$POS_ALT)]
    message(paste('found duplicate positions:', length(dup.pos)))
    # print(snp.df.samp[snp.df.samp$POS_ALT %in% dup.pos, ])
    all.dup.df <- snp.df.samp[snp.df.samp$POS_ALT %in% dup.pos, ]
    snp.df.samp.nodup <- snp.df.samp[!(snp.df.samp$POS_ALT %in% dup.pos), ]
    dup.df.list <- list()
    for (pos in dup.pos){
      this.df <- snp.df.samp[snp.df.samp$POS_ALT == pos, ]
      if(nrow(this.df)> 2){
        print("not sure what to do if >2 duplicates...")
        print("cheating!")
        print(this.df)
        this.df <- this.df[1:2,]
      }
      # do we merge
      if (sum(this.df$AO) <= this.df$DP[1]){
        # merge
        this.df[1, 'AO'] <- sum(this.df$AO)
        this.df <- this.df[1,]
      } else if ((sum(this.df$AO) > this.df$DP[1]) & (this.df$AO[1] ==this.df$AO[2])){
        # duplicate all the way through, elim row
        this.df <- this.df[1,]
      } else {
        warning(paste('Bad duplicate case at', pos))
        this.df <- this.df[2,]
        
      }
      dup.df.list[[pos]] <- this.df
    }
    # message('After deduplication and adding of variants: ')
    # print(snp.df.samp[snp.df.samp$POS_ALT %in% dup.pos, ])
    fixed.dup.df <- do.call(rbind, dup.df.list)
    # print(fixed.dup.df)
    snp.df.samp <- rbind(snp.df.samp.nodup, fixed.dup.df)
    snp.df.samp <- snp.df.samp[order(snp.df.samp$POS),]
    # print(snp.df.samp[snp.df.samp$POS==12357, ])
  }
  
  # areas where there is multiple alleles present
  multi.pos <- grep(',', snp.df.samp$ALT)
  if (length(multi.pos)> 0){
    message('found multiallelic positions')
    print(snp.df.samp[multi.pos,])
    # remove these to a separate dataframe to be delt with later
    multi.df.samp <- snp.df.samp[multi.pos,]
    snp.df.samp <- snp.df.samp[-1*multi.pos, ]
    # multi.df <- rbind(multi.df, multi.df.samp)
  }
  
  # weird cases where theres a comma in AO
  multi.pos <- grep(',', snp.df.samp$AO)
  if (length(multi.pos)> 0){
    message('found multiallelic positions in AO?')
    print(snp.df.samp[multi.pos,])
    # set ao to the max?
    new.ao <- sapply(snp.df.samp[multi.pos,'AO'], function(x) as.numeric(max(strsplit(x,',')[[1]])))
    snp.df.samp[multi.pos,'AO'] <- new.ao
  }
  
  
  snp.df.samp$AO <- as.numeric(snp.df.samp$AO)
  snp.df.samp$RO <- as.numeric(snp.df.samp$RO)
  snp.df.samp$DP <- as.numeric(snp.df.samp$DP)
  snp.df.samp$POS <- as.numeric(snp.df.samp$POS)
  
  # filtering based on counts at each position
  # we're interested in multiallelic sites that have some 
  # significance of reads supporting them. 
  # if we say that must have at least 5 reads at a position
  # and if less than that you don't have a variant?
  # if its a variant that we dont have another entry
  # for (not confident that its multiallelic) then just assign 
  # it as a single variant
  multi.threshold <-5
  snp.df.samp$RO[snp.df.samp$RO < multi.threshold] <- 0
  snp.df.samp$RO[snp.df.samp$AO < multi.threshold] <- 0
  
  # called varaint depth
  called.depth.df <- aggregate(snp.df.samp[, c('RO', 'AO')], by=list(snp.df.samp$POS), sum)
  called.depth.df$all <- rowSums(called.depth.df[,2:3])
  called.depth.vec <- called.depth.df$all
  names(called.depth.vec) <- called.depth.df$Group.1
  
  snp.df.samp$DP.called <- sapply(snp.df.samp$POS, function(x) called.depth.vec[x])
  
  
  # eliminate anything when DP.called falls below 10
  snp.df.samp <- snp.df.samp[which(snp.df.samp$DP.called >=10), ]
  # calculate allelic fraction 
  snp.df.samp$AF.alt <- snp.df.samp$AO / snp.df.samp$DP.called
  snp.df.samp$AF.ref <- snp.df.samp$RO / snp.df.samp$DP.called
  
  
  # change format here
  # major allele and minor allele
  # "reference" no longer used as a term
  dup.pos <- snp.df.samp$POS[duplicated(snp.df.samp$POS)]
  head(snp.df.samp[snp.df.samp$POS %in% dup.pos,])
  
  snp.df.samp.new <- lapply(unique(snp.df.samp$POS), function(pos){
    this.df <- snp.df.samp[snp.df.samp$POS==pos, ]
    if(nrow(this.df) == 1){
      # nothing new
      alleles <- c(this.df$AF.ref, this.df$AF.alt)
      names(alleles) <- c(this.df$REF, this.df$ALT)
      alleles <- sort(alleles, decreasing = T)
      major.allele <-      new.df <- data.frame(SAMPLE=this.df$SAMPLE, POS=this.df$POS, REF=this.df$REF, DP=this.df$DP, DP.called=this.df$DP.called,
                                                major.allele = names(alleles)[1], minor.allele = names(alleles)[2],
                                                major.AF = alleles[1], minor.AF = alleles[2])
    } else if (nrow(this.df)==2){
      # new filtering
      if(any(this.df$AF.ref >0)){
        # print(this.df)
        alleles <- c(this.df$AF.alt, this.df$AF.ref)
        names(alleles) <- c(this.df$ALT, this.df$REF)
        alleles <- sort(alleles, decreasing = T)
        # stop('3 alleles including ref!')
      } else{
        alleles <- this.df$AF.alt
        names(alleles) <- this.df$ALT
        alleles <- sort(alleles, decreasing = T)
      }
      new.df <- data.frame(SAMPLE=this.df$SAMPLE[1], POS=this.df$POS[1], REF=this.df$REF[1], DP=this.df$DP[1], DP.called=this.df$DP.called[1],
                           major.allele = names(alleles)[1], minor.allele = names(alleles)[2],
                           major.AF = alleles[1], minor.AF = alleles[2])
    } else {
      print(this.df)
      stop('3 variants at position!')
    }
    return(new.df)
  })
  snp.df.samp.new <- do.call(rbind, snp.df.samp.new)
  
  return(snp.df.samp.new)
}

# read a depth file 
# from samtools depth
read_depth_file <- function(depth.f){
  if (file.size(depth.f)>0){
    depth.samp <-read.table(depth.f, sep='\t', quote='', header=F)[,3]
    return(depth.samp)
  } else{
    return(NA)
  }
}


# get number and fraction of snps between each 
# given a snp dataframe for sample1 and sample2, 
# depth vector for sample1 and sample2, 
# and the minimum depth to consider positions
get_snp_distance <- function(snp.df1, snp.df2, d1, d2, min.depth){
  # going to exclude sites where one sample is multiallelic
  exclude.s1 <- snp.df1[snp.df1$major.AF < 1, "POS"]
  exclude.s2 <- snp.df2[snp.df2$major.AF < 1, "POS"]
  exclude.both <- unique(c(exclude.s1, exclude.s2))
  # limit to fixes sites 
  snp.df1 <- snp.df1[snp.df1$major.AF == 1,]
  snp.df2 <- snp.df2[snp.df2$major.AF == 1,]
  # remove excluded
  snp.df1 <- snp.df1[!(snp.df1$POS %in% exclude.both),]
  snp.df2 <- snp.df2[!(snp.df2$POS %in% exclude.both),]
  
  # sites where each is covered well enough  
  common.depth <- which(d1>=min.depth & d2>=min.depth)
  # limit to common depth sites
  snp.df1 <- snp.df1[snp.df1$POS %in% common.depth,]
  snp.df2 <- snp.df2[snp.df2$POS %in% common.depth,]
  snp.df1$POS_major.allele <- paste(snp.df1$POS, snp.df1$major.allele, sep='_')
  snp.df2$POS_major.allele <- paste(snp.df2$POS, snp.df2$major.allele, sep='_')
  all.alts <- unique(c(snp.df1$POS_major.allele, snp.df2$POS_major.allele))
  alt.df <- data.frame(site=all.alts, 
                       s1 = all.alts %in% snp.df1$POS_major.allele,
                       s2 = all.alts %in% snp.df2$POS_major.allele)
  snp.distance <- sum(alt.df$s1 != alt.df$s2)
  return(snp.distance)
}

# number of sites covered well enough in each genome
get_common_depth <- function(d1, d2, min.depth){
  common.depth <- which(d1>=min.depth & d2>=min.depth)
  return(length(common.depth))
}


