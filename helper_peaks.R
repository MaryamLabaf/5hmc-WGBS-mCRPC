#source('~/Desktop/DNA_methylation/scripts/helper_peaks.R')

#peakflank <- 0
#windowsize <- 1000

# data assumes chr1, pos1, chr2, pos2, svtype, sample_id
svwindow <- function(data, chr, start, end, window=windowsize, 
                     extra=peakflank, endsonly=T) {
  stopifnot((end-start)%%window==0)
  stopifnot(sum(is.na(data$chr1))==0)
  stopifnot(sum(colnames(data)!=c('chr1','pos1','chr2','pos2','svtype','sample_id'))==0)
  steps <- (end-start)/window
  output <- data.frame()
  for(i in 1:steps) {
    if(i %% 500 == 0) {
      print(paste(i, '/', steps))
    }
    left <- start+((i-1)*window)
    right <- left+window-1
    output[i,'chrom'] <- chr
    output[i,'start'] <- left
    output[i,'end'] <- right
    onechr <- is.na(data$chr2) | is.na(data$pos2)
    samechr <- data$chr1==data$chr2 & !onechr
    diffchr <- data$chr1!=data$chr2 & !onechr
    if(endsonly) {
      diffchr <- diffchr | samechr
      samechr <- rep(F,length(samechr))
    }
    left <- left-extra
    right <- right+extra
    rowsone <- onechr & inrange(data$pos1,left,right) & data$chr1==chr
    rowssame <- samechr & overlap(data$pos1,data$pos2,left,right) & data$chr1==chr
    rowsdiff <- diffchr &
      ((inrange(data$pos1,left,right) & data$chr1==chr) |
         (inrange(data$pos2,left,right) & data$chr2==chr))
    
    datarows <- rowsone | rowssame | rowsdiff
    stopifnot(sum(is.na(datarows))==0)
    output[i,'value'] <- length(unique(data[datarows,'sample_id']))
  }
  return(output)
}


# data assumes chrom, start, end, sample_id
slidingwindow <- function(data, chr, start, end, bysample=T, window=windowsize, 
                          extra=peakflank) {
  stopifnot((end-start)%%window==0)
  if(bysample) {
    stopifnot(sum(colnames(data)!=c('chrom','start','end','sample_id'))==0)
  } else {
    stopifnot(sum(colnames(data)!=c('chrom','start','end','value'))==0)
  }
  
  steps <- (end-start)/window
  output <- data.frame()
  for(i in 1:steps) {
    if(i %% 500 == 0) {
      print(paste(i, '/', steps))
    }
    left <- start+((i-1)*window)
    right <- left+window-1
    output[i,'chrom'] <- chr
    output[i,'start'] <- left
    output[i,'end'] <- right
    
    left <- left-extra
    right <- right+extra
    stopifnot((right-left+1) == (window+extra+extra))
    rows <- overlap(data$start,data$end,left,right) & data$chrom==chr
    stopifnot(sum(is.na(rows))==0)
    if(sum(rows)>0) {
      dataslice <- data[rows,]
      if(bysample) {
        output[i,'value'] <- length(unique(dataslice$sample_id))
      } else {
        #dataslice$start <- pmax(dataslice$start,start)
        #dataslice$end <- pmin(dataslice$start,end)
        #meansize <- mean(dataslice$end-dataslice$start+1)
        #dataslice$weight <- (dataslice$end-dataslice$start+1)/meansize
        #output[i,'value'] <- sum(dataslice$value*dataslice$weight,na.rm=T)
        output[i,'value'] <- sum(dataslice$value,na.rm=T)
      }
    } else {
      output[i,'value'] <- 0
    }
  }
  return(output)
}

findpeaks <- function(sw_bed, threshold, minwidth=10000, within=10000, minsize=3) {
  stopifnot(sum(colnames(sw_bed)!=c('chrom','start','end','value'))==0)
  results <- data.frame()
  inpeak <- F
  currentpeak <- 0
  
  if(max(sw_bed$value)<threshold) {
    return(NULL)
  }
  
  for(i in 1:dim(sw_bed)[1]) {
    stopifnot(!is.na(currentpeak))
    val <- sw_bed[i,'value']
    stopifnot(!is.na(val))
    if(val>=threshold) {
      if(!inpeak) { #new peak?
        inpeak <- T
        newpeakstart <- sw_bed[i,'start']
        if(currentpeak==0) { #first peak
          currentpeak <- currentpeak+1
          results[currentpeak,'start'] <- newpeakstart
          results[currentpeak,'max'] <- val
        } else if(newpeakstart-results[currentpeak,'end']>within) { #new peak!
          currentpeak <- currentpeak+1
          results[currentpeak,'start'] <- newpeakstart
          results[currentpeak,'max'] <- val
        } else {
          #keep current peak, and continue
          results[currentpeak,'end'] <- NA
        }
      } else { 
        #continue peak
        results[currentpeak,'max'] <- max(results[currentpeak,'max'],val)
      }
      if(i==dim(sw_bed)[1]) { #stop if last one
        results[currentpeak,'end'] <- sw_bed[i,'end']
        inpeak <- F
      }
    } else {
      if(inpeak) { #done with peak
        results[currentpeak,'end'] <- sw_bed[i-1,'end']
        inpeak <- F
      } else {
        #do nothing and keep going
      }
    }
    #print(results)
  }
  stopifnot(sum(is.na(results$start))==0)
  stopifnot(sum(is.na(results$end))==0)
  stopifnot(sum(is.na(results$max))==0)
  results$width <- results$end-results$start
  results <- results[results$width>=minwidth & results$max>=minsize,]
  return(results)
}

overlap <- function(x1, x2, y, z) {
  stopifnot(length(x1)==length(x2))
  if(length(x1)==1 && !is.na(x1) && !is.na(x2) && x1>x2) {
    x3 <- x2
    x2 <- x1
    x1 <- x3
  } else if(length(x1)>1) {
    rowswrong <- x1>x2
    rowswrong[is.na(rowswrong)] <- F
    x3 <- x2[rowswrong]
    x2[rowswrong] <- x1[rowswrong]
    x1[rowswrong] <- x3
  }
  stopifnot(y<=z)
  x1in <- x1>=y & x1<=z
  x2in <- x2>=y & x2<=z
  xencompasses <- x1<y & x2>z
  return(x1in | x2in | xencompasses)
}

#####
# Added MS
get_peak_enrich <- function(CHR, START, END, SAMPLES, PEAK_OBJ) {
  #samples2keep <- intersect(SAMPLES, names(PEAK_OBJ))
  theRange <- GRanges(CHR,IRanges(START,END))
  thePeaks <- PEAK_OBJ[SAMPLES]
  subPeaks <- list()
  for(i in 1:length(thePeaks)){
    rows2keep <- countOverlaps(thePeaks[[i]],theRange)>0
    subPeaks[names(thePeaks[i])] <- thePeaks[[i]][rows2keep]
  }
  averageFC <- lapply(subPeaks, FUN = function(x){ max(x$fold_enrichment) })
  return(unlist(averageFC))
}


fivehmcslice_fun <- function(peaks_5hmc,grange){
  all_5hmc <- list()
  for(i in 1:length(peaks_5hmc)) {
    samplename <- names(peaks_5hmc)[i]
    rows2keep <- countOverlaps(peaks_5hmc[[i]],grange)>0
    peaksoverlap <- peaks_5hmc[[i]][rows2keep]
    curtab <- data.frame(peaksoverlap)
    rowindex <- dim(curtab)[1]+1
    curtab[rowindex,1] <- chr
    curtab[rowindex,2] <- end-10
    curtab[rowindex,3] <- end
    curtab$sample <- toupper(samplename)
    cols_to_keep <- c("seqnames","start","end","neglog10q","sample")
    all_5hmc[[samplename]] <- curtab[,cols_to_keep]
  }
  
  fivehmcslice <- rbindlist(all_5hmc)[,c(1:3,5)] %>% as.data.frame()
  return(fivehmcslice)
}

hmr_peaks_fun <- function(chr, dir_methylseqR){
  #dir_methylseqR <- "WCDT_5hmC_additional_data/WGBS_additional_data/2019_05_15_WGBS_MethySeekR/"
  filenames = list.files(dir_methylseqR, pattern=paste(chr,'UMRsLMRs.txt',sep='_'), 
                         full.names=F)
  hmrs = list()
  for(curfile in filenames) {
    nametab = strsplit(curfile,'/',fixed=T)[[1]]
    samplename = strsplit(nametab[length(nametab)],'_',fixed=T)[[1]][1]
    curfile = paste(dir_methylseqR,curfile,sep='/')
    curtab = read.delim(file=curfile,sep='\t',header=TRUE,stringsAsFactors=F)
    curtab$sample = toupper(samplename)
    hmrs[[samplename]] = curtab %>% toGRanges()
  }
  return(hmrs)
}

hmr_fun <- function(hmrs) {
  all_hmrs <- list()
  for(i in 1:length(hmrs)) {
    samplename <- names(hmrs)[i]
    rows2keep <- countOverlaps(hmrs[[i]],grange)>0
    hmroverlap <- hmrs[[i]][rows2keep]
    curtab <- data.frame(hmroverlap)
    rowindex <- dim(curtab)[1]+1
    curtab[rowindex,1] <- chr
    curtab[rowindex,2] <- end-10
    curtab[rowindex,3] <- end
    curtab$sample <- toupper(samplename)
    all_hmrs[[samplename]] <- curtab
  }
  cols_to_keep <- c("seqnames","start","end","sample", "type")
  which(colnames(all_hmrs[[1]]) %in% cols_to_keep)
  hmrslice <- rbindlist(all_hmrs)[,c(1:3,11,6)] %>% as.data.frame()
  return(hmrslice)
}

tpm_exp <- function(ensembl_id, gene2find){
  dir_tpm <- "~/Desktop/DNA_methylation/Felix_Feng_WGBS_2019_05_15_WGBS/RNAseq/"
  file_tpm = list.files(dir_tpm, pattern="tpm")
  stopifnot(length(file_tpm)==1)
  filenamefull = paste(dir_tpm,file_tpm,sep='/')
  tpm = read.delim(filenamefull,sep='',row.names=1,header=TRUE,
                      stringsAsFactors=FALSE,check.names=F)
  #tpm_NS_T <- tpm %>% dplyr::select(matches("-NS-T"))
  #tpm_T <- tpm %>% dplyr::select(!matches("-NS-T"))
  #colnames(tpm_T) <- gsub("-T-RNA", "", colnames(tpm_T))
  #colnames(tpm_NS_T) <- gsub("-NS-T-RNA", "", colnames(tpm_NS_T))
  geneexpr <- tpm[ensembl_id,] %>% t %>% as.data.frame() %>% 
    tibble::rownames_to_column("sample_id")
  colnames(geneexpr)[2] <- gene2find
  geneexpr$sample <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1", geneexpr$sample_id)
  geneexpr$Type <- gsub("-RNA", "", sub("^(?:[^-]*-){3}", "\\1", geneexpr$sample_id))
  return(geneexpr)
}

match.idx = function(A, B, allow.multiple.B=F){
  # return dataframe of indices into A and B restricted to perfect matches
  # between A and B, where idx.A[i] == idx.B[i] for each i in matched pairs
  if( allow.multiple.B ){
    idx.B = which(B %in% A)
    idx.A = match(B[idx.B], A)
  }
  else{
    in.both = intersect(A,B)
    idx.A = match(in.both, A)
    idx.B = match(in.both, B)
  }
  C= data.frame(idx.A, idx.B)
  if( sum( A[ C$idx.A ] != B[ C$idx.B] )>0 )
    stop("ERROR! At least one in idx.A not the same as matched item in idx.B")
  C
}

#gets the mean methylation for a region given the path to the gene-level slice files
meanmethyl = function(rhmrs.df, start, end) {
  methyl <- rhmrlevel
  stopifnot(dim(methyl)[1]>0)
  methylrows = methyl$start >= start & methyl$end <= end
  print(sum(methylrows))
  means = colMeans(methyl[methylrows,-1:-3],na.rm=T)
  return(means)
}
#####
mean5hmc <- function(theRange, SAMPLES, PEAK_OBJ) {
  samples2keep <- intersect(SAMPLES, names(PEAK_OBJ))
  thePeaks <- PEAK_OBJ[samples2keep]
  subPeaks <- list()
  for(i in 1:length(thePeaks)){
    rows2keep <- countOverlaps(thePeaks[[i]],theRange)>0
    print(c(names(thePeaks[i]), sum(rows2keep)))
    subPeaks[names(thePeaks[i])] <- thePeaks[[i]][rows2keep]
  }
  averageFC <- lapply(subPeaks, FUN = function(x){ mean(x$fold_enrichment)})
  return(unlist(averageFC))
}