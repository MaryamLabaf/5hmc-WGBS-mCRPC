# input: a list of gRange peaks abd the target range
# output <- a data frame of overlapped peaks with target range
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
  filenames = list.files(dir_methylseqR, pattern=paste(chr,'UMRsLMRs.txt',sep='_'), full.names=F)
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

tpm_exp <- function(dir_tpm,ensembl_id, gene2find){
  file_tpm = list.files(dir_tpm, pattern="tpm")
  stopifnot(length(file_tpm)==1)
  filenamefull = paste(dir_tpm,file_tpm,sep='/')
  tpm = read.delim(filenamefull,sep='',row.names=1,header=TRUE,stringsAsFactors=FALSE,check.names=F)
  geneexpr <- tpm[ensembl_id,] %>% t %>% as.data.frame() %>% tibble::rownames_to_column("sample_id")
  colnames(geneexpr)[2] <- gene2find
  geneexpr$sample <- sub("^([^-]*-[^-]*-[^-]*).*", "\\1", geneexpr$sample_id)
  geneexpr$Type <- gsub("-RNA", "", sub("^(?:[^-]*-){3}", "\\1", geneexpr$sample_id))
  return(geneexpr)
}

get_peak_enrich <- function(CHR, START, END, SAMPLES, PEAK_OBJ) {
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