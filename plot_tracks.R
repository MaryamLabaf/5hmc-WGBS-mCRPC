use2=FALSE
rangeup = 850000
rangedown = 100000
gene2find = 'AR'

rowensembl <- ensembl2sym$gene_name==gene2find
ensemblid <- rownames(ensembl2sym)[rowensembl]
ensembl_id <- ensembl2sym[ensemblid,"gene_ensembl"]
stopifnot(sum(rowensembl)==1)
ensembl2sym[ensemblid,]

chr <- ensembl2sym[rowensembl,'Chr']
gene_start <- ensembl2sym[rowensembl,'Start']
gene_end <- ensembl2sym[rowensembl,'End']

start <- round(gene_start-rangeup,-3)
end <- round(gene_end+rangedown,-3)
stopifnot(start<end)

grange <- GRanges(chr, IRanges(start, end))

name <- paste(gene2find,'loci',sep='_')
print(paste(name, ' - ', chr, ':', start, '-', end, sep=''))

filename <- paste(outputdir, name,'.pdf',sep='')
########################################################################
# Load in various data 
########################################################################
bedgene <- exons[,c(5,6,7,3,1,4)]
bedgene$strand <- as.numeric(bedgene$strand=='+')
bedgene[bedgene$strand==0,'strand'] <- -1

# Load 5hmC-peaks 
fivehmcslice <- fivehmcslice_fun(peaks_5hmc, grange)
# Add group 
length(peaks_5hmc)
dim(meta_data.flt)
meta_data.flt$sample_id <- toupper(meta_data.flt$sample_id)
common_samples <- intersect(fivehmcslice$sample, meta_data.flt$sample_id)
fivehmcslice <- fivehmcslice %>% left_join(meta_data.flt, by = c("sample" = "sample_id"))
fivehmcslice %>% dplyr::distinct(sample, .keep_all = T) %>% group_by(category) %>% dplyr::summarise(n = n())
fivehmcslice <- fivehmcslice %>% 
  dplyr::mutate(subClass = case_when(sample %in% samples_normal ~ "Benign",
                                     sample %in% samples_localized_tissue ~ "Loc-PCa",
                                     sample %in% samples_wcdt_tissue_nt ~ "N-Adjacent-tissue",
                                     sample %in% tscnc_pids ~ "mCRPC-t-SCNC",
                                     sample %in% union(samples_mcrpc_uniq_all_data, toupper(samples_wcdt_cfdna)) ~ "mCRPC-adeno",
                                     sample %in% toupper(samples_ubc_cfdna) ~ "mCRPC-cfDNA-befor1stARSi"))

fivehmcslice <- fivehmcslice %>% 
  dplyr::mutate(color = case_when(subClass == "Benign" ~ "steelblue3",
                                  subClass == "Loc-PCa" ~ "darkgoldenrod3",
                                  subClass == "N-Adjacent-tissue" ~ "green4", 
                                  subClass == "mCRPC-t-SCNC" ~ "gray2",
                                  subClass == "mCRPC-adeno" ~ "firebrick3",
                                  subClass == "mCRPC-cfDNA-befor1stARSi" ~ "maroon2"))

fivehmcslice %>% dplyr::distinct(sample, .keep_all = T) %>% group_by(subClass) %>% dplyr::summarise(n = n())
unique(fivehmcslice$sample[which(is.na(fivehmcslice$subClass))])
length(unique(fivehmcslice$sample))
fivehmcslice <- fivehmcslice %>% dplyr::filter(!is.na(subClass)) #Remove the NAs from the subClass samples

fivehmcslice$start <- pmax(fivehmcslice$start,start+1)
fivehmcslice$end <- pmin(fivehmcslice$end,end-1)
fivehmcslice$score <- 0
fivehmcslice$strand <- "."
colnames(fivehmcslice)[1] <- "chrom"
fivehmcslice.df <- fivehmcslice %>% dplyr::select(chrom, start, end, sample, score, strand, color)
rows2keep <- fivehmcslice.df$sample %in% c(samples_normal,
                                           samples_localized_tissue,
                                           samples_wcdt_tissue_nt,
                                           tscnc_pids,
                                           union(samples_mcrpc_uniq_all_data,toupper(samples_wcdt_cfdna)),
                                           toupper(samples_ubc_cfdna))
fivehmcslice.df <- fivehmcslice.df[rows2keep,]
fivehmcslice.df %>% group_by(color) %>% dplyr::summarise(n= n())
roworder <- as.numeric(factor(fivehmcslice.df$sample,levels=
                                c(samples_normal,
                                  samples_localized_tissue,
                                  samples_wcdt_tissue_nt,
                                  tscnc_pids,
                                  setdiff(union(samples_mcrpc_uniq_all_data,toupper(samples_wcdt_cfdna)),tscnc_pids),
                                  toupper(samples_ubc_cfdna))))

unique(roworder)

# Load WCDT-UMRs
hmrs <- hmr_peaks_fun(chr,dir_methylseqR = dir_wcdt_UMRLMR)
hmrslice <- hmr_fun(hmrs)
common_samples <- intersect(hmrslice$sample, toupper(sample_ids_wgbs))
#length(unique(hmrslice$sample))
hmrslice <- hmrslice %>% dplyr::filter(!sample %in% samples.rm) %>%
  left_join(meta_data.flt, by = c("sample" = "sample_id")) %>%
  dplyr::left_join(WGBS_sup %>% dplyr::select(Sample, Type), by = c("sample" = "Sample"))
hmrslice <- hmrslice %>% 
  dplyr::mutate(subClass = case_when(!is.na(Type) ~ Type,
                                     sample %in% toupper(samples_upitt_benign_adjacent_prostate) ~ "Benign",
                                     sample %in% toupper(samples_upitt_localized_tumor) ~ "Loc-PCa",
                                     sample %in% toupper(samples_wcdt_tissue_nt) ~ "N-Adjacent-tissue",
                                     str_starts(sample, "NT-") ~ "N-Adjacent-tissue",
                                     .default = NA))

hmrslice$subClass[which(hmrslice$subClass == "mCRPC")] <- "mCRPC-adeno"
hmrslice$subClass[which(hmrslice$subClass == "tSCNC")] <- "mCRPC-tSCNC"
hmrslice <- hmrslice %>% 
  dplyr::mutate(color = case_when(subClass == "Benign" ~ "steelblue3",
                                  subClass == "Loc-PCa" ~ "darkgoldenrod3",
                                  subClass == "N-Adjacent-tissue" ~ "green4", 
                                  subClass == "mCRPC-tSCNC" ~ "gray2",
                                  subClass == "mCRPC-adeno" ~ "firebrick3",
                                  subClass == "mCRPC (CMP)" ~ "tomato"))

hmrslice %>% dplyr::distinct(sample, .keep_all = T) %>%
  group_by(subClass) %>% dplyr::summarise(n())
unique(hmrslice$sample[which(is.na(hmrslice$subClass))])
length(unique(hmrslice$sample))
hmrslice <- hmrslice %>% dplyr::filter(!is.na(subClass)) #Remove the NAs from the subClass samples

hmrslice$start <- pmax(hmrslice$start,start+1)
hmrslice$end <- pmin(hmrslice$end,end-1)

hmrslice$score <- 0
hmrslice$strand <- "."
colnames(hmrslice)[1] <- "chrom"
hmrslice.df <- hmrslice %>% 
  dplyr::select(chrom,start,end,sample, score,strand,color)

loc_ids <- unique(hmrslice$sample[which(hmrslice$subClass == "Loc-PCa")])
Benign_ids <- unique(hmrslice$sample[which(hmrslice$subClass == "Benign")])
N_adjacent_ids <- unique(hmrslice$sample[which(hmrslice$subClass == "N-Adjacent-tissue")])
mCRPC_ids <- unique(hmrslice$sample[which(hmrslice$subClass == "mCRPC-adeno")])
mCRPC_CMP_ids <- unique(hmrslice$sample[which(hmrslice$subClass == "mCRPC (CMP)")])
tSCNC_ids <- unique(hmrslice$sample[which(hmrslice$subClass == "mCRPC-tSCNC")])

hmr_roworder <- as.numeric(factor(hmrslice.df$sample,levels=
                                    c(setdiff(mCRPC_ids,tSCNC_ids), mCRPC_CMP_ids, 
                                      tSCNC_ids,
                                      loc_ids,N_adjacent_ids,Benign_ids)))
unique(hmr_roworder)

## Load UMRs from PSMA samples 
hmrs_psma <- hmr_peaks_fun(chr,dir_methylseqR = dir_bidmc_UMRLMR)
hmrslice_psma <- hmr_fun(hmrs_psma)
hmrslice_psma <- hmrslice_psma %>% 
  dplyr::mutate(color = case_when(sample == "G23" ~ "darkred",
                                  sample == "I24" ~ "coral3",
                                  sample == "G24" ~ "deepskyblue4",
                                  sample == "I23" ~ "blueviolet"))

hmrslice_psma$start <- pmax(hmrslice_psma$start,start+1)
hmrslice_psma$end <- pmin(hmrslice_psma$end,end-1)
hmrslice_psma$score <- 0
hmrslice_psma$strand <- "."
colnames(hmrslice_psma)[1] <- "chrom"
hmrslice_psma.df <- hmrslice_psma %>% 
  dplyr::select(chrom,start,end,sample, score,strand,color)

hmr_psma_roworder <- as.numeric(factor(hmrslice_psma.df$sample, 
                                       levels = c("G23","I24","G24","I23")))

# CpG islands
cpg <- as.data.frame(cpg_islans)[1:3]
cpg$row <- 4

## correlate 5HMC peaks
geneexpr <- tpm_exp(dir_tpm, ensembl_id = ensembl_id, gene2find = gene2find)
rhmrlevel_granges <- makeGRangesFromDataFrame(rhmrlevel, keep.extra.columns = T)
rhmrtrack <- rhmrlevel_granges[countOverlaps(rhmrlevel_granges,grange)>0]
samples_5hmc_expr <- intersect(intersect(samples_mcrpc_uniq_all_data, geneexpr$sample),names(peaks_5hmc))
geneexpr_5hmc <- geneexpr %>% dplyr::filter(sample %in% samples_5hmc_expr) %>% 
  arrange(Type == 'NS-T') %>% distinct(sample, .keep_all = TRUE) %>% 
  dplyr::select(sample, gene2find) %>% 
  tibble::column_to_rownames("sample") %>% t
X <- as.numeric(geneexpr_5hmc)
names(X) <- colnames(geneexpr_5hmc)

samples_wgbs_expr <- intersect(intersect(geneexpr$sample,samples_mcrpc_uniq_all_data),colnames(rhmrlevel))
geneexpr_wgbs <- geneexpr %>% dplyr::filter(sample %in% samples_wgbs_expr) %>% 
  arrange(Type == 'NS-T') %>%  distinct(sample, .keep_all = TRUE) %>%
  dplyr::select(sample, gene2find) %>% 
  tibble::column_to_rownames("sample") %>% t
Y <- as.numeric(geneexpr_wgbs)
names(Y) <- colnames(geneexpr_wgbs)

###########################
around <- 1000
peaks <- data.frame(rhmrtrack)[,1:3]
for(i in 1:length(rhmrtrack)) {
  print(i)
  starti <- peaks[i,'start']
  endi <- peaks[i,'end']
  startiaround <- peaks[i,'start']-around
  endiaround <- peaks[i,'end']+around
  
  hmci <- get_peak_enrich(chr,startiaround,endiaround,samples_5hmc_expr,peaks_5hmc)
  corobj <- cor.test(hmci,X[names(hmci)],method='spearman',na.rm=T)
  peaks[i,'hmc_p'] <- corobj$p.value
  peaks[i,'hmc_cor'] <- corobj$estimate
  
  wgbsi <- as.numeric(data.frame(rhmrtrack[i,samples_wgbs_expr])[-1:-5])
  names(wgbsi) <- samples_wgbs_expr
  corobj <- cor.test(	wgbsi,Y[samples_wgbs_expr],method='spearman',na.rm=T)
  peaks[i,'wgbs_p'] <- corobj$p.value
  peaks[i,'wgbs_cor'] <- corobj$estimate
}

peaks$expr_fdr <- p.adjust(peaks$hmc_p,method='fdr')
peaks$wgbs_fdr <- p.adjust(peaks$wgbs_p,method='fdr')
colnames(peaks)[1] <- 'chr'

#######################################################################
# PLOT TRACKS 
#######################################################################
colsglobal <- c("Benign"="steelblue3", 
                "Loc-PCa"="darkgoldenrod3",
                "N-Adjacent-tissue"="green4", 
                "mCRPC-tSCNC"="gray2",
                "mCRPC-adeno"="firebrick3",
                "mCRPC (CMP)"="tomato", 
                "mCRPC-cfDNA-befor1stARSi"="maroon2")

psam_col <- c("G23"="darkred","I24"="coral3","G24"="deepskyblue4","I23"="blueviolet")
numplots <- 7

pdf(file=filename, onefile=T, width=10, height=numplots/1.2)
matrows <- c(1,1,1,2,2,2,3,3,3,4,4,5,5,6,6,7)
layout(matrix(matrows, length(matrows), 1, byrow=T))

## Plot genes
par(mai=c(0.4,1,0.05,0.5))
exonrows <- bedgene$chr==chr & overlap(bedgene$start, bedgene$end, start, end)
if(sum(exonrows)>0) {
  plotGenes(bedgene[exonrows,],chr,start,end)
  labelgenome(chr,start,end,n=20,scale="Mb")
} else {
  plot.new()
}

par(mai=c(0.05,1,0.05,0.5))

## Plot 5hmc

plotBed(beddata=fivehmcslice.df,chrom=chr,chromstart = start,chromend = end,
        row='supplied',
        rownumber=roworder,
        color=fivehmcslice.df$color, border = T)
axis(side=2,las=2,tcl=.2,lwd = 0.5)
mtext('5hmc',side=2,line=4,cex=0.5)
abline(v = c(tss_start,tss_end), col="gray80", lwd=1, lty=1)

## Plot HMR-WCDT
plotBed(beddata=hmrslice.df,chrom=chr,chromstart = start,chromend = end,
        row='supplied',
        rownumber=hmr_roworder,
        color=hmrslice.df$color)

axis(side=2,las=2,tcl=.2)
mtext('HMR',side=2,line=4,cex=0.5)
abline(v = c(tss_start,tss_end), col="gray80", lwd=1, lty=1)

## Plot HMR-BIDMC samples
plotBed(beddata=hmrslice_psma.df,chrom=chr,chromstart = start,chromend = end,
        row='supplied',
        rownumber=hmr_psma_roworder,
        color=hmrslice_psma.df$color)

axis(side=2,las=2,tcl=.2)
mtext('HMR-PSMA-data',side=2,line=4,cex=0.5)
abline(v = c(tss_start,tss_end), col="gray80", lwd=1, lty=1)


## Plot rHMR
plotBed(beddata=rhmrslice.df,chrom=chr,chromstart = start,chromend = end,
        row="supplied",
        rownumber=rhmr_roworder,
        color=rhmrslice.df$color)

axis(side=2,las=2,tcl=.2)
mtext("rHMR",side=2,line=4,cex=0.5)

## plot CpGs

plotBed(cpg,chr,start,end,type="region",row='supplied',rownumber=1)
axis(side=2,las=2,tcl=.2)
mtext('CPG',side=2,line=4,cex=0.5)
abline(v = c(tss_start,tss_end), col="gray80", lwd=1, lty=1)

## Plot correlation with expression
cor_range <- c(-1,1)
colscor <- c('darkred','dodgerblue3')
plotBedgraph(peaks[,c('chr','start','end','hmc_cor')],chr,start,end,color=colscor[2],range=cor_range,transparency=0.6)
plotBedgraph(peaks[,c('chr','start','end','wgbs_cor')],chr,start,end,color=colscor[1],range=cor_range,transparency=0.6,overlay=T)
legend("topleft",inset=0,legend=c('WGBS','5hMC'),
       fill=opaque(colscor),border=colscor,cex=0.5)
axis(side=2,las=2,tcl=.2)
mtext('Rho',side=2,line=4,cex=0.5)
abline(h=0, col='black')

dev.off()