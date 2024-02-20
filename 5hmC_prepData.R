setwd("~/Desktop/DNA_methylation/Felix_Feng_5hmC_CancerRes_2022/")
#______________________________________________________________________#
ensembl2sym <- read.delim("WCDT_5hmC_additional_data/annotation_data/gene_models/gencode.v28.mapping.txt")
gene_model_hg38 <- read.delim("WCDT_5hmC_additional_data/annotation_data/gene_models/hg38_models.txt")
exons <- read.delim("WCDT_5hmC_additional_data/annotation_data/gene_models/gencode28_exons.txt")
#hypermethylated regions
wcdt_cfdna_incl <- readRDS("WCDT_5hmC_additional_data/hmc_secondary_data/peaks/peaks_wcdt_cfdna_incl.rds")
class(wcdt_cfdna_incl)
length(wcdt_cfdna_incl)
peaks_wcdt_cfdna_incl <- lapply(wcdt_cfdna_incl, 
                                ChIPQC:::GetGRanges, simple = TRUE)
head(peaks_wcdt_cfdna_incl[[1]])
localized_incl <- readRDS("WCDT_5hmC_additional_data/hmc_secondary_data/peaks/peaks_localized_incl.rds")
class(localized_incl)
length(localized_incl)
peaks_localized_incl <- lapply(localized_incl, 
                               ChIPQC:::GetGRanges, simple = TRUE)
head(peaks_localized_incl[[1]])

ubc_cfdna_incl <- readRDS("WCDT_5hmC_additional_data/hmc_secondary_data/peaks/peaks_ubc_cfdna_incl.rds")
class(ubc_cfdna_incl)
length(ubc_cfdna_incl)
peaks_ubc_cfdna_incl <- lapply(ubc_cfdna_incl, 
                               ChIPQC:::GetGRanges, simple = TRUE)
head(peaks_ubc_cfdna_incl[[1]])

wcdt_tissue_incl <- readRDS("WCDT_5hmC_additional_data/hmc_secondary_data/peaks/peaks_wcdt_tissue_incl.rds")
class(wcdt_tissue_incl)
length(wcdt_tissue_incl)
peaks_wcdt_tissue_incl <- lapply(wcdt_tissue_incl, 
                               ChIPQC:::GetGRanges, simple = TRUE)
head(peaks_wcdt_tissue_incl[[1]])

# Normal 5hmc
dir_normal_5hmc <- "Benign_PC_GSE144530_5hmc_peaks_hg19/"
N5hmc_filenames = list.files(dir_normal_5hmc, pattern='.narrowPeak', 
                             full.names=F)
Nr_5hmc.files <- dir("Benign_PC_GSE144530_5hmc_peaks_hg19/.", pattern="*.narrowPeak$")
Nr_5hmc.peaks <- sapply(paste0(dir_normal_5hmc,Nr_5hmc.files), function(.ele) 
  toGRanges(.ele, format="narrowPeak"))
names(Nr_5hmc.peaks) <- gsub(".*/(GSM\\d+)_.*", "\\1", names(Nr_5hmc.peaks))
chain <- import.chain("Benign_PC_GSE144530_5hmc_peaks_hg19/hg19ToHg38.over.chain")
Nr_5hmc.peaks.hg38 = lapply(Nr_5hmc.peaks, 
                            function(x) liftOver(x,chain) %>% 
                              as.data.frame() %>% 
                              dplyr::select(-c(group, group_name)) %>%
                              dplyr::mutate(length = end-start) %>% 
                              dplyr::mutate(neglog10q = -log10(qValue)) %>%
                              toGRanges())

peaks_5hmc <- c(peaks_wcdt_cfdna_incl, peaks_localized_incl, 
                peaks_ubc_cfdna_incl, peaks_wcdt_tissue_incl, 
                Nr_5hmc.peaks.hg38)

meta_data <- read_excel("supplementry_tables/can-22-1123_supplementary_table_2_suppst2.xlsx")
peak_name <- data.frame(sample_id = names(peaks_5hmc))
peak_name <- peak_name %>% 
  dplyr::mutate(category = case_when(sample_id %in% names(peaks_wcdt_cfdna_incl) ~ "15-matchedTissue-cfDNA",
                          sample_id %in% names(peaks_localized_incl) ~ "64-localized-PCa",
                          sample_id %in% names(peaks_wcdt_tissue_incl) ~ "103-tissue",
                          sample_id %in% names(peaks_ubc_cfdna_incl) ~ "64-cfDNA-before1stARSI"))
meta_data.flt <- meta_data %>% dplyr::select(sample_id, sample_type, disease_stage)
meta_data.flt <- meta_data.flt %>% dplyr::left_join(peak_name, by = "sample_id")
meta_data.flt[is.na(meta_data.flt)] <- "5-Benign-prostate"
#----- CpG islands
cpg_islans <- read.delim("../scripts/hg38_CpG.bed", header = F)
colnames(cpg_islans) <- c("chr", "start", "end", "id")
###########################################################################################################
###########################################################################################################
# Load WGBS data
###########################################################################################################
###########################################################################################################
# Load data used in the Zhao et al. Nat Gen 2020. Created with /data1/projects/WCDT_WGBS_2019/WCDT_WGBS/scripts/2019_05_15_prepare_environment.r on 2021-12-21
# Data and script available on https://github.com/DavidQuigley/WCDT_WGBS
## WGBS data
#fn_wgbs_environment <- createFn("/data/WGBS_rdata_NO_UPLOAD/2019_05_15_prepare_environment.RData")
#fn_hmrseg_upd <- createFn("/data/WGBS_additional_data/hmrseg_upd.txt") 
hmrseg_upd <- "WCDT_5hmC_additional_data/WGBS_additional_data/hmrseg_upd.txt"
#fn_hmrlevels <- createFn("/data/WGBS_additional_data/rhmr.tsv")
rhmrlevels <- "WCDT_5hmC_additional_data/WGBS_additional_data/rhmr.tsv" #rHMR: recurent HMRs
#fn_promoter_methylation <- createFn("/data/WGBS_additional_data/promoter_methylation_wcdt_tissue.txt")
promoter_methylation <- "WCDT_5hmC_additional_data/WGBS_additional_data/promoter_methylation_wcdt_tissue.txt"

hmrseg_gene <- "~/Desktop/DNA_methylation/Felix_Feng_WGBS_2019_05_15_WGBS/recurrent_HMRs/hmrseg_gene.txt"

tracks <- list()
# load the hmrs-gene data
HMRseggene <- read.delim(hmrseg_gene, stringsAsFactors = F, header = T)
HMRseggene[1:5,]
tracks[["HMRseggene"]] <- makeGRangesFromDataFrame(HMRseggene, keep.extra.columns = T)
# Updated rHMR file (hypomethylated regions called by whole-genome bisulfite sequencing per sample)
hmrseg <- read.table(hmrseg_upd, stringsAsFactors = F, header = T)
head(hmrseg)
tracks[["HMRsegs"]] <- makeGRangesFromDataFrame(hmrseg)
# Methylation levels per sample in rHMRs
rhmrlevel <- read.delim(rhmrlevels, header=T, stringsAsFactors=F,
                        sep='\t', strip.white=T,check.names=F)
rhmrlevel[1:5,1:6]
tracks[["rHMR"]] <- makeGRangesFromDataFrame(rhmrlevel, keep.extra.columns = T)
#head(tracks[["rHMR"]])
# Promoter methylation
promoter_methylation <- read.table(promoter_methylation, header = T, 
                                   stringsAsFactors = F, check.names = F)
promoter_methylation[1:5,1:5]
## generate peaks from rhmrs
sample_name <- colnames(rhmrlevel[-c(1:3)])
rhmr <- list()
for (i in 1:length(sample_name)){
  indx <- i+3
  samplename <- colnames(rhmrlevel[indx])
  df <- rhmrlevel[,c(1,2,3,indx)]
  colnames(df)[4] <- "Beta"
  #df$name <- paste0(samplename, "_", seq(length(rownames(df))))
  gr <- GRanges(df)
  rhmr[[samplename]] <- gr
}

# hmrs
'dir_methylseqR <- "WCDT_5hmC_additional_data/WGBS_additional_data/2019_05_15_WGBS_MethySeekR/"
filenames = list.files(dir_methylseqR, pattern=paste(chr,"UMRsLMRs.txt",sep="_"), 
                       full.names=F)
hmrs = list()
for(curfile in filenames) {
  nametab = strsplit(curfile,'/',fixed=T)[[1]]
  samplename = strsplit(nametab[length(nametab)],"_",fixed=T)[[1]][1]
  curfile = paste(dir_methylseqR,curfile,sep="/")
  curtab = read.delim(file=curfile,sep="\t",header=TRUE,stringsAsFactors=F)
  curtab$sample = toupper(samplename)
  hmrs[[samplename]] = curtab %>% toGRanges()
}'
#hmrdf = rbindlist(hmrs)
#rowshmr = overlap(hmrdf$start,hmrdf$end,start,end) #& hmrdf$sample %in% 
#c(sample_ids_wgbs,samples_normal,samples_upitt_benign_prostate,toupper(samples_upitt_localized_tumor))

##########################################################################
# tpm gene expression
dir_tpm <- "~/Desktop/DNA_methylation/Felix_Feng_WGBS_2019_05_15_WGBS/RNAseq/"
file_tpm = list.files(dir_tpm, pattern="tpm")
stopifnot(length(file_tpm)==1)
filenamefull = paste(dir_tpm,file_tpm,sep='/')
tpm = read.delim(filenamefull,sep='',row.names=1,header=TRUE,
                 stringsAsFactors=FALSE,check.names=F)
tpm_NS_T <- tpm %>% dplyr::select(matches("-NS-T"))
tpm_T <- tpm %>% dplyr::select(!matches("-NS-T"))
colnames(tpm_T) <- gsub("-T-RNA", "", colnames(tpm_T))
colnames(tpm_NS_T) <- gsub("-NS-T-RNA", "", colnames(tpm_NS_T))
 

