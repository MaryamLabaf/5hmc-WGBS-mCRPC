###########################################################################################################
# Prepare environment for analysis
# WGBS/RNA-seq/5hmc for WCDT and BIDMC mCRPC patients 
###########################################################################################################
# Libraries
libraries_to_load <- c("stringr",
                       "stringi",
                       "dplyr",
                       "plyr",
                       "tidyverse",
                       "factoextra",
                       "readxl",
                       "data.table",
                       "ggplot2",
                       "ggrepel",
                       "RColorBrewer",
                       "viridis",
                       "ggthemes",
                       "ggformula",
                       "reshape2",
                       "ggsci",
                       "ggpubr",
                       "cowplot",
                       "gridExtra",
                       "ggfortify",
                       "gridExtra",
                       "circlize",
                       "rtracklayer",
                       "Sushi",
                       "rCGH",
                       "plotrix",
                       "leaflet",
                       "rstatix",
                       "matrixStats",
                       "ComplexHeatmap",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene",
                       "GenomicFeatures",
                       "abind",
                       "foreach",
                       "scales",
                       "DESeq2",
                       "fgsea",
                       "msigdbr",
                       "MASS",
                       "qdapTools",
                       "survival",
                       "rms",
                       "survminer",
                       "ChIPseeker",
                       "ChIPpeakAnno")

lapply(libraries_to_load,require,character.only = T)

#BASE_DIR = '~/mnt/data1/projects/WCDT_WGBS_2019/WCDT_WGBS'
BASE_DIR = '~/Balk Lab Sequencing/Maryam/WGBS_5hmc_WCDT_mCRPC'
make_path = function(BASE_DIR, s){
  paste(BASE_DIR, s, sep='/' )
}

#######################################################################################################
# Load all custom function calls
#######################################################################################################

source(make_path(BASE_DIR, 'scripts/helper_peaks.R') )

#######################################################################################################
# File locations
#######################################################################################################
dir_gencode.v28.mapping <- make_path(BASE_DIR,"meta_data/gencode.v28.mapping.txt")
dir_hg38_models <- make_path(BASE_DIR,"meta_data/hg38_models.txt")
dir_gencode28_exons <- make_path(BASE_DIR,"meta_data/gencode28_exons.txt")
dir_cpg <- make_path(BASE_DIR,"meta_data/hg38_CpG.bed")
dir_meta <- make_path(BASE_DIR,"meta_data/can-22-1123_supplementary_table_2_suppst2.xlsx")

dir_normal_5hmc <- make_path(BASE_DIR,"5hmc_data/Benign_PC_GSE144530_5hmc_peaks_hg19/") # 5hmc data dir
dir_hg19ToHg38 <- make_path(BASE_DIR,"5hmc_data/Benign_PC_GSE144530_5hmc_peaks_hg19/hg19ToHg38.over.chain")
dir_peaks_wcdt_cfdna <- make_path(BASE_DIR,"5hmc_data/WCDT_hmc_peaks/peaks_wcdt_cfdna_incl.rds")
dir_peaks_localized <- make_path(BASE_DIR,"5hmc_data/WCDT_hmc_peaks/peaks_localized_incl.rds")
dir_peaks_ubc_cfdna <- make_path(BASE_DIR,"5hmc_data/WCDT_hmc_peaks/peaks_ubc_cfdna_incl.rds")
dir_peaks_wcdt_tissue <- make_path(BASE_DIR,"5hmc_data/WCDT_hmc_peaks/peaks_wcdt_tissue_incl.rds")

dir_wcdt_UMRLMR <- make_path(BASE_DIR,"wgbs_data/WCDT_WGBS_MethySeekR/")
dir_bidmc_UMRLMR <- make_path(BASE_DIR, "wgbs_data/UMRLMR_methylseq_bidmc/")
dir_tpm <- make_path(BASE_DIR, "RNAseq/")
dir_rhmrlevels <- make_path(BASE_DIR, "wgbs_data/rhmr.tsv") #rHMR (recurent HMRs)
#######################################################################################################
# Load data
#######################################################################################################
ensembl2sym <- read.delim(dir_gencode.v28.mapping)
gene_model_hg38 <- read.delim(dir_hg38_models)
exons <- read.delim(dir_gencode28_exons)

cpg_islans <- read.delim(dir_cpg, header = F) # cpg islands
colnames(cpg_islans) <- c("chr", "start", "end", "id")

meta_data <- read_excel(dir_meta)
# Normal 5hmc
N5hmc_filenames = list.files(dir_normal_5hmc, pattern='.narrowPeak',full.names=F)
Nr_5hmc.files <- dir("Benign_PC_GSE144530_5hmc_peaks_hg19/.", pattern="*.narrowPeak$")
Nr_5hmc.peaks <- sapply(paste0(dir_normal_5hmc,N5hmc_filenames), function(.ele) toGRanges(.ele, format="narrowPeak"))
names(Nr_5hmc.peaks) <- gsub(".*/(GSM\\d+)_.*", "\\1", names(Nr_5hmc.peaks))
chain <- import.chain(dir_hg19ToHg38)
Nr_5hmc.peaks.hg38 = lapply(Nr_5hmc.peaks, 
                            function(x) liftOver(x,chain) %>% 
                              as.data.frame() %>% 
                              dplyr::select(-c(group, group_name)) %>%
                              dplyr::mutate(length = end-start) %>% 
                              dplyr::mutate(neglog10q = -log10(qValue)) %>%
                              toGRanges())
# mCRPC 5hmc
wcdt_cfdna <- readRDS(dir_peaks_wcdt_cfdna)
peaks_wcdt_cfdna <- lapply(wcdt_cfdna, ChIPQC:::GetGRanges, simple = TRUE)
#head(peaks_wcdt_cfdna[[1]])

localized <- readRDS(dir_peaks_localized)
peaks_localized <- lapply(localized,ChIPQC:::GetGRanges, simple = TRUE)
#head(peaks_localized[[1]])

ubc_cfdna <- readRDS(dir_peaks_ubc_cfdna)
peaks_ubc_cfdna <- lapply(ubc_cfdna,ChIPQC:::GetGRanges, simple = TRUE)
#head(peaks_ubc_cfdna[[1]])

wcdt_tissue <- readRDS(dir_peaks_wcdt_tissue)
peaks_wcdt_tissue <- lapply(wcdt_tissue,ChIPQC:::GetGRanges, simple = TRUE)
#head(peaks_wcdt_tissue[[1]])

# merge all the 5hmc peaks from all the samples 
peaks_5hmc <- c(peaks_wcdt_cfdna, peaks_localized,peaks_ubc_cfdna, 
                peaks_wcdt_tissue,Nr_5hmc.peaks.hg38)

peak_name <- data.frame(sample_id = names(peaks_5hmc))
peak_name <- peak_name %>% 
  dplyr::mutate(category = case_when(sample_id %in% names(peaks_wcdt_cfdna) ~ "15-matchedTissue-cfDNA",
                                     sample_id %in% names(peaks_localized) ~ "64-localized-PCa",
                                     sample_id %in% names(peaks_wcdt_tissue) ~ "103-tissue",
                                     sample_id %in% names(peaks_ubc_cfdna) ~ "64-cfDNA-before1stARSI"))

meta_data.flt <- meta_data %>% 
  dplyr::select(sample_id, sample_type, disease_stage) %>% 
  dplyr::left_join(peak_name, by = "sample_id")
meta_data.flt[is.na(meta_data.flt)] <- "5-Benign-prostate"

###########################################################################################################
# Load WGBS data
###########################################################################################################
# Methylation levels per sample in rHMRs
rhmrlevel <- read.delim(dir_rhmrlevels, header=T, stringsAsFactors=F,sep='\t', strip.white=T,check.names=F)



