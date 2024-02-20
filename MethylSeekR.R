#perfomr MethylSeekR for WGBS-Bismark outut data
## Load the required libraries
#BiocManager::install("BSgenome")
library(BSgenome)
available.genomes()
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
require("BSgenome.Hsapiens.UCSC.hg38")
sLengths <- seqlengths(Hsapiens)
head(sLengths)
#BiocManager::install("MethylSeekR")
library(MethylSeekR)
library(dplyr)
library(rtracklayer)
library(IRanges)
library(parallel)
detectCores()
set.seed(123)

round(memory.limit()/2^20, 2)

#import the bismart output 
setwd('~/Desktop/DNA_methylation/Bismark_output/Bismark_CpG_report_mrgd_output/')
bismart_CpG_report_dir <- "~/Desktop/DNA_methylation/Bismark_output/Bismark_CpG_report_mrgd_output/"
filenames <- list.files(bismart_CpG_report_dir, pattern='merged_CpG_evidence.cov', 
                       full.names=F)
bismark_CpG_list <- lapply(filenames, 
                           function(x) read.delim(x, header = F,
                                            col.names=c("chr", "start", "end", 
                                                        "methylation_prcnt", 
                                                        "count_methylated", 
                                                        "count_unmethylated")))
str(bismark_CpG_list)
names(bismark_CpG_list) <- gsub(".bismark.*", "", filenames)
#----- convert to MethylSeeKR format: 
# chr, pos, T, M
# chromosome, genomic coordinate, total number of reads, number of reads showing methylation
chrs <- c(paste0("chr", seq(1,22)), "chrX", "chrY")
ls_flt <- lapply(bismark_CpG_list, FUN = function(x) x %>% 
                   dplyr::filter(chr %in% chrs) %>%
                   dplyr::mutate(T = count_methylated+count_unmethylated) %>%
                   dplyr::rename(M = count_methylated) %>% 
                   dplyr::select(chr, start, T, M))
str(ls_flt)
# create gragne for each dataset
ls_gr <- list()
for (i in 1:length(ls_flt)){
  name <- names(ls_flt[i])
  gr <- GRanges(seqnames = ls_flt[[i]]$chr,
                ranges = IRanges(start = ls_flt[[i]]$start),
                strand = "*", 
                T = as.numeric(ls_flt[[i]]$T), 
                M = as.numeric(ls_flt[[i]]$M))
  ls_gr[[name]] <- gr
}

# Segmentation of partially methylated domains (PMD)
#The distribution of α-values for each chromosome (PMD quantification)
plot_alpha_dits <- function(gr_list, pdf_file, chroms){
  for (chr in chroms){
    for (i in 1:length(ls_gr)){
      gr_name <- names(ls_gr[i])
      gr_obj <- ls_gr[[i]]
      chr_gr <- subset(gr_obj, seqnames(gr_obj) == chr)
      plotAlphaDistributionOneChr(m=chr_gr,chr.sel=chr,num.cores=9, 
                                  pdfFilename=paste0(pdf_file,gr_name,
                                                       "_alpha_dis_",chr,".pdf"))
    }
  }
}

pdf_file <- "~/Desktop/DNA_methylation/Bismark_output/MethySeekR/"
#chrs <- c("chr11","chr14","chr19","chr22")
plot_alpha_dits(ls_gr, pdf_file, chrs)

train_pmd_on_chroms <- function(gr_list, dir_to_save, chroms, sLengths){
  for (chr in chroms){
    for (i in 1:length(ls_gr)){
      gr_name <- names(ls_gr[i])
      print(gr_name)
      gr_obj <- ls_gr[[i]]
      chr_gr <- subset(gr_obj, seqnames(gr_obj) == chr)
      PMDsegments.gr <- segmentPMDs(m=chr_gr, chr.sel=chr,
                                    seqLengths=sLengths, num.cores=9)
      savePMDSegments(PMDs=PMDsegments.gr,
                      GRangesFilename=paste0(dir_to_save,gr_name,"_MethylSeekR_",
                                             chr,"_PMD.gr.rds"))
    }
  }
}

train_pmd_on_chroms(ls_gr, pdf_file, chrs, sLengths)
plot_alpha_dits_byOverlap <- function(gr_list, pmd_gr, pdf_file, chroms){
  for (chr in chroms){
    for (i in 1:length(ls_gr)){
      gr_name <- names(ls_gr[i])
      gr_obj <- ls_gr[[i]]
      chr_gr <- subset(gr_obj, seqnames(gr_obj) == chr)
      plotAlphaDistributionOneChr(m=subsetByOverlaps(chr_gr,pmd_gr[values(pmd_gr)$type=="notPMD"]), 
                                  chr.sel=chr,
                                  num.cores=9,
                                  pdfFilename = paste0(pdf_file,gr_name,
                                                       "_alpha_dis_noPMD",
                                                       chr,".pdf"))
    }
  }
}

G23_chr22 <- subset(ls_gr$G23, seqnames(ls_gr$G23) == "chr22")
G23_pmd_chr22 <- readRDS("../MethySeekR/G23_MethylSeekR_chr22_PMD.gr.rds")

I24_chr22 <- subset(ls_gr$I24, seqnames(ls_gr$G24) == "chr22")
I24_pmd_chr22 <- readRDS("../MethySeekR/I24_MethylSeekR_chr22_PMD.gr.rds")
plotAlphaDistributionOneChr(m=subsetByOverlaps(I24_chr22,I24_pmd_chr22[values(I24_pmd_chr22)$type=="notPMD"]), 
                            chr.sel="chr22",
                            num.cores=9,
                            pdfFilename = paste0(pdf_file,"I24_alpha_dis_noPMD_chr22.pdf"))

plotPMDSegmentation(m=G23_chr22, segs=G23_pmd_chr22,
                    pdfFilename = paste0(pdf_file,"G23_PMDplot_chr22.pdf"))

# Calculate FDRs with MethylSeekR
#----- CpG islands
#session <- browserSession()
#genome(session) <- "hg38"
#query <- ucscTableQuery(session, table = "cpgIslandExt")
#CpGislands.gr <- track(query)
#genome(CpGislands.gr) <- NA
cpg_gr <- import("../../scripts/hg38_CpG.bed", format = "bed")
# Extend the CpG islands to 5kb to make sure that all unmethylated CpGs lying 
#  in CpG islands are excluded from the FDR calculation.
cpg_gr <- suppressWarnings(resize(cpg_gr, 5000, fix="center"))
# calculate FDRs and run UMRLMR
calc_URLMR <- function(ls_gr,path_to_pmd, chroms, cpg_gr, path2save_plt,FDR.cutoff,m.sel){
  #stat_list <- list()
  sLengths <- seqlengths(Hsapiens)
  for (chr in chroms){
    for (i in 1:length(ls_gr)){
      gr_name <- names(ls_gr[i])
      gr_obj <- ls_gr[[i]]
      chr_gr <- subset(gr_obj, seqnames(gr_obj) == chr)
      chr_pmd <- readRDS(paste0(path_to_pmd,gr_name,"_MethylSeekR_",chr,"_PMD.gr.rds"))
      print(paste0("processing ",gr_name," ",chr))
      cur_stats <- calculateFDRs(m=chr_gr, CGIs=cpg_gr,PMDs=chr_pmd, num.cores=9, 
                                 pdfFilename = paste0(path2save_plt,gr_name,"_FRD_",chr,".pdf"))
      n.sel <- as.integer(names(cur_stats$FDRs[as.character(m.sel), ]
                                [cur_stats$FDRs[as.character(m.sel), ]<FDR.cutoff])[1])
      print(paste0("n.sel=",n.sel))
      #stat_list[[paste0(gr_name,"_",chr)]] <- cur_stats
      UMRLMRsegments.gr <- segmentUMRsLMRs(m=chr_gr, meth.cutoff=m.sel,
                                           nCpG.cutoff=n.sel, PMDs=chr_pmd,
                                           num.cores=9, myGenomeSeq=Hsapiens,
                                           seqLengths=sLengths)
      saveUMRLMRSegments(segs=UMRLMRsegments.gr,
                         TableFilename=paste0(path2save_plt,gr_name,"_methylseekr_",chr,"UMRsLMRs.txt"))
    }
  }
  #return(stat_list)
}

#stats <- calculateFDRs(m=I24_chr22, CGIs=cpg_gr,PMDs=I24_pmd_chr22, num.cores=1)
path_to_pmd <- "../MethySeekR/"
path2save_plt <- "../MethySeekR/"
FDR.cutoff <- 5
m.sel <- 0.5
calc_URLMR(ls_gr,path_to_pmd, chrs, cpg_gr, path2save_plt,FDR.cutoff,m.sel)

