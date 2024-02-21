# Sample identifiers
setwd("~/Desktop/DNA_methylation/Felix_Feng_5hmC_CancerRes_2022/")
###########################################################################################################
###########################################################################################################
# Samples
samples_normal <- c("GSM4290243","GSM4290244","GSM4290245","GSM4290246","GSM4290247")

samples_wcdt_tissue_nt <- c("NT-DTB-003-BL","NT-DTB-021-BL","NT-DTB-037-BL","NT-DTB-074-BL","NT-DTB-077-PRO","NT-DTB-188-BL","NT-DTB-194-PRO")

samples_mcrpc_uniq_all_data <- c("DTB-003-BL","DTB-005-BL","DTB-008-BL","DTB-011-BL","DTB-018-BL","DTB-019-PRO","DTB-020-BL","DTB-021-BL","DTB-022-BL","DTB-023-BL","DTB-024-PRO","DTB-034-BL", "DTB-036-BL","DTB-037-BL","DTB-040-BL","DTB-042-BL","DTB-053-BL","DTB-055-PRO","DTB-059-BL","DTB-060-BL","DTB-061-BL","DTB-063-BL","DTB-064-BL","DTB-067-PRO","DTB-069-BL","DTB-071-BL","DTB-074-BL","DTB-077-PRO","DTB-080-BL","DTB-083-BL","DTB-085-BL","DTB-089-BL","DTB-090-PRO","DTB-091-BL","DTB-092-BL","DTB-094-BL","DTB-097-PRO","DTB-098-PRO2","DTB-100-BL","DTB-101-BL","DTB-104-BL","DTB-111-PRO","DTB-119-PRO","DTB-121-BL","DTB-124-BL","DTB-126-BL","DTB-127-PRO","DTB-128-BL","DTB-129-BL","DTB-132-BL","DTB-137-PRO","DTB-138-BL","DTB-140-BL","DTB-141-BL","DTB-143-BL","DTB-146-BL","DTB-149-BL","DTB-151-BL","DTB-156-BL","DTB-159-BL","DTB-165-PRO","DTB-167-PRO","DTB-170-BL","DTB-172-BL","DTB-173-BL","DTB-175-BL","DTB-176-BL","DTB-186-BL","DTB-187-BL","DTB-188-BL","DTB-190-BL","DTB-194-PRO","DTB-201-PRO","DTB-202-BL","DTB-204-BL","DTB-205-BL","DTB-206-BL","DTB-210-BL","DTB-213-BL","DTB-214-BL","DTB-216-PRO","DTB-222-BL","DTB-223-BL","DTB-232-PRO","DTB-234-BL","DTB-251-BL","DTB-252-BL","DTB-255-BL","DTB-258-BL","DTB-260-BL","DTB-261-BL","DTB-265-PRO","DTB-266-BL")

samples_wcdt_cfdna <-  c("DTB-119-PRO_cfdna","DTB-127-PRO_cfdna","DTB-140-BL_cfdna","DTB-149-BL_cfdna","DTB-149-PRO_cfdna","DTB-165-PRO_cfdna","DTB-202-BL_cfdna","DTB-214-BL_cfdna","DTB-216-BL_cfdna","DTB-216-PRO_cfdna","DTB-234-BL_cfdna","DTB-252-BL_cfdna","DTB-258-BL_cfdna","DTB-261-BL_cfdna","DTB-265-BL_cfdna")

samples_wcdt_tissue_mcrpc <- unique(c(samples_mcrpc_uniq_all_data,gsub("_cfdna","",samples_wcdt_cfdna))) 

intersect(samples_mcrpc_uniq_all_data, samples_wcdt_cfdna)
samples_localized_tissue <- c("CPCG0196","CPCG0208","CPCG0210","CPCG0235","CPCG0236","CPCG0238","CPCG0241","CPCG0246","CPCG0248","CPCG0249","CPCG0250","CPCG0256","CPCG0258","CPCG0260" ,"CPCG0263","CPCG0266","CPCG0324","CPCG0331","CPCG0336","CPCG0341","CPCG0342","CPCG0344","CPCG0345","CPCG0346","CPCG0348","CPCG0349","CPCG0350","CPCG0352","CPCG0353","CPCG0354","CPCG0355","CPCG0356","CPCG0357","CPCG0358","CPCG0360","CPCG0361","CPCG0362","CPCG0368","CPCG0369","CPCG0371","CPCG0372","CPCG0373","CPCG0374","CPCG0375","CPCG0377","CPCG0379","CPCG0380","CPCG0382","CPCG0391","CPCG0394","CPCG0458","CPCG0545")

samples_ubc_cfdna <- c("ubc1","ubc10","ubc11","ubc12","ubc13","ubc14","ubc15","ubc17","ubc18","ubc19","ubc2","ubc20","ubc21","ubc22","ubc23","ubc25","ubc26","ubc27","ubc28","ubc29","ubc3","ubc31","ubc32","ubc33","ubc34","ubc35","ubc36","ubc37","ubc38","ubc39","ubc4","ubc40","ubc41","ubc42","ubc43","ubc44","ubc45","ubc46","ubc47","ubc48","ubc5","ubc51","ubc52","ubc53","ubc54","ubc55","ubc56","ubc57","ubc58","ubc59","ubc6","ubc60","ubc61","ubc62","ubc63","ubc64","ubc65","ubc66","ubc67","ubc68","ubc7","ubc70","ubc8","ubc9")

# Get patient identifiers from WGBS
samples_upitt_benign_adjacent_prostate = c('Bis159AT', 'Bis165AT', 'Bis171AT', 'Bis49AT')
samples_upitt_localized_tumor = c('Bis158T', 'Bis159T', 'Bis165T', 'Bis171T', 'Bis49T')
tscnc_pids = c('DTB-003-BL','DTB-036-BL', 'DTB-040-BL', 'DTB-135-PRO', 'DTB-205-BL')

samples.rm <- c("DTB-053-RP", "DTB-265-BL", "DTB-193-BL") # Remove the samples not included in WGBS
samples.rm <- c(samples.rm, c('DTB-083-BL','DTB-126-BL')) # Remove hypermutated samples
#------------------------------------------------------------------------------------------------------
sample_ids_wgs = c("DTB-003-BL", "DTB-005-BL", "DTB-008-BL", "DTB-009-BL", "DTB-011-BL", "DTB-018-BL", "DTB-019-PRO", "DTB-020-BL", "DTB-021-BL", 
                   "DTB-022-BL", "DTB-023-BL", "DTB-024-PRO", "DTB-034-BL", "DTB-035-BL", "DTB-036-BL", "DTB-037-BL", "DTB-040-BL", "DTB-042-BL", 
                   "DTB-053-BL", "DTB-055-PRO", "DTB-059-BL", "DTB-060-BL", "DTB-061-BL", "DTB-063-BL", "DTB-064-BL", "DTB-067-PRO", "DTB-069-BL", 
                   "DTB-071-BL", "DTB-074-BL", "DTB-077-PRO", "DTB-080-BL", "DTB-083-BL", "DTB-085-BL", "DTB-089-BL", "DTB-090-PRO", "DTB-091-BL", 
                   "DTB-092-BL", "DTB-094-BL", "DTB-097-PRO", "DTB-098-PRO2", "DTB-100-BL", "DTB-101-BL", "DTB-102-PRO", "DTB-104-BL", "DTB-111-PRO", 
                   "DTB-112-BL", "DTB-119-PRO", "DTB-121-BL", "DTB-124-BL", "DTB-126-BL", "DTB-127-PRO", "DTB-128-BL", "DTB-129-BL", "DTB-132-BL", 
                   "DTB-135-PRO", "DTB-137-PRO", "DTB-138-BL", "DTB-140-BL", "DTB-141-BL", "DTB-143-BL", "DTB-146-BL", "DTB-149-BL", "DTB-151-BL", 
                   "DTB-156-BL", "DTB-159-BL", "DTB-165-PRO", "DTB-167-PRO", "DTB-170-BL", "DTB-172-BL", "DTB-173-BL", "DTB-175-BL", "DTB-176-BL", 
                   "DTB-183-BL", "DTB-186-BL", "DTB-187-BL", "DTB-188-BL", "DTB-190-BL", "DTB-193-BL", "DTB-194-PRO", "DTB-201-PRO", "DTB-202-BL", 
                   "DTB-204-BL", "DTB-205-BL", "DTB-206-BL", "DTB-210-BL", "DTB-213-BL", "DTB-214-BL", "DTB-216-PRO", "DTB-220-BL", "DTB-222-BL", 
                   "DTB-223-BL", "DTB-232-PRO", "DTB-234-BL", "DTB-251-BL", "DTB-252-BL", "DTB-255-BL", "DTB-258-BL", "DTB-260-BL", "DTB-261-BL", 
                   "DTB-265-PRO", "DTB-266-BL")
sample_ids_wgbs = setdiff(sample_ids_wgs, "DTB-193-BL")

sample_blacklist = c('DTB-004-BL', 'DTB-053-RP', 'DTB-265-BL')
samples_hypermut = c('DTB-083-BL','DTB-126-BL')

samples_sra_benign_prostate = c('SRR6156035','SRR6156036','SRR6156037')
samples_sra_localized_tumor = c('SRR6156038','SRR6156039','SRR6156040',
                                'SRR6156041','SRR6156042','SRR6156043','SRR6156044',
                                'SRR6156045','SRR6156046','SRR6156047','SRR6156048')
#samples_upitt_benign_prostate = toupper(samples_upitt_benign_prostate)
#samples_upitt_localized_tumor = toupper(samples_upitt_localized_tumor)

samples_normal = c('NT-DTB-003-BL','NT-DTB-021-BL','NT-DTB-035-BL','NT-DTB-037-BL',
                   'NT-DTB-074-BL','NT-DTB-077-Pro','NT-DTB-104-BL','NT-DTB-124-BL', 
                   'NT-DTB-188-BL','NT-DTB-194-Pro')
samples_normal_filecase = samples_normal
samples_normal = toupper(samples_normal)
#locs_normal = metastasis_locations[gsub('NT-','',samples_normal)]
#names(locs_normal) = paste('NT-',names(locs_normal),sep='')

samples_upitt_benign_prostate = c('Bis49AT','Bis159AT','Bis165AT','Bis171AT')
samples_upitt_localized_tumor <- c('Bis49T','Bis158T','Bis159T', 'Bis165T', 'Bis171T')
samples_upitt_benign_prostate <- toupper(samples_upitt_benign_prostate)
#samples_upitt_localized_tumor <- toupper(samples_upitt_localized_tumor)

samples_braf = "DTB-022-BL"
samples_idh1 = "DTB-112-BL"
samples_tet2 <- c("DTB-023-BL", "DTB-188-BL", "DTB-202-BL", "DTB-252-BL") #252 is the one out of CMP
samples_dnmt3b <- c("DTB-024-PRO", "DTB-094-BL")


samples_tscnc = c('DTB-003-BL','DTB-036-BL','DTB-040-BL','DTB-135-PRO','DTB-205-BL')

samples_cmp = c(
  "DTB-018-BL",
  "DTB-021-BL",
  "DTB-022-BL",  #BRAF
  "DTB-023-BL",  #TET2
  "DTB-024-PRO", #DNMT3B
  "DTB-037-BL",
  "DTB-042-BL",
  "DTB-074-BL",
  "DTB-089-BL",
  "DTB-094-BL",  #DNMT3B
  "DTB-102-PRO",
  "DTB-112-BL",  #IDH1
  "DTB-124-BL",
  "DTB-137-PRO",
  "DTB-141-BL",
  "DTB-151-BL",
  "DTB-188-BL",  #TET2
  "DTB-190-BL",
  "DTB-194-PRO",
  "DTB-202-BL",# TET2
  "DTB-204-BL",
  "DTB-260-BL")

WGBS_sup <- read_excel("WCDT_5hmC_additional_data/WGBS_additional_data/WGBS_supplementry.xlsx",
                       sheet = 1, skip = 1)