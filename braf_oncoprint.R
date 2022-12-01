## Pinky Langat 2020
## This code is built from nmf_oncoprint by Eleanor Woodward
## cBio Portal example (https://jokergoo.github.io/ComplexHeatmap-reference/book/oncoprint.html#apply-to-cbioportal-dataset)

#install.packages("viridis")
#install.packages("alluvial")
#install.packages("ggalluvial")
#install.packages("hrbrthemes")
#install.packages("tidytidbits")
#install.packages("survivalAnalysis")
#install.packages("segmented")

library(ComplexHeatmap)
library(lubridate)
library(viridis)

library(alluvial)
library(ggalluvial)
library(ggsci)

library(survminer)
library(survival)

library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
library(RColorBrewer)
library(segmented)

## Load in required data
source("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/2-Code/MafFunctions_BRAF_PL.R")

load("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/braf_oncoprint_data.RData")

braf_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/alterations_braf_samples.txt", stringsAsFactors = F)

jhh_braf_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/JHH_alterations_braf_samples.txt", stringsAsFactors = F)

braf_muts_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/braf_mutations.txt", stringsAsFactors = F)

jhh_braf_muts_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/JHH_braf_mutations.txt", stringsAsFactors = F)

braf_arm_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/braf_arm-level_calls_chr.txt", stringsAsFactors = F)

## *** check these fusion data files
fusion_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/braf_fusion_samples.txt", stringsAsFactors = F)

dfci_fusion_data <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/braf_dfci_rearrangements_sv.csv", stringsAsFactors = F) # compiled from REQ_WLB13_91353_GENOMIC_SPECIMEN.csv

genie_fusion_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/GENIE/data_fusions_5.0-public.txt", stringsAsFactors = F)

master.sheet.dfci <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/DFCI/master_sheet_dfci_updated.txt", stringsAsFactors = F)

path_subtypes_who <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/classifications/path_subtypes_who.2021-07-09.txt", stringsAsFactors = F)

dfci.req.genomic.specimen  <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/Oncopanel/REQ_WLB13_91353_GENOMIC_SPECIMEN.csv", stringsAsFactors = F)

dfci.req.cancer.diagnosis.careg <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/Oncopanel/REQ_WLB13_91353_CANCER_DIAGNOSIS_CAREG.csv", stringsAsFactors = F)

dfci.req.treatment.plan <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/Oncopanel/REQ_WLB13_91353_TREATMENT_PLAN.csv", stringsAsFactors = F)

dfci.req.demographics.reg <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/Oncopanel/REQ_WLB13_91353_DEMOGRAPHICS_REGISTRATION.csv", stringsAsFactors = F)

dfci.req.genomic.mutations <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/Oncopanel/REQ_WLB13_91353_GENOMIC_MUTATION_RESULTS.csv", stringsAsFactors = F)

dfci.req.genomic.sv <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/Oncopanel/REQ_WLB13_91353_GENOMIC_SV_RESULTS.csv", stringsAsFactors = F)

adult.outcomes.verified <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/outcomes/braf_adult_outcomes_verified.txt", stringsAsFactors = F)

peds.outcomes.verified <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/peds/braf_peds_outcomes_verified.2021-04-09.txt", stringsAsFactors = F)

tcga.lgg.outcomes <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/LGG_Outcomes_full.txt", stringsAsFactors = F)

tcga.gbm.outcomes <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/GBM_Outcomes.txt", stringsAsFactors = F)

tcga.lgg.mutations <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/LGG_Mutations.txt", stringsAsFactors = F, comment.char = "#")

tcga.lgg.cnv <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/LGG_focal_CNV.txt", stringsAsFactors = F, comment.char = "#")

tcga.lgg.cnv.broad <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/LGG_broad_CNV.txt", stringsAsFactors = F, comment.char = "#")

tcga.lgg.svs <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/LGG_WGS_SVs.txt", stringsAsFactors = F, comment.char = "#")

tcga.lgg.svs.exome <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/LGG_Exome_SVs.txt", stringsAsFactors = F, comment.char = "#")

tcga.lgg.svs.low <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/LGG_LowPass_SVs.txt", stringsAsFactors = F, comment.char = "#")

tcga.lgg.indels <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/LGG_Indels.txt", stringsAsFactors = F, comment.char = "#")

tcga.gbm.mutations <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/GBM_Mutations.txt", stringsAsFactors = F, comment.char = "#")

tcga.gbm.cnv <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/GBM_focal_CNV.txt", stringsAsFactors = F, comment.char = "#")

tcga.gbm.cnv.broad <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/GBM_broad_CNV.txt", stringsAsFactors = F, comment.char = "#")

tcga.gbm.svs <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/TCGA/GBM_SVs.txt", stringsAsFactors = F, comment.char = "#")

tcga.fromJR <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/rna-seq/tcga_braf.fromJR.txt", stringsAsFactors = F, comment.char = "#")

genie.outcomes <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/GENIE/combined_study_clinical_data.tsv", stringsAsFactors = F)

genie.demographics <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/GENIE/data_clinical_patient_5.0-public.txt", stringsAsFactors = F)

genie.mutations <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/GENIE/data_mutations_extended_5.0-public_edit.txt", stringsAsFactors = F)
#genie.mutations <- NULL
#genie.CNA <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/raw/GENIE/data_CNA_5.0-public.txt", stringsAsFactors = F, comment.char = "#")

ras_path_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/ras-pathway/ras-pathway-gene-names.txt", stringsAsFactors = F)

jhh_genes_v2 <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/JHH/SolidTumorPanel-II_GeneList_v2.0_JHH.txt", stringsAsFactors = F)
jhh_genes_v3 <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/JHH/SolidTumorPanel-II_GeneList_v3.0_JHH.txt", stringsAsFactors = F)
jhh_genes_v3.2 <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/JHH/SolidTumorPanel-II_GeneList_v3.2_JHH.txt", stringsAsFactors = F)
jhh_genes_v4 <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/JHH/SolidTumorPanel-II_GeneList_v4.0_JHH.txt", stringsAsFactors = F)
jhh_genes_npp <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/JHH/NPP_GeneList_JHH.txt", stringsAsFactors = F)
jhh_genes_stpl <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/JHH/STPL_GeneList_JHH.txt", stringsAsFactors = F)

#os_data <- read.delim("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/clinical-master-sheets/braf_all_overall_survival.2021-03-17.txt", stringsAsFactors = F)

jhh_os_data <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/JHH/JHH_DeIDdata_formatted_revised.2021-07-09.csv", stringsAsFactors = F)

jhh_panel_data <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/institutional/JHH/jhh_gene_panel_info.csv", stringsAsFactors = F)

censored.adult.braf.mastersheet <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/clinical-master-sheets/censored_braf_adult_glioma_mastersheet.csv", stringsAsFactors = F)

censored.peds.braf.mastersheet <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/clinical-master-sheets/censored_braf_peds_glioma_mastersheet.2021-07-09.csv", stringsAsFactors = F)

 #### ONCOPRINT START #### 
## Quick revising input data files
# Quick add of TCGA mutation column for sample accessions
tcga.gbm.mutations$SAMPLE_ACCESSION_NBR <- substr(tcga.gbm.mutations$Tumor_Sample_Barcode, 1, 12)
unique(tcga.gbm.mutations$SAMPLE_ACCESSION_NBR)

tcga.lgg.mutations$SAMPLE_ACCESSION_NBR <-substr(tcga.lgg.mutations$Tumor_Sample_Barcode, 1, 12)
unique(tcga.lgg.mutations$SAMPLE_ACCESSION_NBR)

# Quick add of GENIE os for sample accessions
genie.outcomes$SAMPLE_ACCESSION_NBR <- paste("GENIE-MSK-", genie.outcomes$Sample.ID, sep = '')

# Quick subset of genie mutations df ince it is so large
genie.mutations <- genie.mutations[genie.mutations$Tumor_Sample_Barcode %in% sub_samps,]

# Quick revise some columns for JHH OS data
jhh_os_data$SAMPLE_ACCESSION_NBR <- paste("JHH-", jhh_os_data$redcap_patient_id, sep = "")
jhh_os_data$PATIENT_ID <- jhh_os_data$redcap_patient_id

jhh_os_data$Age_Interval[jhh_os_data$age_at_dx <18] <- "Pediatric"
jhh_os_data$Age_Interval[jhh_os_data$age_at_dx >= 18 & jhh_os_data$age_at_dx < 35] <- "Young Adult"
jhh_os_data$Age_Interval[jhh_os_data$age_at_dx > 34] <- "Middle"
jhh_os_data$Age_Interval[jhh_os_data$age_at_dx > 50] <- "Old"
jhh_os_data$age_int2 <- jhh_os_data$Age_Interval

jhh_os_data$age_int2[jhh_os_data$age_int2 == 'Pediatric'] <- "<18"
jhh_os_data$age_int2[jhh_os_data$age_int2 == 'Young Adult'] <- "18-34"
jhh_os_data$age_int2[jhh_os_data$age_int2 == 'Middle'] <- "35-50"
jhh_os_data$age_int2[jhh_os_data$age_int2 == 'Old'] <- ">50"

## Combining DFCI and JHH genomic data
jhh_braf_data
jhh_braf_data[c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")]
all_braf_data <- rbind(braf_data, jhh_braf_data[c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "variant_classification")])
all_braf_data

## This is just getting the data into the raw input format from the reference
all_samples <- unique(all_braf_data$SAMPLE_ACCESSION_NBR)
all_genes <- names(sort(table(all_braf_data$BEST_EFF_GENE), decreasing = T))
all_samples
all_genes

jhh_gen_samples <- unique(all_braf_data$SAMPLE_ACCESSION_NBR[grepl("JHH", all_braf_data$SAMPLE_ACCESSION_NBR)])
jhh_gen_samples

nmf_alts <- all_braf_data[all_braf_data$SAMPLE_ACCESSION_NBR %in% all_samples & all_braf_data$BEST_EFF_GENE %in% all_genes,]
unique(nmf_alts$SAMPLE_ACCESSION_NBR)
unique(nmf_alts$BEST_EFF_GENE)

## Combine some genes that make sense together
nmf_alts$BEST_EFF_GENE[nmf_alts$BEST_EFF_GENE %in% c("CDK4", "CDK6")] <- "CDK4/6"
nmf_alts$BEST_EFF_GENE[nmf_alts$BEST_EFF_GENE %in% c("CDK4", "CDK6")] <- "CDK4/6"
nmf_alts$BEST_EFF_GENE[nmf_alts$BEST_EFF_GENE %in% c("MDM2", "MDM4")] <- "MDM2/4"
nmf_alts$BEST_EFF_GENE[nmf_alts$BEST_EFF_GENE %in% c("PIK3CA", "PIK3R1")] <- "PIK3CA/3R1"
nmf_alts$BEST_EFF_GENE[nmf_alts$BEST_EFF_GENE %in% c("IDH1", "IDH2")] <- "IDH1/2"
nmf_alts$BEST_EFF_GENE[nmf_alts$BEST_EFF_GENE %in% c("CDKN2A", "CDKN2B")] <- "CDKN2A/B"

## Load arm-level CNVs
braf_arm_data <- data.frame(t(braf_arm_data))
rownames(braf_arm_data) <- gsub("\\.", "-", rownames(braf_arm_data))
braf_arm_data$SAMPLE_ACCESSION_NBR <- rownames(braf_arm_data)

## Make dataframe of arm-level CNVs and add to nmf_alts
chr7p <- data.frame(SAMPLE_ACCESSION_NBR = braf_arm_data$SAMPLE_ACCESSION_NBR, BEST_EFF_GENE = "chr7p", variant_classification = braf_arm_data$chr7p)
chr7q <- data.frame(SAMPLE_ACCESSION_NBR = braf_arm_data$SAMPLE_ACCESSION_NBR, BEST_EFF_GENE = "chr7q", variant_classification = braf_arm_data$chr7q)
chr10p <- data.frame(SAMPLE_ACCESSION_NBR = braf_arm_data$SAMPLE_ACCESSION_NBR, BEST_EFF_GENE = "chr10p", variant_classification = braf_arm_data$chr10p)
chr10q <- data.frame(SAMPLE_ACCESSION_NBR = braf_arm_data$SAMPLE_ACCESSION_NBR, BEST_EFF_GENE = "chr10q", variant_classification = braf_arm_data$chr10q)
chr1p <- data.frame(SAMPLE_ACCESSION_NBR = braf_arm_data$SAMPLE_ACCESSION_NBR, BEST_EFF_GENE = "chr1p", variant_classification = braf_arm_data$chr1p)
chr19q <- data.frame(SAMPLE_ACCESSION_NBR = braf_arm_data$SAMPLE_ACCESSION_NBR, BEST_EFF_GENE = "chr19q", variant_classification = braf_arm_data$chr19q)

braf_arm_calls <- rbind(chr7p, chr7q, chr10p, chr10q, chr1p, chr19q)
braf_arm_calls <- braf_arm_calls[!(braf_arm_calls$variant_classification == "0"),]

braf_arm_calls$variant_classification[braf_arm_calls$variant_classification == "-1"] <- "arm_level_loss"
braf_arm_calls$variant_classification[braf_arm_calls$variant_classification == "1"] <- "arm_level_gain"
braf_arm_calls
nmf_alts <- rbind(nmf_alts, braf_arm_calls)

## Reclassify the variant types existing in the data into a few simpler categories
nmf_alts$variant_classification[nmf_alts$variant_classification == "-2"] <- "2DEL"
nmf_alts$variant_classification[nmf_alts$variant_classification == "2"] <- "HA"
nmf_alts$variant_classification[nmf_alts$variant_classification == "2DEL"] <- "loss"
nmf_alts$variant_classification[nmf_alts$variant_classification == "HA"] <- "gain"
table(nmf_alts$variant_classification)

nmf_alts$variant_classification_simple <- "other"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "loss"] <- "focal loss"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "gain"] <- "focal gain"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "arm_level_loss"] <- "arm level loss"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "arm_level_gain"] <- "arm level gain"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "rearrangement"] <- "rearrangement"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "frameshift_indel"] <- "damaging mutation"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "nonsense"] <- "damaging mutation"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "stop_codon"] <- "damaging mutation"
nmf_alts$variant_classification_simple[nmf_alts$variant_classification == "indel"] <- "damaging mutation"
unique(nmf_alts[,c("variant_classification", "variant_classification_simple")])
table(nmf_alts$variant_classification)
table(nmf_alts$variant_classification_simple)

## Not really necessary in current code, meant to ignore alts of any not included samples
#quiet.samples <- all_samples[!(all_samples %in% nmf_alts$SAMPLE_ACCESSION_NBR)]
#nmf_alts <- rbind(nmf_alts, data.frame(SAMPLE_ACCESSION_NBR = quiet.samples, BEST_EFF_GENE = "quiet.mut", variant_classification = "altered", variant_classification_simple = "altered"))

## Make df of just BRAF-related alts
nmf_alts_braf_only <- nmf_alts[nmf_alts$BEST_EFF_GENE == "BRAF",]
nmf_alts_braf_only

## Clarify BRAF fusions
braf_rearrangement_samps <- nmf_alts_braf_only$SAMPLE_ACCESSION_NBR[nmf_alts_braf_only$variant_classification == "rearrangement"]
braf_rearrangement_samps
nmf_alts[nmf_alts$SAMPLE_ACCESSION_NBR %in% braf_rearrangement_samps & nmf_alts$variant_classification=="rearrangement",]

# DFCI rearrangements
dfci_braf_genomic_specimen <- dfci.req.genomic.specimen[dfci.req.genomic.specimen$SAMPLE_ACCESSION_NBR %in% classifier$SAMPLE_ACCESSION_NBR,]
temp_df <- master.sheet.dfci[grepl("fusion", master.sheet.dfci$REPORT_COMMENT.x),]
temp_df2 <- master.sheet.dfci[grepl("fusion", master.sheet.dfci$REPORT_COMMENT.y),]
table(classifier$braf_class[classifier$SAMPLE_ACCESSION_NBR %in% temp_df$SAMPLE_ACCESSION_NBR])

temp_df$SAMPLE_ACCESSION_NBR[temp_df$SAMPLE_ACCESSION_NBR %in% peds.braf.mastersheet$SAMPLE_ACCESSION_NBR]
peds.braf.mastersheet$SAMPLE_ACCESSION_NBR[peds.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% fusion_list]
temp_df$SAMPLE_ACCESSION_NBR[temp_df$SAMPLE_ACCESSION_NBR %in% adult.braf.mastersheet$SAMPLE_ACCESSION_NBR]
adult.braf.mastersheet$SAMPLE_ACCESSION_NBR[adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% fusion_list]
fusion_list[grepl("BL", fusion_list)]

# Most clear approach to DFCI fusion checks
temp_df <- dfci.req.genomic.specimen[dfci.req.genomic.specimen$SAMPLE_ACCESSION_NBR %in% all_samples,]
temp_df <- temp_df[,c("SAMPLE_ACCESSION_NBR", "SV_COMMENT")]
temp_df$BRAF_Class <- classifier$braf_class[match(temp_df$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
temp_df$BRAF_Mutation <- classifier$braf_mutation[match(temp_df$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
temp_df$Subtype <- classifier$subtype[match(temp_df$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
master.sheet.dfci$REPORT_COMMENT.x[master.sheet.dfci$SAMPLE_ACCESSION_NBR == "BL-16-K20293"]

# GENIE
genie_fusions_data$Fusion[genie_fusions_data$Tumor_Sample_Barcode %in% braf_rearrangement_samps & grepl("BRAF", genie_fusions_data$Fusion)]
temp_df <- genie_fusions_data[genie_fusions_data$Tumor_Sample_Barcode %in% (classifier$SAMPLE_ACCESSION_NBR[classifier$braf_mut_spec == "rearrangement" & classifier$braf_class == "Other"]),]
temp_df <- genie_fusions_data[genie_fusions_data$Tumor_Sample_Barcode %in% (classifier$SAMPLE_ACCESSION_NBR[classifier$subtype == "GBM"]) & grepl("BRAF", genie_fusions_data$Fusion),]
temp_df <- genie_fusions_data[genie_fusions_data$Tumor_Sample_Barcode %in% classifier$SAMPLE_ACCESSION_NBR,]

table(classifier$braf_class[classifier$SAMPLE_ACCESSION_NBR %in% temp_df$Tumor_Sample_Barcode])

table(classifier$subtype[classifier$SAMPLE_ACCESSION_NBR %in% genie_fusions_data$Tumor_Sample_Barcode])

## Load known BRAF fusions 
fusion_list <- fusion_data[,1]

fusion_list

## Finding the missing fusions
#nmf_alts$SAMPLE_ACCESSION_NBR[nmf_alts$BEST_EFF_GENE=='PTPRZ1' & nmf_alts$variant_classification=='rearrangement']
#'BL-17-N11461' %in% fusion_list
#nmf_alts$SAMPLE_ACCESSION_NBR[nmf_alts$BEST_EFF_GENE=='BCAS1' & nmf_alts$variant_classification=='rearrangement']
#'BL-17-X49573' %in% fusion_list
fusion_list_addons <- unique(braf_rearrangement_samps[!(braf_rearrangement_samps %in% fusion_list)]) # which rearranged samps are not in fusion_list
fusion_list <- c(fusion_list, fusion_list_addons)
fusion_list

## Classify ALL BRAF alteration categories
nmf_alts_braf_only$BRAF_CLASS <- "Other"
nmf_alts_braf_only$BRAF_CLASS[nmf_alts_braf_only$variant_classification == "gain"] <- "Gain"
#nmf_alts_braf_only$BRAF_CLASS[nmf_alts_braf_only$variant_classification == "rearrangement"] <- "Rearrangement"
nmf_alts_braf_only$BRAF_CLASS[nmf_alts_braf_only$SAMPLE_ACCESSION_NBR %in% fusion_list] <- "Fusion"
table(nmf_alts_braf_only$BRAF_CLASS)
unique(nmf_alts_braf_only[,c("variant_classification", "BRAF_CLASS")])

# Checking for double classified
unique(nmf_alts_braf_only$SAMPLE_ACCESSION_NBR)
table(nmf_alts_braf_only$BRAF_CLASS)
temp_df <- nmf_alts_braf_only[nmf_alts_braf_only$SAMPLE_ACCESSION_NBR %in% nmf_alts_braf_only$SAMPLE_ACCESSION_NBR[duplicated(nmf_alts_braf_only$SAMPLE_ACCESSION_NBR)],]
table(temp_df$SAMPLE_ACCESSION_NBR)

# SNVs
braf_muts_data <- rbind(braf_muts_data, jhh_braf_muts_data)
braf_muts_data
braf_muts_data$BRAF_CLASS <- "Other"
# Class I
braf_muts_data$BRAF_CLASS[grepl("p.V600[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class I"
# Class II
braf_muts_data$BRAF_CLASS[grepl("p.P367[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
braf_muts_data$BRAF_CLASS[grepl("p.G464[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
braf_muts_data$BRAF_CLASS[grepl("p.G469[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
braf_muts_data$BRAF_CLASS[grepl("p.L485.",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
braf_muts_data$BRAF_CLASS[grepl("p.N486.",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
braf_muts_data$BRAF_CLASS[grepl("p.E586[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
braf_muts_data$BRAF_CLASS[grepl("p.L597[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
braf_muts_data$BRAF_CLASS[grepl("p.T599.",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
braf_muts_data$BRAF_CLASS[grepl("p.K601.",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class II"
# Class III
braf_muts_data$BRAF_CLASS[grepl("p.D287[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"
braf_muts_data$BRAF_CLASS[grepl("p.V459[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"
braf_muts_data$BRAF_CLASS[grepl("p.G466[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"
braf_muts_data$BRAF_CLASS[grepl("p.S467[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"
braf_muts_data$BRAF_CLASS[grepl("p.G469E",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"
braf_muts_data$BRAF_CLASS[grepl("p.N581[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"
braf_muts_data$BRAF_CLASS[grepl("p.D594[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"
braf_muts_data$BRAF_CLASS[grepl("p.F595[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"
braf_muts_data$BRAF_CLASS[grepl("p.G596[A-Z]",braf_muts_data$BEST_EFF_PROTEIN_CHANGE)] <- "Class III"

## Just taking stock of the BRAF classifications
unique(braf_muts_data[,c("BEST_EFF_PROTEIN_CHANGE","VARIANT_TYPE", "BRAF_CLASS")])
table(braf_muts_data$BRAF_CLASS)

## Multi mutated BRAF samples
braf_muts_data$SAMPLE_ACCESSION_NBR[duplicated(braf_muts_data$SAMPLE_ACCESSION_NBR)]
braf_muts_data[braf_muts_data$SAMPLE_ACCESSION_NBR %in% c("TCGA-DU-6392", "BL-16-J49955", "BL-14-W13819", "JHH-18"),]
nmf_alts_braf_only$BRAF_CLASS[nmf_alts_braf_only$SAMPLE_ACCESSION_NBR %in% c("TCGA-DU-6392", "BL-16-J49955", "BL-14-W13819", "JHH-18")]

## Compile a basic phenotypes data frame 'classifier'
classifier <- master.sheet.glioma.unique[
  master.sheet.glioma.unique$SAMPLE_ACCESSION_NBR %in% all_samples,]

## Correct for the two samples missing age info
classifier$SAMPLE_ACCESSION_NBR[is.na(classifier$Age)]
classifier$Age[classifier$SAMPLE_ACCESSION_NBR == "BL-15-G14321"] <-
  master.sheet.dfci$Age[master.sheet.dfci$SAMPLE_ACCESSION_NBR == "BL-15-G14321"]
classifier$age_num[classifier$SAMPLE_ACCESSION_NBR == "BL-15-G14321"] <-
  as.numeric(classifier$Age[classifier$SAMPLE_ACCESSION_NBR == "BL-15-G14321"])

classifier$Age[classifier$SAMPLE_ACCESSION_NBR == "BL-16-W42664"] <-
  master.sheet.dfci$Age[master.sheet.dfci$SAMPLE_ACCESSION_NBR == "BL-16-W42664"]
classifier$age_num[classifier$SAMPLE_ACCESSION_NBR == "BL-16-W42664"] <-
  as.numeric(classifier$Age[classifier$SAMPLE_ACCESSION_NBR == "BL-16-W42664"])

## Initialize JHH samples for classifier
jhh.classifier <- data.frame(matrix(ncol = 27, nrow = length(jhh_gen_samples)))
colnames(jhh.classifier) <- colnames(classifier)
jhh.classifier$SAMPLE_ACCESSION_NBR <- jhh_gen_samples
colnames(jhh.classifier)
colnames(jhh_os_data)

## Fill in JHH classifier data
jhh.classifier$PATIENT_ID <- jhh_os_data$PATIENT_ID[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$PANEL_VERSION <- jhh_panel_data$PANEL_VERSION[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_panel_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$cohort <- "JHH"
jhh.classifier$Age <- jhh_os_data$age_at_dx[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$Age_Interval <- jhh_os_data$Age_Interval[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$Grade <- jhh_os_data$who_grade[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$Primary <- jhh_os_data$Primary[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$age_num <- as.integer(jhh.classifier$Age)
jhh.classifier$Gender <- jhh_os_data$patient_gender[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$WHO_Diagnosis <- jhh_os_data$who_dx[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$IDH1.2 <- jhh_os_data$IDH[match(jhh.classifier$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]
jhh.classifier$IDH1.2
jhh.classifier$IDH1.2[jhh.classifier$SAMPLE_ACCESSION_NBR %in% IDH_muts] <- 1
jhh.classifier$IDH1.2

jhh_os_data$SAMPLE_ACCESSION_NBR[jhh_os_data$IDH == 1]

## Combine jhh with dfci classifier
classifier <- rbind(classifier, jhh.classifier)

## Sample names as 'membership'
classifier$membership <- classifier$SAMPLE_ACCESSION_NBR

## Reclassify age_nums
classifier$age_num[classifier$Age == "<18"] <- 17.9

## Reclassify age_ints
table(classifier$Age_Interval)

classifier$Age_Interval[classifier$age_num >= 18 & classifier$age_num < 35] <- "Young Adult"
classifier$Age_Interval[classifier$age_num > 34] <- "Middle"
classifier$Age_Interval[classifier$age_num > 50] <- "Old"
classifier$age_int <- classifier$Age_Interval
classifier$Pediatric[classifier$age_num >= 18.0] <- "Adult"

table(classifier$age_num[classifier$age_int == "Pediatric"])
table(classifier$age_num[classifier$age_int == "Young Adult"])
table(classifier$age_num[classifier$age_int == "Middle"])
table(classifier$age_num[classifier$age_int == "Old"])
table(classifier$age_int)

sort(unique(classifier$Age[classifier$Age_Interval == "Pediatric"]))
sort(unique(classifier$Age[classifier$Age_Interval == "Young Adult"]))
table(classifier$Age_Interval)

## Rename age intervals
classifier$age_int[classifier$age_int == 'Pediatric'] <- "<18"
classifier$age_int[classifier$age_int == 'Young Adult'] <- "18-34"
classifier$age_int[classifier$age_int == 'Middle'] <- "35-50"
classifier$age_int[classifier$age_int == 'Old'] <- ">50"
#classifier$age_int[classifier$age_int == '>50'] <- "51+"
classifier$age_int <- factor(classifier$age_int, levels = c("<18", "18-34", "35-50", ">50"))
table(classifier$age_int)

## Classify broad glioma pathological subtypes
classifier$subtype_path <- classifier$WHO_Diagnosis
classifier$subtype_path <- path_subtypes_who$Subtype[match(classifier$subtype_path, path_subtypes_who$WHO_Diagnosis)]

table(classifier$subtype_path)
path_subtypes_who

## IDH Mutations
# Check IDH mutants position if it is the classical IDH mutations

# Collect all samples with IDH1 or IDH2 mutation
IDH_muts <- unique(nmf_alts$SAMPLE_ACCESSION_NBR[grepl("IDH", nmf_alts$BEST_EFF_GENE) & !(grepl("IDH3", nmf_alts$BEST_EFF_GENE))])
IDH_muts

# Collect IDH mutation info and combine into a df
temp_df <- (dfci.req.genomic.mutations[c("SAMPLE_ACCESSION_NBR", "CANONICAL_GENE", "CANONICAL_PROTEIN_CHANGE")])[dfci.req.genomic.mutations$SAMPLE_ACCESSION_NBR %in% sub_samps & grepl("IDH", dfci.req.genomic.mutations$CANONICAL_GENE),]
names(temp_df) <- c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "AA_CHANGE")
IDH_muts_df <- temp_df

temp_df <- (genie.mutations[c("Tumor_Sample_Barcode", "Hugo_Symbol", "HGVSp_Short")])[genie.mutations$Tumor_Sample_Barcode %in% sub_samps & grepl("IDH", genie.mutations$Hugo_Symbol),]
names(temp_df) <- c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "AA_CHANGE")
IDH_muts_df <- rbind(IDH_muts_df, temp_df)

temp_df <- (tcga.gbm.mutations[c("SAMPLE_ACCESSION_NBR", "Hugo_Symbol", "HGVSp_Short")])[tcga.gbm.mutations$SAMPLE_ACCESSION_NBR %in% IDH_muts & grepl("IDH", tcga.gbm.mutations$Hugo_Symbol),]
names(temp_df) <- c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "AA_CHANGE")
IDH_muts_df <- rbind(IDH_muts_df, temp_df)

temp_df <- (tcga.lgg.mutations[c("SAMPLE_ACCESSION_NBR", "Hugo_Symbol", "HGVSp_Short")])[tcga.lgg.mutations$SAMPLE_ACCESSION_NBR %in% sub_samps & grepl("IDH", tcga.lgg.mutations$Hugo_Symbol),]
names(temp_df) <- c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "AA_CHANGE")
IDH_muts_df <- rbind(IDH_muts_df, temp_df)

temp_df <- (jhh_braf_data[c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "AA_change")])[jhh_braf_data$SAMPLE_ACCESSION_NBR %in% all_samples & grepl("IDH", jhh_braf_data$BEST_EFF_GENE),]
names(temp_df) <- c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE", "AA_CHANGE")
IDH_muts_df <- rbind(IDH_muts_df, temp_df)

# Cleanup IDH mutations df
IDH_muts_df <- IDH_muts_df[!grepl("IDH3", IDH_muts_df$BEST_EFF_GENE),]
IDH_muts_df <- IDH_muts_df[!(IDH_muts_df$AA_CHANGE == ''),]

IDH_muts_df
unique(IDH_muts_df$SAMPLE_ACCESSION_NBR)
IDH_muts_df$SAMPLE_ACCESSION_NBR[!(IDH_muts_df$SAMPLE_ACCESSION_NBR %in% IDH_muts)]
IDH_muts[!(IDH_muts %in% IDH_muts_df$SAMPLE_ACCESSION_NBR)]

# Double check IDH1.2 TRUE/FALSE status
table(classifier$IDH1.2[classifier$SAMPLE_ACCESSION_NBR %in% IDH_muts])
table(classifier$IDH1.2[!(classifier$SAMPLE_ACCESSION_NBR %in% IDH_muts)])
table(classifier$IDH1.2)

# Classify IDH canonical vs other mutation
IDH_muts_df$MUT_TYPE <- "other"
IDH_muts_df$MUT_TYPE[grepl("132",IDH_muts_df$AA_CHANGE)] <- "canonical"
IDH_muts_df$MUT_TYPE[grepl("172",IDH_muts_df$AA_CHANGE)] <- "canonical"

IDH_muts
IDH_muts_canonical <- unique(IDH_muts_df$SAMPLE_ACCESSION_NBR[IDH_muts_df$MUT_TYPE == 'canonical'])
IDH_muts_canonical
IDH_muts_other <- unique(IDH_muts_df$SAMPLE_ACCESSION_NBR[IDH_muts_df$MUT_TYPE == 'other'])
IDH_muts_other <- unique(IDH_muts_other[!(IDH_muts_other %in% IDH_muts_canonical)])
IDH_muts_other

classifier$IDH_mutation <- "Wildtype"
classifier$IDH_mutation[classifier$SAMPLE_ACCESSION_NBR %in% IDH_muts_other] <- "Other"
classifier$IDH_mutation[classifier$SAMPLE_ACCESSION_NBR %in% IDH_muts_canonical] <- "Canonical"

unique(classifier[c("IDH1.2", "IDH_mutation")])

## Checking for cIMPACT GBM classifications
# 1p19q codeletion
chr1p_loss <- nmf_alts$SAMPLE_ACCESSION_NBR[grepl("chr1p", nmf_alts$BEST_EFF_GENE) & nmf_alts$variant_classification_simple == "arm level loss"]
chr19q_loss <-nmf_alts$SAMPLE_ACCESSION_NBR[grepl("chr19q", nmf_alts$BEST_EFF_GENE) & nmf_alts$variant_classification_simple == "arm level loss"]
c1p19qloss <- chr19q_loss[chr19q_loss %in% chr1p_loss]

# Get rid of sample(s) that have multiple whole chromosome losses
c1p19qloss <- braf_arm_data$SAMPLE_ACCESSION_NBR[braf_arm_data$SAMPLE_ACCESSION_NBR %in% c1p19qloss & !(braf_arm_data$chr1q == "-1")]
c1p19qloss
classifier$WHO_Diagnosis[classifier$SAMPLE_ACCESSION_NBR %in% c1p19qloss]

# c7gain c10loss
chr7gains <- nmf_alts$SAMPLE_ACCESSION_NBR[grepl("chr7p", nmf_alts$BEST_EFF_GENE) & nmf_alts$variant_classification_simple == "arm level gain"]
chr7gains <- nmf_alts$SAMPLE_ACCESSION_NBR[grepl("chr7q", nmf_alts$BEST_EFF_GENE) & nmf_alts$variant_classification_simple == "arm level gain" & nmf_alts$SAMPLE_ACCESSION_NBR %in% chr7gains]
chr7gains

chr10loss <- nmf_alts$SAMPLE_ACCESSION_NBR[grepl("chr10p", nmf_alts$BEST_EFF_GENE) & nmf_alts$variant_classification_simple == "arm level loss"]
chr10loss <- nmf_alts$SAMPLE_ACCESSION_NBR[grepl("chr10q", nmf_alts$BEST_EFF_GENE) & nmf_alts$variant_classification_simple == "arm level loss" & nmf_alts$SAMPLE_ACCESSION_NBR %in% chr10loss]
chr10loss

c7gain10loss <- chr10loss[chr10loss %in% chr7gains]
c7gain10loss
"BL-15-E44047" %in% c7gain10loss

# EGFR gain
egfrgain <- nmf_alts$SAMPLE_ACCESSION_NBR[grepl("EGFR", nmf_alts$BEST_EFF_GENE) & nmf_alts$variant_classification == "gain"]
egfrgain
unique_egfrgain <- egfrgain[!(egfrgain %in% c7gain10loss)]
unique_egfrgain
#braf_data$variant_classification[grepl("EGFR", braf_data$BEST_EFF_GENE) & braf_data$SAMPLE_ACCESSION_NBR == "GENIE-MSK-P-0004119-T01-IM5"] #testing what the raw data look like

# TERT Promoter
tertpmut <- nmf_alts$SAMPLE_ACCESSION_NBR[grepl("TERT_Promoter", nmf_alts$BEST_EFF_GENE)]
tertpmut
tertpmut[tertpmut %in% c7gain10loss]
unique_tertpmut <- tertpmut[!(tertpmut %in% c7gain10loss)]
unique_tertpmut <- unique_tertpmut[!(unique_tertpmut %in% egfrgain)]
unique_tertpmut

## Combine cIMPACT and filter appropriately
cIMPACT_samps <- c(c7gain10loss, unique_egfrgain, unique_tertpmut)
cIMPACT_samps <- unique(cIMPACT_samps)
cIMPACT_samps ## n=55

# Make sure cIMPACT are 1p19q retained (avoids some oligos)
table(classifier$x1p_19q[classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps]) 
cIMPACT_samps <- classifier$SAMPLE_ACCESSION_NBR[classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps & (is.na(classifier$x1p_19q) | classifier$x1p_19q == 'FALSE')]
cIMPACT_samps ## n=53

# Make sure cIMPACT are IDH wt
table(classifier$IDH1.2[classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps])
classifier$SAMPLE_ACCESSION_NBR[classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps & classifier$IDH1.2==TRUE]
cIMPACT_samps <- classifier$SAMPLE_ACCESSION_NBR[classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps & (is.na(classifier$IDH1.2) | classifier$IDH1.2 == '0')]
cIMPACT_samps ## n=51

# Check what kinds of subtypes they are
table(classifier$subtype_path[classifier$membership %in% cIMPACT_samps])

# Only keep GBM, HGG, LGG (i.e. NOT PA or PXA)
cIMPACT_samps <- classifier$SAMPLE_ACCESSION_NBR[classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps & (grepl('GBM', classifier$subtype_path) | grepl('Astro', classifier$subtype_path) | grepl('HGG', classifier$subtype_path) | grepl('LGG', classifier$subtype_path))]
cIMPACT_samps ## n=45

# Assess characteristics of cIMPACT defined samples
table(classifier$subtype_path[classifier$membership %in% cIMPACT_samps])
table(classifier$IDH1.2[classifier$membership %in% cIMPACT_samps])
table(classifier$x1p_19q[classifier$membership %in% cIMPACT_samps])
table(classifier$H3K27M[classifier$membership %in% cIMPACT_samps])
table(classifier$Grade[classifier$membership %in% cIMPACT_samps])
table(classifier$WHO_Diagnosis[classifier$membership %in% cIMPACT_samps])

# Categorize cIMPACT3 status
classifier$cIMPACT <- ifelse(classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps, "1", "0")

# Assess which GBM are also cIMPACT or not
table(classifier$subtype[classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps])
classifier$WHO_Diagnosis[classifier$subtype == "GBM, IDH-wt" & !(classifier$SAMPLE_ACCESSION_NBR %in% cIMPACT_samps)]

# Classify subtypes by molecular definitions (ie GBM if met cIMPACT-3 criteria)
classifier$subtype <- classifier$subtype_path
table(classifier$subtype)
classifier$subtype[classifier$membership %in% cIMPACT_samps] <- "GBM, IDH-wt"
table(classifier$grade[classifier$cIMPACT == 1])
classifier$grade[classifier$membership %in% cIMPACT_samps] <- 4
table(classifier$grade[classifier$cIMPACT == 1])
table(classifier$subtype)

# If GBM, IDH-mt: make Astro, IDH-mt
classifier$subtype[classifier$SAMPLE_ACCESSION_NBR %in% (classifier$SAMPLE_ACCESSION_NBR[classifier$subtype_path == "GBM, IDH-wt" & classifier$IDH1.2 == 1])] <- "Astro, IDH-mt"

classifier$IDH_mutation[classifier$SAMPLE_ACCESSION_NBR %in% (classifier$SAMPLE_ACCESSION_NBR[classifier$subtype_path == "GBM, IDH-wt" & classifier$IDH1.2 == 1])]

classifier$IDH_mutation[classifier$subtype == "Astro, IDH-mt"]
classifier$IDH_mutation[classifier$IDH1.2 == '1']
                        
# If Oligo, IDH-wt: make HGG
classifier$subtype[grepl("ligo", classifier$WHO_Diagnosis) & classifier$IDH1.2 == '0'] # it already is

## Checking on grades and revising as needed
classifier$grade <- classifier$Grade
table(classifier$Cancer_Type_Specific[classifier$grade == "Low"])
table(classifier$Cancer_Type_Specific[classifier$grade == "Unknown"])

classifier$grade[classifier$SAMPLE_ACCESSION_NBR == "BL-15-G14321"]
classifier$grade[classifier$SAMPLE_ACCESSION_NBR == "BL-15-G14321"] <- "1" ## Based on Dr Bi's discussion with pathology; see BRAF_data_mastersheet.*.xlsx for details
classifier$grade[classifier$grade == "Low"] <- "2" 
classifier$grade[classifier$grade == "Unknown"] <- "2" 
table(classifier$grade)

## Classifying BRAF classes
table(nmf_alts_braf_only$BRAF_CLASS)
table(braf_muts_data$BRAF_CLASS)
unique(nmf_alts_braf_only$SAMPLE_ACCESSION_NBR[nmf_alts_braf_only$BRAF_CLASS == 'Other'])
unique(braf_muts_data$SAMPLE_ACCESSION_NBR[braf_muts_data$BRAF_CLASS == 'Other'])

braf_classes <- rbind(unique(nmf_alts_braf_only[c("SAMPLE_ACCESSION_NBR", "BRAF_CLASS")]), unique(braf_muts_data[c("SAMPLE_ACCESSION_NBR", "BRAF_CLASS")]))
braf_classes
table(braf_classes$BRAF_CLASS)
braf_classes[braf_classes$SAMPLE_ACCESSION_NBR %in% braf_classes$SAMPLE_ACCESSION_NBR[duplicated(braf_classes$SAMPLE_ACCESSION_NBR)],]

braf_classes_other <- unique(braf_classes[braf_classes$BRAF_CLASS == 'Other',])
braf_classes_noother <- unique(braf_classes[!(braf_classes$BRAF_CLASS == 'Other'),])
table(braf_classes_other$SAMPLE_ACCESSION_NBR %in% braf_classes_noother$SAMPLE_ACCESSION_NBR)
table(braf_classes_noother$SAMPLE_ACCESSION_NBR %in% braf_classes_other$SAMPLE_ACCESSION_NBR)

braf_classes_noother$SAMPLE_ACCESSION_NBR[(braf_classes_noother$SAMPLE_ACCESSION_NBR %in% braf_classes_other$SAMPLE_ACCESSION_NBR)]

# Collecting BRAF specific mutations as df's
braf_muts_spec_other <- unique(braf_muts_data[braf_muts_data$BRAF_CLASS == 'Other',])
braf_muts_spec_noother <- unique(braf_muts_data[!(braf_muts_data$BRAF_CLASS == 'Other'),])
braf_muts_spec_other[braf_muts_spec_other$SAMPLE_ACCESSION_NBR %in% braf_muts_spec_noother$SAMPLE_ACCESSION_NBR,]
braf_muts_spec_noother[braf_muts_spec_other$SAMPLE_ACCESSION_NBR %in% braf_muts_spec_noother$SAMPLE_ACCESSION_NBR,]

genie_braf_fusion_data <- genie_fusion_data[genie_fusion_data$Tumor_Sample_Barcode %in% all_samples,]
genie_braf_fusion_data$BRAF_Class <- classifier$braf_class[match(genie_braf_fusion_data$Tumor_Sample_Barcode, classifier$SAMPLE_ACCESSION_NBR)]
genie_braf_fusion_data <- genie_braf_fusion_data[!(grepl("Archer", genie_braf_fusion_data$Fusion)),]
genie_braf_fusion_data <- genie_braf_fusion_data[grepl("BRAF", genie_braf_fusion_data$Fusion),]

# Specifying specific BRAF mutation
classifier$braf_mutation <- braf_muts_spec_noother$BEST_EFF_PROTEIN_CHANGE[match(classifier$SAMPLE_ACCESSION_NBR, braf_muts_spec_noother$SAMPLE_ACCESSION_NBR)]
classifier$braf_mutation <- ifelse(is.na(classifier$braf_mutation), braf_muts_spec_other$BEST_EFF_PROTEIN_CHANGE[match(classifier$SAMPLE_ACCESSION_NBR, braf_muts_spec_other$SAMPLE_ACCESSION_NBR)],classifier$braf_mutation)

# Specifying BRAF rearrangements
classifier$braf_mutation <- ifelse(is.na(classifier$braf_mutation), dfci_fusion_data$BRAF_SV[match(classifier$SAMPLE_ACCESSION_NBR, dfci_fusion_data$SAMPLE_ACCESSION_NBR)], classifier$braf_mutation)
classifier$braf_mutation <- ifelse(is.na(classifier$braf_mutation), genie_braf_fusion_data$Fusion[match(classifier$SAMPLE_ACCESSION_NBR, genie_braf_fusion_data$Tumor_Sample_Barcode)], classifier$braf_mutation)

classifier$braf_mutation <- sub("p.", "", classifier$braf_mutation)
classifier$braf_mutation[classifier$braf_mutation == 'null'] <- NA

classifier$braf
## Old BRAF classification method
#classifier$braf_class <- nmf_alts_braf_only$BRAF_CLASS[match(classifier$SAMPLE_ACCESSION_NBR, nmf_alts_braf_only$SAMPLE_ACCESSION_NBR)]
#classifier$braf_class <- ifelse(classifier$SAMPLE_ACCESSION_NBR %in% braf_muts_data$SAMPLE_ACCESSION_NBR, braf_muts_data$BRAF_CLASS[match(classifier$SAMPLE_ACCESSION_NBR, braf_muts_data$SAMPLE_ACCESSION_NBR)], classifier$braf_class)

# New BRAF classification method
classifier$braf_class <- braf_classes_other$BRAF_CLASS[match(classifier$SAMPLE_ACCESSION_NBR, braf_classes_other$SAMPLE_ACCESSION_NBR)]
classifier$braf_class <- ifelse(classifier$SAMPLE_ACCESSION_NBR %in% braf_classes_noother$SAMPLE_ACCESSION_NBR, braf_classes_noother$BRAF_CLASS[match(classifier$SAMPLE_ACCESSION_NBR, braf_classes_noother$SAMPLE_ACCESSION_NBR)], classifier$braf_class)

classifier$braf_mutation <- ifelse(is.na(classifier$braf_mutation) & classifier$braf_class == "Other", nmf_alts_braf_only$variant_classification[match(classifier$SAMPLE_ACCESSION_NBR, nmf_alts_braf_only$SAMPLE_ACCESSION_NBR)], classifier$braf_mutation)

# Check BRAF class numbers
table(classifier$braf_class)
table(classifier$braf_class[!(classifier$age_int == '<18')]) # adults
table(classifier$braf_class[(classifier$age_int == '<18')]) # peds 

## Add overall survival info from updated censored clinical mastersheets
classifier$deceased <- censored.adult.braf.mastersheet$Alive[match(classifier$SAMPLE_ACCESSION_NBR, censored.adult.braf.mastersheet$SAMPLE_ACCESSION_NBR)]
classifier$deceased <- ifelse(is.na(classifier$deceased), censored.peds.braf.mastersheet$Alive[match(classifier$SAMPLE_ACCESSION_NBR, censored.peds.braf.mastersheet$SAMPLE_ACCESSION_NBR)], classifier$deceased)
classifier$deceased[classifier$deceased == "0"] <- "deceased"
classifier$deceased[classifier$deceased == "1"] <- "alive"
classifier$deceased[classifier$deceased == "deceased"] <- "1"
classifier$deceased[classifier$deceased == "alive"] <- "0"

table(is.na(classifier$deceased)) # how many we do or don't have alive/death info from
table(classifier$deceased)

classifier$os <- censored.adult.braf.mastersheet$OS[match(classifier$SAMPLE_ACCESSION_NBR, censored.adult.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

classifier$os <- ifelse(is.na(classifier$os), censored.peds.braf.mastersheet$OS[match(classifier$SAMPLE_ACCESSION_NBR, censored.peds.braf.mastersheet$SAMPLE_ACCESSION_NBR)], classifier$os)

classifier$os <- ifelse(is.na(classifier$os), genie.outcomes$Overall.Survival..Months.[match(classifier$SAMPLE_ACCESSION_NBR, genie.outcomes$SAMPLE_ACCESSION_NBR)], classifier$os)

table(is.na(classifier$os)) # how many we do or don't have os info from
classifier$deceased[is.na(classifier$os)]
classifier$os[is.na(classifier$deceased)]

table(classifier$age_int[!is.na(classifier$os)]) # n=179 --> n=187
 
### Old method of gene subsampling
### Remove genes that appear mutated in only one cohort
tcga_samples = classifier$membership[classifier$cohort == "TCGA"]
tcga_genes = unique(nmf_alts$BEST_EFF_GENE[nmf_alts$SAMPLE_ACCESSION_NBR %in% tcga_samples])

genie_samples = classifier$membership[classifier$cohort == "GENIE"]
genie_genes = unique(nmf_alts$BEST_EFF_GENE[nmf_alts$SAMPLE_ACCESSION_NBR %in% genie_samples])

dfci_samples = classifier$membership[classifier$cohort == "DFCI"]
dfci_genes = unique(nmf_alts$BEST_EFF_GENE[nmf_alts$SAMPLE_ACCESSION_NBR %in% dfci_samples])

jhh_samples = classifier$membership[classifier$cohort == "JHH"]
jhh_genes = unique(nmf_alts$BEST_EFF_GENE[nmf_alts$SAMPLE_ACCESSION_NBR %in% jhh_samples])

tcga_only = tcga_genes[!(tcga_genes %in% genie_genes) & !(tcga_genes %in% dfci_genes) & !(tcga_genes %in% jhh_genes)]
genie_only = genie_genes[!(genie_genes %in% tcga_genes) & !(genie_genes %in% dfci_genes) & !(genie_genes %in% jhh_genes)]
dfci_only = dfci_genes[!(dfci_genes %in% tcga_genes) & !(dfci_genes %in% genie_genes) & !(dfci_genes %in% jhh_genes)]
jhh_only = jhh_genes[!(jhh_genes %in% tcga_genes) & !(jhh_genes %in% genie_genes) & !(jhh_genes %in% dfci_genes)]

tcga_only
genie_only
dfci_only
jhh_only

cohort_spec_genes = c(tcga_only, genie_only, dfci_only, jhh_only)
cohort_spec_genes


genes_subset <- unique(nmf_alts$BEST_EFF_GENE)
unique(genes_subset)
genes_subset <- genes_subset[!(genes_subset %in% cohort_spec_genes)]
unique(genes_subset)

# unique(GENIE.gene.list.extended.trimmed$Hugo_Symbol)

# 
# # How many shared genes that are mutated
# unique(genes_subset[(genes_subset %in% tcga_genes) & (genes_subset %in% genie_genes) & (genes_subset %in% dfci_genes)]) # n=223
# 
# unique(genes_subset[(genes_subset %in% tcga_genes) & (genes_subset %in% genie_genes) & (genes_subset %in% dfci_genes) & (genes_subset %in% jhh_genes)]) # n=70
# 
# nmf_alts_subset <- nmf_alts[nmf_alts$BEST_EFF_GENE %in% genes_subset, ]
# nmf_alts_subset
# table(nmf_alts_subset$BEST_EFF_GENE[nmf_alts_subset$BEST_EFF_GENE %in% braf_arm_calls$BEST_EFF_GENE])
# 
# ##Genes in all 3 cohorts at BWH
# genes_in_all <- all_genes
# genes_in_all <- genes_in_all[(genes_in_all %in% tcga_genes) & (genes_in_all %in% genie_genes) & (genes_in_all %in% dfci_genes)]
# genes_in_all
# 
# ## Define a subset of samples from the above data that you can to make an oncoprint for
# #sub_samps <- classifier$membership[classifier$age_int %in% c("Young Adult", "Middle", "Old")]
# sub_samps <- classifier$membership
# #sub_samps <- unique(nmf_alts_subset$SAMPLE_ACCESSION_NBR)
# sub_samps
# 
# nmf_alts_subset <- nmf_alts
# nmf_alts_subset <- nmf_alts_subset[nmf_alts_subset$SAMPLE_ACCESSION_NBR %in% sub_samps, ]
# unique(nmf_alts_subset$SAMPLE_ACCESSION_NBR)
# unique(nmf_alts_subset$BEST_EFF_GENE)
# 
# ## Keep only genes that are altered in 5% of the samples (5% of 158 > 8 samples) or 5% of 248 >12
# genes_table <- nmf_alts_subset[,c("SAMPLE_ACCESSION_NBR", "BEST_EFF_GENE")]
# genes_table <- unique(genes_table)
# genes_table <- data.frame(table(genes_table$BEST_EFF_GENE))
# genes_table
# #genes_table <- genes_table[genes_table$Freq > 7,]
# genes_table <- genes_table[genes_table$Freq > 12,]
# genes_subset <- unique(genes_table$Var1)
# genes_subset
# unique(nmf_alts_subset$BEST_EFF_GENE[nmf_alts_subset$variant_classification_simple == "rearrangement" & !(nmf_alts_subset$BEST_EFF_GENE %in% genes_subset)])
# nmf_alts_subset <- nmf_alts_subset[nmf_alts_subset$BEST_EFF_GENE %in% genes_subset, ]
# unique(nmf_alts_subset$SAMPLE_ACCESSION_NBR)
# unique(nmf_alts_subset$BEST_EFF_GENE)
# 
# genes_subset
# genes_subset[!(genes_subset %in% IntersectGeneList)]
# IntersectGeneList[!(IntersectGeneList %in% genes_subset)]
# unique(nmf_alts$BEST_EFF_GENE[grepl("RAS",nmf_alts$BEST_EFF_GENE)])


### Optional: If you want to do smaller oncoprint of just the RAS-RAF-ERK pathway genes
#"BRAF" %in% ras_path_data$Gene.name
#genes_subset
#genes_subset <- genes_subset[genes_subset %in% ras_path_data$Gene.name]
#nmf_alts_subset <- nmf_alts_subset[nmf_alts_subset$BEST_EFF_GENE %in% genes_subset, ]
#unique(nmf_alts_subset$SAMPLE_ACCESSION_NBR)
#unique(nmf_alts_subset$BEST_EFF_GENE)

## Add genes to gene.lists that are missing
dfci_genes[!(dfci_genes %in% gene.list$Gene)]
dfci_missing_genes[!(dfci_missing_genes %in% GENIE.gene.list.extended.trimmed$Hugo_Symbol)]
dfci_missing_genes <- dfci_genes[!(dfci_genes %in% gene.list$Gene)]
dfci_missing_genes

gene.list$Gene
gene.list.extended <- gene.list
#gene.list.extended[dfci_missing_genes,] <- NA
gene.list.extended[dfci_missing_genes,] <- 0
gene.list.extended[dfci_missing_genes[grepl("chr",dfci_missing_genes)],] <- 1
gene.list.extended$Gene[rownames(gene.list.extended) %in% dfci_missing_genes] <- rownames(gene.list.extended[dfci_missing_genes,])

## Checking which panels for gene list extended dfci missing genes
dfci_missing_genes
for(g in dfci_missing_genes){
  #print(g)
  print(paste("Gene:", g))
  for(s in nmf_alts$SAMPLE_ACCESSION_NBR[nmf_alts$BEST_EFF_GENE == g]){
    #print(s)
    print(paste("...", (classifier$PANEL_VERSION[classifier$SAMPLE_ACCESSION_NBR == s])))
  }
  #print(nmf_alts$SAMPLE_ACCESSION_NBR[nmf_alts$BEST_EFF_GENE == g])
}
gene.list.extended[dfci_missing_genes,]

# Correct DFCI gene.list.extended genes
# From manual check
dfci.ext.genes.v2 = c("UTP23", "KIAA1549", "FAM131B", "MAP4K5", "AGK", "RALA", "lincRNA")
dfci.ext.genes.v3 = c("C17orf70", "ERCC6-PGBD3", "UTP23", "RAD51L3-RFFL", "MEF2BNB-MEF2B", "KIAA1549", "PDE7A", "EIF2SL3", "CNTN1", "NFAM", "INHBC", "FAM131B", "MAP4K5", "AGK", "CRLF1", "BCAS1", "PTPRZ1", "RALA", "lincRNA", "GIT2", "NEAT1", "GNAI1")

gene.list.extended$OncoPanel.v2 <- ifelse(gene.list.extended$Gene %in% dfci.ext.genes.v2, 1, gene.list.extended$OncoPanel.v2)
gene.list.extended$OncoPanel.v3 <- ifelse(gene.list.extended$Gene %in% dfci.ext.genes.v3, 1, gene.list.extended$OncoPanel.v3)

# Add JHH gene panel genes to gene.list.extended
gene.list.extended$STPv2 <- ifelse(gene.list.extended$Gene %in% jhh_genes_v2$Gene, 1, 0)
gene.list.extended$STPv3 <- ifelse(gene.list.extended$Gene %in% jhh_genes_v3$Gene, 1, 0)
gene.list.extended$STPv3.2 <- ifelse(gene.list.extended$Gene %in% jhh_genes_v3.2$Gene, 1, 0)
gene.list.extended$STPv4 <- ifelse(gene.list.extended$Gene %in% jhh_genes_v4$Gene, 1, 0)
gene.list.extended$STPL <- ifelse(gene.list.extended$Gene %in% jhh_genes_stpl$Gene, 1, 0)
gene.list.extended$NPP <- ifelse(gene.list.extended$Gene %in% jhh_genes_npp$Gene, 1, 0)

gene.list.extended[dfci_missing_genes,]

# GENIE genes to add
genie_missing_genes <- genie_genes[!(genie_genes %in% GENIE.gene.list.trimmed$Hugo_Symbol)]
genie_missing_genes
GENIE.gene.list.extended.trimmed <- GENIE.gene.list.trimmed
GENIE.gene.list.extended.trimmed[genie_missing_genes,] <- NA
GENIE.gene.list.extended.trimmed$Hugo_Symbol[rownames(GENIE.gene.list.extended.trimmed) %in% genie_missing_genes] <- rownames(GENIE.gene.list.extended.trimmed[genie_missing_genes,])
GENIE.gene.list.extended.trimmed$SEQ_ASSAY_ID[rownames(GENIE.gene.list.extended.trimmed) %in% genie_missing_genes] 
unique(GENIE.gene.list.extended.trimmed$SEQ_ASSAY_ID)
unique(master.sheet.glioma_all$PANEL_VERSION)
genie_missing_genes
nmf_alts_subset[nmf_alts_subset$BEST_EFF_GENE %in% genie_missing_genes, c("BEST_EFF_GENE","SAMPLE_ACCESSION_NBR")]
master.sheet.glioma_all$PANEL_VERSION[master.sheet.glioma_all$SAMPLE_ACCESSION_NBR %in% (nmf_alts_subset$SAMPLE_ACCESSION_NBR[nmf_alts$BEST_EFF_GENE %in% genie_missing_genes])]


## Compare gene panels
gene.list.extended$Gene[gene.list.extended$OncoPanel.v1 == 1 & gene.list.extended$OncoPanel.v2 == 1 & gene.list.extended$OncoPanel.v3 == 1 & gene.list.extended$STPv2 == 1 & gene.list.extended$STPv3 == 1 & gene.list.extended$STPv3.2 == 1 & gene.list.extended$STPv4 == 1 & gene.list.extended$STPL == 1]
dfci_genes[!dfci_genes %in% gene.list.extended$Gene[gene.list.extended$OncoPanel.v1 == 1 & gene.list.extended$OncoPanel.v2 == 1 & gene.list.extended$OncoPanel.v3 == 1]]

dfci_genes

## JHH gene panels
table(classifier$PANEL_VERSION)
jhh_genes[(jhh_genes %in% gene.list.extended$Gene)]
jhh_genes[!(jhh_genes %in% gene.list.extended$Gene)]

table(jhh_genes %in% unique(GENIE.gene.list.extended.trimmed$Hugo_Symbol))

unique(master.sheet.glioma_all$PANEL_VERSION[master.sheet.glioma_all$SAMPLE_ACCESSION_NBR %in% (nmf_alts_subset$SAMPLE_ACCESSION_NBR[nmf_alts$BEST_EFF_GENE %in% genie_missing_genes])])

## *** Look into not.covered genes further ***
## Not covered genes evaluation and correction
table(gene.list.extended$Gene[gene.list.extended$OncoPanel.v1 == 0] %in% not.covered$`1`)
gene.list.extended$Gene[(gene.list.extended$Gene[gene.list.extended$OncoPanel.v1 == 0] %in% not.covered$`1`)]
gene.list.extended$Gene[!(gene.list.extended$Gene[gene.list.extended$OncoPanel.v1 == 0] %in% not.covered$`1`)]

not.covered$`1`[grepl("IDH",not.covered$`1`)]

not.covered
not.covered.updated <- not.covered
"AGK" %in% not.covered.updated$`1`
not.covered.updated$`1` <- append("AGK", not.covered.updated$`1`)


## Subsample
all_samples
sub_samps <- all_samples
nmf_alts_subset <- nmf_alts
genes_subset
genes_subset[!(genes_subset %in% poor_coverage_genes)]
nmf_alts_subset <- nmf_alts_subset[nmf_alts_subset$BEST_EFF_GENE %in% genes_subset,]
nmf_alts_subset <- nmf_alts_subset[!(nmf_alts_subset$BEST_EFF_GENE %in% poor_coverage_genes),]

unique(nmf_alts_subset$BEST_EFF_GENE)
unique(nmf_alts_subset$SAMPLE_ACCESSION_NBR)

## This makes a long list of alterations in the form of a data frame with columns: Sample ID, Gene, variant classification (including "wt" for wildtype and "nc" for no coverage)
nmf_alts_input <- PlotComut(
  mut.maf1 = nmf_alts_subset, 
  samples = sub_samps,
  input.samples = "SAMPLE_ACCESSION_NBR", 
  input.genes = "BEST_EFF_GENE", input.mutations = "variant_classification_simple",
  dimensions = c(18,10), 
  col.vector = c("Unknown" = "white", "damaging mutation" = "skyblue", "focal loss" = "blue", "focal gain" = "red", "rearrangement" = "green", "missense" = "skyblue", "other mutation" = "orange",
                 "arm level gain" = "red", "arm level loss" = "blue","Glioblastoma" = "black", "Anaplasticoligo" = "green",
                 "Angiocentric" = "orange", "Astro" = "green", "Diffuseoligo" = "green", "Oligo" = "yellow", "oligooligo" = "yellow", "Pilocyticoligo" = "green", 
                 "Pediatric" = brewer.pal(4, "OrRd")[1], "Young Adult" = brewer.pal(4, "OrRd")[2], "Middle" = brewer.pal(4, "OrRd")[3], "Old" = brewer.pal(4, "OrRd")[4], 
                 "Other LGG" = "magenta", "ncc" == "white",
                 "Class I" = brewer.pal(9, "Set3")[1], "Class II" = brewer.pal(9, "Set3")[2], 
                 "Class III" = brewer.pal(9, "Set3")[3], "other BRAF" = brewer.pal(9, "Set3")[9]
                 ),
  #fixed.order = sample.order,
  #manual.order = c(gene.order),
  #file.path = "../Figures/Incidence Figures/Glioma Mega Comut Genes.pdf",
  y.axis.font.size = 8,
  coverage = TRUE,
  return.matrix = TRUE
  #phenotypes = classifier[, c("membership", "age_int", "subtype", "grade")], 
  #pheno_top = TRUE
)

## Define a subset of samples from the above data that you can to make an oncoprint for
#sub_samps <- classifier$membership[classifier$age_int %in% c("Young Adult", "Middle", "Old")]
sub_samps <- unique(nmf_alts_input$samples)
sub_genes <- unique(nmf_alts_input$genes)
#sub_genes <- IntersectGeneList
sub_genes
nmf_alts_input$genes

## Filter genes with low coverage
gene_coverage_df <- data.frame(table(nmf_alts_input$genes, nmf_alts_input$mutations == "not covered"))
gene_coverage_df <- gene_coverage_df[gene_coverage_df$Var2 == "TRUE",]
gene_coverage_df
names(gene_coverage_df) <- c("gene", "ignore", "not_covered")
nmf_alts
unique(gene_coverage_df$gene)
#poor_coverage_genes <- gene_coverage_df$gene[gene_coverage_df$not_covered > 148]
#poor_coverage_genes <- gene_coverage_df$gene[gene_coverage_df$not_covered > 148 | gene_coverage_df$not_covered == 135 | gene_coverage_df$not_covered == 65]
poor_coverage_genes <- gene_coverage_df$gene[gene_coverage_df$not_covered > 1]
poor_coverage_genes <- poor_coverage_genes[!(poor_coverage_genes %in% fusion_partners)]
poor_coverage_genes <- poor_coverage_genes[!(poor_coverage_genes %in% dfci_missing_genes)]
poor_coverage_genes <- poor_coverage_genes[!(poor_coverage_genes %in% genie_missing_genes)]
poor_coverage_genes <- poor_coverage_genes[!(poor_coverage_genes %in% c("TERT_Promoter"))]
#gene_coverage_df[gene_coverage_df$gene %in% ras_path_data$Gene.name,]
"BCAS1" %in% fusion_partners

table(nmf_alts_input$mutations)
table(gene_coverage_df$not_covered)
gene_coverage_df$gene[gene_coverage_df$not_covered == 0]

## Sub nmf alts input
sub_genes
#as.character(sub_genes[1:180])
as.character(sub_genes[1:45])
sub_nmf_alts_input <- nmf_alts_input
sub_nmf_alts_input$samples
sub_nmf_alts_input$genes
unique(sub_nmf_alts_input$genes)
sub_genes
select_genes <- c(as.character(sub_genes[1:45]), "chr1p")
select_genes <- select_genes[select_genes != "CNTNAP2"]
select_genes
#sub_nmf_alts_input <- sub_nmf_alts_input[sub_nmf_alts_input$samples %in% sub_samps & sub_nmf_alts_input$genes %in% sub_genes, ]
#sub_nmf_alts_input <- sub_nmf_alts_input[sub_nmf_alts_input$samples %in% sub_samps & sub_nmf_alts_input$genes %in% as.character(sub_genes[1:180]),]
sub_nmf_alts_input <- sub_nmf_alts_input[sub_nmf_alts_input$samples %in% sub_samps & sub_nmf_alts_input$genes %in% select_genes,]
sub_nmf_alts_input$samples <- sub_nmf_alts_input$samples[sub_nmf_alts_input$samples %in% sub_samps, drop = TRUE]
#sub_nmf_alts_input$genes <- sub_nmf_alts_input$genes[sub_nmf_alts_input$genes %in% as.character(sub_genes[1:180]), drop = TRUE]
sub_nmf_alts_input$genes <- sub_nmf_alts_input$genes[sub_nmf_alts_input$genes %in% select_genes, drop = TRUE]
#sub_nmf_alts_input <- nmf_alts_input[nmf_alts_input$genes %in% sub_genes, ]
sub_genes[!(sub_genes %in% sub_nmf_alts_input$genes)]
unique(sub_nmf_alts_input$samples)
unique(sub_nmf_alts_input$genes)

genes_order[!(genes_order %in% unique(sub_nmf_alts_input$genes))]

## Checking some things in the classifier
#table(classifier$braf_class[classifier$subtype == "GBM" & classifier$membership %in% sub_samps])
table(classifier$subtype[classifier$membership %in% sub_samps])
table(classifier$braf_class[classifier$membership %in% sub_samps])
table(classifier$age_int[classifier$membership %in% sub_samps])

## Check that you have 18 AND over and subset of samples is correct
adults <- classifier$membership[!(classifier$age_int == "<18")]
adults
classifier$SAMPLE_ACCESSION_NBR[is.na(classifier$age_num)]
classifier$Age_Interval[classifier$membership %in% adults[!(adults %in% sub_samps)]]
classifier$Age_Interval[classifier$membership %in% sub_samps[!(sub_samps %in% adults)]]
classifier$Age_Interval[classifier$membership %in% sub_samps]

## Get overall age info
median(classifier$age_num[classifier$SAMPLE_ACCESSION_NBR %in% adults])
min(classifier$age_num[classifier$SAMPLE_ACCESSION_NBR %in% adults])
max(classifier$age_num[classifier$SAMPLE_ACCESSION_NBR %in% adults])

median(all.braf.mastersheet$Age)
min(all.braf.mastersheet$Age)
max(all.braf.mastersheet$Age)

## Testing if there is a discrepancy in rearrangement alterations before/after PlotCoMut
table(nmf_alts_subset$variant_classification_simple)
table(nmf_alts_input$mutations)

## Find those classified as other
classifier$SAMPLE_ACCESSION_NBR[classifier$braf_class == 'Other']
unique(braf_muts_data$BEST_EFF_PROTEIN_CHANGE[braf_muts_data$BRAF_CLASS == 'Other'])
unique(braf_muts_data$BEST_EFF_VARIANT_CLASS[braf_muts_data$BRAF_CLASS == 'Other'])
unique(braf_muts_data$VARIANT_TYPE)

#nmf_alts_subset[nmf_alts_subset$SAMPLE_ACCESSION_NBR == "GENIE-MSK-P-0019039-T01-IM6" & nmf_alts_subset$BEST_EFF_GENE == "UBE2H",]
#nmf_alts_input[nmf_alts_input$samples == "GENIE-MSK-P-0019039-T01-IM6" & nmf_alts_input$genes == "UBE2H",]
#sub_nmf_alts_input[sub_nmf_alts_input$samples == "GENIE-MSK-P-0019039-T01-IM6" & sub_nmf_alts_input$genes == "UBE2H",]
#nmf_alts[nmf_alts$BEST_EFF_GENE == "UBE2H",]
#nmf_alts_input[nmf_alts_input$genes == "UBE2H",]

#before_func <- nmf_alts_subset[nmf_alts_subset$variant_classification_simple == "rearrangement",]
#after_func <- nmf_alts_input[nmf_alts_input$mutations == "rearrangement",]
###

## Format a wide version of the comut data
## Note: df.wide.mutations will be used for covariate analysis specifically, excluding genes with <5 events

# Check between input genes and Taibo's intersect gene list
IntersectGeneList[!(IntersectGeneList %in% colnames(df.wide.nmf_alts_input))]
temp_df <- nmf_alts[nmf_alts$BEST_EFF_GENE %in% IntersectGeneList[!(IntersectGeneList %in% colnames(df.wide.nmf_alts_input))],]


# Format df, incl as 'wide' horizontal format
df.wide.nmf_alts_input <- reshape(nmf_alts_input, v.names = "mutations", idvar = "samples", timevar = "genes", direction = "wide")
df.wide.mutations <- reshape(nmf_alts_input, v.names = "mutations", idvar = "samples", timevar = "genes", direction = "wide")

# Check types of mutations
#table(df.wide.nmf_alts_input$mutations.TP53)
table(df.wide.nmf_alts_input$TP53)
df.wide.mutations[, -1]
table(df.wide.mutations$mutations.TP53)

# Set 1 for mutation 0 for wt and NA for not covered genes
df.wide.nmf_alts_input <- as.data.frame(cbind(as.character(df.wide.nmf_alts_input[, 1]), apply(df.wide.nmf_alts_input[, -1], 2, helper)), stringsAsFactors = FALSE)
#df.wide.nmf_alts_input <- as.data.frame(cbind(as.character(df.wide.nmf_alts_input[, 1])), apply(df.wide.nmf_alts_input[, -1], 2, helper), stringsAsFactors = FALSE)
#df.wide.nmf_alts_input <- as.data.frame(cbind(as.numeric(df.wide.nmf_alts_input[, 1]), apply(df.wide.nmf_alts_input[, -1], 2, helper)), stringsAsFactors = FALSE)

df.wide.mutations <- as.data.frame(cbind(as.character(df.wide.mutations[, 1]), apply(df.wide.mutations[, -1], 2, helper)), stringsAsFactors = FALSE)
is.numeric(df.wide.mutations$mutations.TP53)

## Make numeric if needed
#df.wide.nmf_alts_input<-as.data.frame(cbind(as.data.frame(df.wide.nmf_alts_input[,1]), as.data.frame(sapply(df.wide.nmf_alts_input[,-1], as.numeric))))
df.wide.mutations <-as.data.frame(cbind(as.data.frame(df.wide.mutations[,1]), as.data.frame(sapply(df.wide.mutations[,-1], as.numeric))))

is.numeric(df.wide.mutations$mutations.TP53)
table(df.wide.mutations$mutations.TP53)

## Calculate column sums and subsetting df to exlcude anything <5 events
colSums(df.wide.mutations[,-1], na.rm=TRUE)

df.wide.muts <- df.wide.mutations[,-1][colSums(df.wide.mutations[,-1], na.rm=TRUE) > 4]
colSums(df.wide.muts, na.rm=TRUE)
df.wide.muts <- as.data.frame(cbind(as.data.frame(df.wide.mutations[,1]), as.data.frame(df.wide.muts)))

df.wide.muts <- df.wide.mutations[,(colSums(df.wide.mutations[,-1], na.rm=TRUE) > 4)]
colSums(df.wide.muts[,-1], na.rm=TRUE)
colnames(df.wide.mutations)[(colSums(df.wide.mutations[,-1], na.rm=TRUE) > 4)]

table(df.wide.mutations$mutations.ZRSR2)
colSums(df.wide.mutations[,2:ncol(df.wide.mutations)], na.rm=TRUE)
df.wide.muts <- df.wide.mutations[,(colSums(df.wide.mutations[,2:ncol(df.wide.mutations)], na.rm=TRUE) > 4)]

colSums(df.wide.mutations[,2:ncol(df.wide.mutations)], na.rm=TRUE) > 4

df.wide.mutations[,colSums(df.wide.mutations[,2:ncol(df.wide.mutations)], na.rm=TRUE) > 4]
colSums(df.wide.muts[,-1], na.rm=TRUE)
df.wide.muts[,1:(ncol(df.wide.muts)-3)]
ncol(df.wide.muts)

## Change colnames
colnames(df.wide.nmf_alts_input) <- c(gsub("mutations.", "", colnames(df.wide.nmf_alts_input)))
#colnames(df.wide.nmf_alts_input)[1] <- "SAMPLE_ACCESSION_NBR"
colnames(df.wide.nmf_alts_input)
table(df.wide.nmf_alts_input$KIAA1549)
IntersectGeneList[IntersectGeneList %in% colnames(df.wide.nmf_alts_input)]

colnames(df.wide.muts) <- c(gsub("mutations.", "", colnames(df.wide.muts)))
colnames(df.wide.muts)[1] <- "sample"


# This turns the above data frame into a list of matrices, each one representing a different type of alteration
mat_list <- NULL

for (i in 1:length(unique(sub_nmf_alts_input$mutations))){
  alt <- unique(sub_nmf_alts_input$mutations)[i]
  print(alt)
  alts_input <- sub_nmf_alts_input
  alts_input$mutations[alts_input$mutations != alt | alts_input$genes %in% c("quiet.mut", "cnv.quiet", "cnv.middle", "cnv.disrupted")] <- "na"
  alts_input <- rbind(alts_input[alts_input$mutations == alt, ], alts_input[alts_input$mutations == "na", ])
  test1 <- sapply(unique(alts_input$genes), function(x){length(unique(alts_input$samples[alts_input$genes == x & alts_input$mutations == alt]))})
  names(test1) <- unique(alts_input$genes)
  test1 <- sort(test1)
  test1 <- test1[!(names(test1) %in% c("CDKN2A/B", "PIK3CA/3R1", "MDM2/4","IDH1/2", "CDK4/6", "samples"))]
  #print(test1)
  print(table(alts_input$mutations[alts_input$genes == "BRAF"]))
  wide_alts <- reshape(alts_input, idvar = "samples", timevar = "genes",  direction = "wide")
  wide_alts[wide_alts == "na"] <- 0
  wide_alts[, -1][wide_alts[, -1] != 0] <- 1
  colnames(wide_alts) <- gsub("(mutations.)(*)", "\\2", colnames(wide_alts))
  print(table(wide_alts$BRAF))
  test2 <- sort(colSums(wide_alts[ ,!(colnames(wide_alts) %in% c("CDKN2A/B", "PIK3CA/3R1", "MDM2/4","IDH1/2", "CDK4/6", "samples"))]==1))
  print(test2)
  print(sum(test2==test1))
  print(test2[test2 != test1])
  wide_alts <- t(wide_alts)
  colnames(wide_alts) <- wide_alts[1, ]
  wide_alts <- wide_alts[-1, ]
  row.labels <- rownames(wide_alts)
  wide_alts <- apply(wide_alts, 2, as.numeric)
  rownames(wide_alts) <- row.labels
  if (i == 1){
    genes <- rownames(wide_alts)
    print(head(genes))
    samp.list <- colnames(wide_alts)
    print(head(samp.list))
  }
  wide_alts <- wide_alts[genes, samp.list]
  mat_list <- c(mat_list, list(wide_alts))
  names(mat_list)[i] <- alt
}

names(mat_list) <- gsub(" ", "_", names(mat_list))

## Assess mutational load & define hypermutated
table(nmf_alts_input$mutations)
#mut_num_data <- as.data.frame(table(sub_nmf_alts_input$samples[!(sub_nmf_alts_input$mutations == "wt") & !(sub_nmf_alts_input$mutations == "not covered")]))

#mut_num_data <- as.data.frame(table(sub_nmf_alts_input$samples[(sub_nmf_alts_input$mutations == "other") | (sub_nmf_alts_input$mutations == "damaging mutation")]))

mut_num_data <- as.data.frame(table(nmf_alts_input$samples[(nmf_alts_input$mutations == "other") | (nmf_alts_input$mutations == "damaging mutation")]))

mut_num_data

mut_num_data2 <- as.data.frame(table(nmf_alts_input$samples[!(nmf_alts_input$mutations == "not covered") & !(nmf_alts_input$mutations == "wt")]))

table(nmf_alts_input$mutations)
mut_num_data2
#mut_num_data<- as.data.frame(table(sub_nmf_alts_input$samples[(sub_nmf_alts_input$mutations == "other")]))

data.frame(table(mut_num_data$Freq)) %>%
  ggplot(aes(x=Var1, y=Freq)) +
  geom_point()
  #geom_bar(stat = 'identity', fill="#69b3a2") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
data.frame(table(mut_num_data$Freq))
quantile(mut_num_data$Freq, c(0.9, 0.925, 0.95, 0.975))

# Using segmented
dati <- data.frame(table(mut_num_data$Freq))
dati
dati$Var1 <- as.numeric(dati$Var1)
dati <- dati[2:30,]
dati$logTMB <- log10(as.numeric(dati$Var1))
dati <- dati[dati$logTMB > 0, ]
#dati$logFreq <- log10(dati$Freq)

# Plot
#ggplot(dati, aes(x = Var1, y = Freq)) +
ggplot(aes(y = Freq, x = Var1), data= dati) +
  geom_point() #+
#geom_line(data = dat2, color = 'blue')

# Fits
lfit <- lm(Freq ~ logTMB, data = dati)
lfit

sfit <- segmented(lfit, seg.Z = ~ logTMB)
exp10 <- function(x)10^x
exp10(sfit$psi)
summary(sfit)

bpoints <- 10^(sfit$psi)[, 2]
newx <- c(0.2, bpoints, 50)
newy <- predict(sfit, data.frame(logTMB = log10(newx)))

ggplot(aes(y = Freq, x = as.numeric(Var1)), data = dati) + 
  geom_path(data = data.frame(Var1 = newx, Freq = newy), colour = "grey 50", size = 2) +
  geom_point(size = 3) + 
  scale_x_log10()

## Define hypermutated samples

mut_num_data$Var1[mut_num_data$Freq >= 56]
mut_num_data$Var1[mut_num_data$Freq >= quantile(mut_num_data$Freq, 0.95)]
#hypermut_samps <- mut_num_data$Var1[mut_num_data$Freq >= 56]
hypermut_samps <- mut_num_data$Var1[mut_num_data$Freq >= quantile(mut_num_data$Freq, 0.95)]
hypermut_samps
classifier$braf_class[classifier$membership %in% hypermut_samps]

mut_num_data$Freq[match(classifier$SAMPLE_ACCESSION_NBR, mut_num_data$Var1)] 
classifier$mut_num <- mut_num_data$Freq[match(classifier$SAMPLE_ACCESSION_NBR, mut_num_data$Var1)] 
#classifier$mut_num_all <- mut_num_data2$Freq[match(classifier$SAMPLE_ACCESSION_NBR, mut_num_data2$Var1)] 
classifier$hypermutated <- "FALSE"
classifier$hypermutated[classifier$membership %in% hypermut_samps] <- "TRUE"
classifier$braf_class[classifier$hypermutated == "TRUE"]

quantile(classifier$mut_num, c(0.9, 0.925, 0.95, 0.975))
mut_num_data$Var1[mut_num_data$Freq >= quantile(mut_num_data$Freq, 0.95)]

## Going back to hypermutated genes and subsampling genes
hypermut_samps
(unique(nmf_alts_input$genes[nmf_alts_input$samples %in% hypermut_samps & nmf_alts_input$mutations == "other"]))

gene_coverage_df$not_covered[gene_coverage_df$gene %in% dfci_missing_genes]

gene_coverage_df[gene_coverage_df$gene %in% (unique(nmf_alts_input$genes[nmf_alts_input$samples %in% hypermut_samps & nmf_alts_input$mutations == "other"])),]

##Sorting samples by TMB
#mut_order <- mut_num_data$Var1[order(-mut_num_data$Freq)]
#mut_order
#sort(table(sub_nmf_alts_input$samples[!(sub_nmf_alts_input$mutations == "wt") & !(sub_nmf_alts_input$mutations == "not covered")]), decreasing = T)
#table(sub_nmf_alts_input$samples[!(sub_nmf_alts_input$mutations == "wt") & !(sub_nmf_alts_input$mutations == "not covered")])
#mut_order_old <- mut_order
#mut_order <- names(sort(table(sub_nmf_alts_input$samples[!(sub_nmf_alts_input$mutations == "wt") & !(sub_nmf_alts_input$mutations == "not covered")]), decreasing = T))
#mut_order <- names(sort(table(sub_nmf_alts_input$samples[(sub_nmf_alts_input$mutations == "other")]), decreasing = T))

### End hypermutation definition ###

## Define order of samples first by chromosomes
#sample_order_df <- sub_nmf_alts_input[sub_nmf_alts_input$samples %in% sub_samps,]
#sample_order_df <- sample_order_df[order(sub_nmf_alts_input$mutations[sub_nmf_alts_input$genes == "chr7q"]),]
#sample_order <- sample_order_df$samples
#sample_order

## Define order of samples by BRAF class 
#braf_order_df <- classifier[classifier$SAMPLE_ACCESSION_NBR %in% sub_samps,]
#braf_order_df$chr7q <- (sub_nmf_alts_input$mutations[sub_nmf_alts_input$genes == "chr7q"])[match(braf_order_df$SAMPLE_ACCESSION_NBR, sub_nmf_alts_input$samples)]
#braf_order_df$chr7q[is.na(braf_order_df$chr7q)] <- "x"

#braf_order_df <- braf_order_df[order(braf_order_df$chr7q, braf_order_df$braf_class),]
#classifier$braf_class <- factor(classifier$braf_class, levels = c("Class I", "Class II", "Class III", "Fusion", "Rearrangement", "Gain", "Other"))

## Keep original braf_class just in case!
# classifier$braf_class_original <- classifier$braf_class

## Add hypermutated classification for BRAF
classifier$braf_class <- ifelse(classifier$hypermutated == TRUE, "Hypermutated", classifier$braf_class)
table(classifier$braf_class)
table(classifier$braf_class_original)
table(classifier$subtype)

#braf_order_df <- classifier[classifier$SAMPLE_ACCESSION_NBR %in% sub_samps,]

## Order the samples
braf_order_df <- classifier
braf_order_df$braf_class <- factor(classifier$braf_class, levels = c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other"))
braf_order_df$subtype <- factor(classifier$subtype, levels = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA"))
braf_order_df$Pediatric <- factor(braf_order_df$Pediatric, levels = c("Pediatric", "Adult"))
braf_order_df$Pediatric[order(braf_order_df$Pediatric)]
braf_order_df$cIMPACT <- factor(braf_order_df$cIMPACT, levels = c("TRUE", "FALSE"))
#braf_order_df$hypermutated <- factor(braf_order_df$hypermutated, levels = c("TRUE", "FALSE"))
#braf_order_df$hypermutated[order(braf_order_df$hypermutated)]
#braf_order_df <- braf_order_df[with(braf_order_df, order(braf_order_df$braf_class, braf_order_df$Pediatric, -braf_order_df$mut_num, braf_order_df$cIMPACT, braf_order_df$subtype)),]
braf_order_df <- braf_order_df[with(braf_order_df, order(braf_order_df$braf_class, braf_order_df$Pediatric, -braf_order_df$mut_num_all, braf_order_df$subtype)),]
#braf_order_df <- braf_order_df[with(braf_order_df, order(braf_order_df$braf_class, braf_order_df$Pediatric, braf_order_df$cIMPACT, -braf_order_df$mut_num, braf_order_df$subtype)),]

braf_order <- braf_order_df$SAMPLE_ACCESSION_NBR 
braf_order
#braf_order_old <- braf_order
## Define order of genes by frequency mutated
#genes_order <- names(sort(table(nmf_alts$BEST_EFF_GENE), decreasing = T))
#genes_order <- names(sort(table(nmf_alts_input$genes), decreasing = T))
#names(sort(table(nmf_alts_input$genes), decreasing = T))

## Ordering genes

# Gains
nmf_alts_braf_gains <- nmf_alts[nmf_alts$SAMPLE_ACCESSION_NBR %in% classifier$membership[classifier$braf_class == "Gain"],]
sort(table(sub_nmf_alts_input$genes[sub_nmf_alts_input$mutations == "focal gain"]))
gain_genes <- c("SMO", "EZH2", "AGK", "FAM131B", "UBE2H", "CNTNAP2", "CDK4/6", "MET")
#altered_genes <- table(sub_nmf_alts_input$genes[!(sub_nmf_alts_input$mutations == "wt") & !(sub_nmf_alts_input$mutations == "not covered")])
#start_order <- c("BRAF", "chr7p", "chr7q", "chr10p", "chr10q", "chr1p", "chr19q")

start_order <- c("BRAF", "chr7p", "chr7q", "chr10p", "chr10q", "chr1p", "chr19q", "TERT_Promoter", "EGFR", "IDH1/2")
unique(sub_nmf_alts_input$genes)
#fusion_partners <- unique(sub_nmf_alts_input$genes[sub_nmf_alts_input$mutations == "rearrangement" & !(sub_nmf_alts_input$genes %in% start_order)])
fusion_partners <- unique(fusion_sankey_df$Fusion_Partner)
fusion_partners <- as.character(fusion_partners)
gain_genes[gain_genes %in% sub_nmf_alts_input$genes]
fusion_partners[fusion_partners %in% sub_nmf_alts_input$genes]
start_order <- c(start_order, fusion_partners[fusion_partners %in% unique(sub_nmf_alts_input$genes, gain_genes[gain_genes %in% unique(sub_nmf_alts_input$genes)])])
start_order
#ras_path_genes <- ras_path_data$Gene.name[!(ras_path_data$Gene.name %in% start_order)]
#start_order <- c(start_order, ras_path_genes)

ordering_genes <- unique(sub_nmf_alts_input$genes[!(sub_nmf_alts_input$genes %in% start_order), drop = TRUE])
ordering_genes <- levels(ordering_genes)[ordering_genes]
ordering_genes
#ordering_genes[1:164]
#altered_genes <- table(sub_nmf_alts_input$genes[sub_nmf_alts_input$genes %in% ordering_genes & !(sub_nmf_alts_input$mutations == "wt") & !(sub_nmf_alts_input$mutations == "not covered"), drop = TRUE])
#altered_genes
#names(sort(altered_genes,decreasing = T))
#genes_order <- c(start_order, names(sort(altered_genes, decreasing = T)))
genes_order <- c(start_order, ordering_genes)
genes_order <- unique(genes_order)
genes_order[(!genes_order %in% unique(sub_nmf_alts_input$genes))]
unique(sub_nmf_alts_input$genes)
genes_order
## Define all the colors and shapes for each type of alteration
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "beige", col = NA))
  },
    focal_loss = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                         gp = gpar(fill = col["focal_loss"], col = NA)),
  damaging_mutation = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                              gp = gpar(fill = col["damaging_mutation"], col = NA)),
  other = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                     gp = gpar(fill = col["other"], col = NA)),
  missense = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                              gp = gpar(fill = col["missense"], col = NA)),
  focal_gain = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                              gp = gpar(fill = col["focal_gain"], col = NA)),
  not_covered = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                              gp = gpar(fill = col["not_covered"], col = NA)),
  wt = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                      gp = gpar(fill = col["wt"], col = NA)),
  arm_level_gain = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                              gp = gpar(fill = col["arm_level_gain"], col = NA)),
  arm_level_loss = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["arm_level_loss"], col = NA)),
  test = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["test"], col = NA)),
  cnv_info = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                                  gp = gpar(fill = col["cnv_info"], col = NA)),
  rearrangement = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, 
                                                 gp = gpar(fill = col["rearrangement"], col = NA)),
  BRAFother = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, 
                                        gp = gpar(fill = col["BRAFother"], col = NA))
  )

## Ignore this -- This is where I defined a predetermined order for the genes in the oncoprint
# siggenes_sig <- siggenes[siggenes$qvalue < 0.1,]
# siggenes_sig$V1 <- gsub("([A-Z0-9])\\.", "\\1/", siggenes_sig$V1)
# siggenes_sig <- siggenes_sig[order(siggenes_sig$qvalue), ]
# siggenes_sig$V1 <- gsub("_[gl][ainoss]*", "", siggenes_sig$V1)
# siggenes_sig$V1 <- gsub("([A-Z0-9])(\\.)([A-Z0-9])", "\\1/\\3", siggenes_sig$V1)
# siggenes_sig$V1[siggenes_sig$V1 == "MDM2/MDM4"] <- "MDM2/4"
# siggenes_sig <- rbind(siggenes_sig[siggenes_sig$sign_ratio == "Enriched", ],siggenes_sig[siggenes_sig$sign_ratio == "Exclusive", ])
# gene.order <- siggenes_sig[!duplicated(siggenes_sig$V1), ]
# cluster.genes <- gene.order$cluster[order(gene.order$cluster)]
# gene.order <- gene.order$V1[order(gene.order$cluster)]
# gene.order[gene.order == "7"] <- "7p"
# gene.order <- c(gene.order[1:(which(gene.order=="7p"))], "7q", gene.order[(which(gene.order=="7p")+1):length(gene.order)])
# cluster.genes <- c(cluster.genes[1:(which(gene.order=="7p"))], cluster.genes[which(gene.order=="7p")], cluster.genes[(which(gene.order=="7p")+1):length(cluster.genes)])
# 
# gene.order <- c(gene.order, unique(nmf_alts_input$genes)[!(unique(nmf_alts_input$genes) %in% gene.order)])
# gene.order <- data.frame(gene = gene.order, cluster = cluster.genes)


## Define legend colors, and which traits you want to appear in the barplot on the left side
col = c("blue", "red", "green", "skyblue", "orange", "beige", "white", "red", "blue", "purple","orange")
names(col) <- c("focal_loss", "focal_gain", "rearrangement", "damaging_mutation", "other", "wt", "not_covered", "arm_level_gain","arm_level_loss", "cnv_info", "test")

barplot_traits <- c("focal_loss", "focal_gain", "rearrangement", "arm_level_gain","arm_level_loss", "damaging_mutation", "other")
#barplot_traits <- c("focal_loss", "focal_gain", "rearrangement","damaging_mutation", "other")

## Below is the code for making clinical annotations along the top or bottom of the oncoprint
#ano_col = c("GBM" = "black", "Oligo" = "yellow", "Astro" = "green", "Other astro" = "lightsalmon", "Other" = "gray", "PA" = "skyblue", "Other High Grade" = "red", "PXA" = "violet", "Other glioma"= "purple")
#ano_col = c("GBM" = "black", "Oligo" = "yellow", "HGG" = "lightsalmon", "LGG" = "green", "Other" = "gray", "PA" = "skyblue", "Other High Grade" = "red", "PXA" = "violet", "Other glioma"= "purple")
ano_col = c("GBM, IDH-wt" = "black", "Oligo" = "lightsalmon", "HGG" = "red", "Astro, IDH-mt" = "green", "LGG" = "yellow", "Other" = "gray", "PA" = "skyblue", "Other High Grade" = "red", "PXA" = "violet", "Other glioma"= "purple")

table(classifier$subtype)
df <- classifier[match(colnames(mat_list[[1]]), classifier$membership),]
df$age <- df$age_int
df$class <- df$braf_class
#df$age[df$age == "?"] <- NA
#df <- df[ , c("primary", "grade", "subtype", "class")] #, "primary")]
#df <- df[ , c("age", "grade", "subtype", "class")] #, "primary")]
df <- df[ , c("age", "subtype", "class")] #, "primary")]
row.labels <- df$membership
rownames(df) <- row.labels
#df <- df[rownames(df) %in% sub_samps, ]
df

age_col <- brewer.pal(8, "OrRd")[c(1,3,5,7)]
names(age_col) <- c("<18", "18-34", "35-50", ">50")
#hyper_col <- brewer.pal(10, "Paired")[10:9]
#names(hyper_col) <- c("TRUE", "FALSE")
grade_col <- brewer.pal(7, "YlGnBu")[c(1,3,5,7)]
names(grade_col) <- c("1", "2", "3", "4")

#braf_col <- c(brewer.pal(9,"Set1")[1:3], brewer.pal(9, "Set1")[4], brewer.pal(9,"Set1")[9]) 
#names(braf_col) <- c("Class I", "Class II", "Class III", "fusion", "other BRAF")

#braf_col <- c("Class I" = "#F8766D", "Class II" = "#619CFF", "Class III" = "#00BA38", "Fusion" = "#C77CFF", "Other" = "#D9D9D9")

#braf_col <- c(brewer.pal(9,"Set1")[1:6], brewer.pal(9,"Set1")[9])
#names(braf_col) <- c("Gain", "Class I", "Class II", "Class III", "Fusion", "Fusion2","Other")

#braf_col <- c("Class I" = "#C84E4C", "Class II" = "#E8852C", "Class III" = "#CD9B27", "Fusion" = "#449F43", "Fusion2" = "#2383B2", "Gain" = "#945B9D", "Other" = "#D9D9D9")

#braf_col <- c("Class I" = "#4D8CE8", "Class II" = "#6ABFB6", "Class III" = "#9BD478", "Fusion" = "#D7D452", "Rearrangement" = "#FFB341", "Gain" = "#FF6131", "Other" = "#D9D9D9")

#Classic
#braf_col <- c("Class I" = "#4D8CE8", "Class II" = "#9BD478", "Class III" = "#D7D452", "Fusion" = "#FFB341", "Gain" = "#FF6131", "Other" = "gray")

# Spectral
braf_col <- c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray50", "gray")
names(braf_col) <- c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other")

#Sankey Taibo greens
#braf_col <- c("Class I" = "#00441b", "Class II" = "#238b45", "Class III" = "#66c2a4", "Fusion" = "#fb8072", "Gain" = "#80b1d3", "Other" = "#ffffb3")

colnames(df) <- capitalize(colnames(df))
df
## Order the plots
df$Age <- factor(df$Age, levels = c("<18", "18-34", "35-50", ">50"))
#df$Hypermutated <- factor(df$Hypermutated, levels = c("TRUE", "FALSE"))
#df$Grade <- factor(df$Grade, levels = c("1", "2", "3", "4"))
df$Subtype <- factor(df$Subtype, levels = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA"))
df$Class <- factor(df$Class, levels = c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other"))
df
ha_top = HeatmapAnnotation(
  df = df,
  #col = list(molecular = ano_col, grade = grade_col, primary = prim_col),
  #col = list(Molecular = ano_col, Grade = grade_col, Age = age_cols),
  #col = list(Primary = prim_col, Grade = grade_col, Subtype = ano_col, Class = braf_col),
  #col = list(Age = age_col, Grade = grade_col, Subtype = ano_col, Class = braf_col),
  col = list(Age = age_col, Subtype = ano_col, Class = braf_col),
  #na_col = "grey"
  na_col = "white"
  )
ha_top

df2 <- classifier[match(colnames(mat_list[[1]]), classifier$membership),]
#df2$age <- classifier$Age_Interval[match(df2$membership, classifier$SAMPLE_ACCESSION_NBR)]
#df2$age <- classifier$age_int[match(df2$membership, classifier$SAMPLE_ACCESSION_NBR)]
df2$primary <- classifier$Primary[match(df2$membership, classifier$SAMPLE_ACCESSION_NBR)]
df2$primary[df2$primary == 0 & !is.na(df2$primary)] <- "FALSE"
df2$primary[df2$primary == 1 & !is.na(df2$primary)] <- "TRUE"
df2$cohort <- classifier$cohort[match(df2$membership, classifier$SAMPLE_ACCESSION_NBR)]
df2 <- df2[,c("primary", "cohort")]
df2
#df2$age[df2$age == "?"] <- NA
#age_cols <- brewer.pal(6, "OrRd")[c(1,3,5)]
#names(age_cols) <- c("Young Adult", "Middle", "Old")
#cohort_cols <- c("DFCI" = "slateblue", "TCGA" = "orangered", "GENIE"="springgreen")
#cohort_cols <- c("DFCI" = "#F8766D", "TCGA" = "#00BA38", "GENIE"="#619CFF")
#cohort_cols <- c("DFCI" = "#E7B800", "TCGA" = "#2E9FDF", "GENIE"="#FC4E07")
prim_col <- brewer.pal(4, "Paired")[3:4]
names(prim_col) <- c("TRUE", "FALSE")
cohort_cols <- brewer.pal(4, "Set2")
names(cohort_cols) <- c("GENIE", "DFCI", "TCGA", "JHH")
rownames(df2) <- row.labels
#df2 <- df2[rownames(df2) %in% sub_samps, ]
colnames(df2) <- capitalize(colnames(df2))
df2$Primary <- factor(df2$Primary, levels = c("TRUE", "FALSE"))
#df2$Age <- factor(df2$Age, levels = c("Young Adult", "Middle", "Old"))

ha_bottom = HeatmapAnnotation(
  df = df2[ , c("Primary", "Cohort")],
  col = list(Primary = prim_col, Cohort = cohort_cols),
  #na_col = "grey"
  na_col = "white"
 )
## Make the oncoprint
mat_list[c("focal_gain", "focal_loss", "rearrangement", "damaging_mutation", "other", "arm_level_gain", "arm_level_loss")]
mat_list

#ht <- oncoPrint(mat_list[c("focal_gain", "focal_loss", "rearrangement", "damaging_mutation", "other", #"arm_level_gain", "arm_level_loss")], show_pct = F,
ht <- oncoPrint(mat_list, show_pct = F,
                row_names_side = "left", pct_side = "right",
                row_names_gp = gpar(size = 8, fontsize = 8),
          alter_fun = alter_fun, col = col, 
          column_split = classifier$braf_class[match(samp.list, classifier$membership)], 
          #column_split = braf_order_df$braf_class[match(samp.list, $membership)], 
          column_km = 1, column_title = NULL,
          #column_split = classifier$basis[match(samp.list, classifier$membership)], 
          #row_order = gene.order$gene, 
        #cluster_columns = hc_order,
        column_order = braf_order,
        row_order = genes_order,
        top_annotation = ha_top, 
            #c(HeatmapAnnotation(
            #column_barplot = anno_oncoprint_barplot(barplot_traits[barplot_traits != "cnv_info"]#, ylim=c(0,50) #border = TRUE, # only MUT
                                                    #height = unit(4, "cm")
            #                                        )), ha_top),
          #row_split = gene.order$cluster[match(genes, gene.order$gene)],
          right_annotation = NULL,
          #  rowAnnotation(
          #row_barplot = anno_oncoprint_barplot(barplot_traits,  # only AMP and HOMDEL
          #                                     border = TRUE, #height = unit(4, "cm"),
          #                                axis_param = list(side = "bottom", labels_rot = 45))),
          bottom_annotation = ha_bottom
          )

print(ht)

##Save as a pdf, if desired
pdf("NMF/NMF some info/NMF k5 oncoprint molec.pdf", width = 16, height = 8)
ht
dev.off()

## Files to write
#write.csv(classifier, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/clinical-master-sheets/braf_all_glioma_classifier_datasheet.censored.csv")

#write.csv(dfci_braf_genomic_specimen, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/dfci_braf_genomic_specimen.csv")
#write.csv(classifier[classifier$age_int %in% c("Young Adult", "Middle", "Old"),], "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/clinical-master-sheets/braf_adult_glioma_classifier_datasheet.censored.csv")

#classifier$braf_class <- braf_classes_other$BRAF_CLASS[match(classifier$SAMPLE_ACCESSION_NBR, braf_classes_other$SAMPLE_ACCESSION_NBR)]
#classifier$braf_class <- ifelse(classifier$SAMPLE_ACCESSION_NBR %in% braf_classes_noother$SAMPLE_ACCESSION_NBR, braf_classes_noother$BRAF_CLASS[match(classifier$SAMPLE_ACCESSION_NBR, braf_classes_noother$SAMPLE_ACCESSION_NBR)], classifier$braf_class)

#write.csv(temp_df, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/braf_fusion_partners.csv")

#write.csv(temp_df, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/molecular/braf_dfci_fusion_sv_comments.csv")


## Just a place to check things

# Checking fusions (again)
dfci.req.genomic.sv[dfci.req.genomic.sv$SAMPLE_ACCESSION_NBR %in% all_samples,] # none of the study samples are in this SV datafile!
dfci.req.genomic.mutations[dfci.req.genomic.mutations$SAMPLE_ACCESSION_NBR %in% all_samples,] # all SNVs only, no SVs

temp_df <- adult.braf.mastersheet[adult.braf.mastersheet$BRAF_Class == "Hypermutated",]
temp_df$SAMPLE_ACCESSION_NBR

hypermut_samps

table(temp_df$BRAF_Mutation_Specific)
table(temp_df$IDH_Mutation)
adult_gbm_fusion_samps <- temp_df$SAMPLE_ACCESSION_NBR
adult_gbm_fusion_samps
temp_df <- df.wide.muts[df.wide.muts$sample %in% adult_classI_samps,]
table(temp_df$IDH_Mutation)

temp_df <- df.wide.muts[df.wide.muts$sample %in% adult_gbm_fusion_samps,]
classifier$braf_class[classifier$SAMPLE_ACCESSION_NBR %in% df.wide.muts$sample[df.wide.muts$HRAS == "1"]]
table(all.braf.mastersheet$BRAF_Class)
table(peds.braf.mastersheet$BRAF_Mutation_Specific[peds.braf.mastersheet$BRAF_Class == "Fusion"])


adult.braf.mastersheet$WHO_Grade
adult.braf.mastersheet[adult.braf.mastersheet$BRAF_Class == "Fusion", c("BRAF_Mutation_Specific", "WHO_Grade")]

peds.braf.mastersheet[peds.braf.mastersheet$BRAF_Class == "Fusion", c("BRAF_Mutation_Specific", "WHO_Grade")]


mean(adult.braf.mastersheet$Age)
table(peds.braf.mastersheet$Age[peds.braf.mastersheet$Cohort == "DFCI"])
peds.braf.mastersheet$Age
table(peds.braf.mastersheet$Cohort)

classIII_samps.adult <- adult.braf.mastersheet$SAMPLE_ACCESSION_NBR[adult.braf.mastersheet$BRAF_Class == "Class III"]
classIII_samps <- all.braf.mastersheet$SAMPLE_ACCESSION_NBR[all.braf.mastersheet$BRAF_Class == "Class III"]
classIII_samps
temp_df <- nmf_alts[nmf_alts$SAMPLE_ACCESSION_NBR %in% classIII_samps.adult,]
temp_df2 <- nmf_alts[nmf_alts$SAMPLE_ACCESSION_NBR %in% classIII_samps,]

fusion_samps <- all.braf.mastersheet$SAMPLE_ACCESSION_NBR[all.braf.mastersheet$BRAF_Class == "Fusion"]
temp_list <- unique(nmf_alts$SAMPLE_ACCESSION_NBR[nmf_alts$SAMPLE_ACCESSION_NBR %in% fusion_samps & nmf_alts$BEST_EFF_GENE == "CDKN2A/B"])

adult.braf.mastersheet$Subtype[adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% temp_list]
