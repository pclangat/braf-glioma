## Pinky Langat 2021
## This code is for covariate analysis: univariate and multivariate Cox analysis 
## HRs for overall survival; will plot forest plot
## This code is designed to run after braf_oncoprint.R
## Has supporting analyses

library(tidyverse)
library(tidytidbits)
library(survivalAnalysis) 

## Setup integrated outcomes/genomic dataframe ##
df_forest <- classifier[!is.na(classifier$os),]
df_forest
df_forest$deceased <- as.numeric(df_forest$deceased)
df_forest$grade <- as.numeric(df_forest$grade)
df_forest$IDH1.2 <- ifelse(df_forest$IDH_mutation == "Canonical", 1, 0)

## Add other features from mastersheet
df_forest$alive <- all.braf.mastersheet$Alive[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]
df_forest$alive <- as.numeric(df_forest$alive)

df_forest$Tumor_Location <- all.braf.mastersheet$Tumor_Location[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

df_forest$Tumor_Laterality <- all.braf.mastersheet$Tumor_Laterality[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

df_forest$Extent_Resection <- all.braf.mastersheet$Extent_Resection[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

df_forest$BRAF_Inhibitor <- all.braf.mastersheet$BRAF_Inhibitor[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

df_forest$MGMT <- all.braf.mastersheet$MGMT[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

df_forest$c7Gain_10Loss <- all.braf.mastersheet$c7Gain_10Loss[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

df_forest$EGFR_Gain <- all.braf.mastersheet$EGFR_Gain[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

df_forest$TERTp_Mut <- all.braf.mastersheet$TERTP_Mut[match(df_forest$SAMPLE_ACCESSION_NBR, all.braf.mastersheet$SAMPLE_ACCESSION_NBR)]

df_forest$CDKN2AB_Loss <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "CDKN2A/B",]$variant_classification)[match(df_forest$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "CDKN2A/B",]$SAMPLE_ACCESSION_NBR)]
df_forest$CDKN2AB_Loss <- ifelse(df_forest$CDKN2AB_Loss == "loss", 1, 0)
df_forest$CDKN2AB_Loss[is.na(df_forest$CDKN2AB_Loss)] <- 0
df_forest$CDKN2AB_Loss

## Setup remaining binary gene alteration data
df_forest <- cbind(df_forest, df.wide.muts[match(df_forest$SAMPLE_ACCESSION_NBR, df.wide.muts$sample),])
#df_forest <- cbind(df_forest, df.wide.nmf_alts_input[match(df_forest$SAMPLE_ACCESSION_NBR, df.wide.nmf_alts_input$V1),])
#df_forest <- cbind(df_forest, df.wide.nmf_alts_input[match(df_forest$SAMPLE_ACCESSION_NBR, df.wide.nmf_alts_input$V1),])
#df_forest <- df_forest[, !duplicated(colnames(df_forest))]

colnames(df_forest)

## Do the filtering again for genes altered in at least 5 patients of this integrated analysis:
## Calculate column sums and subsetting df to exlcude anything <5 events (or at least just make a list of these genes)
colnames(df_forest)
colSums(df_forest[,c("EGFR", "PDGFRB")])
colSums(df_forest[,52:ncol(df_forest)])
colSums(df_forest[,52:ncol(df_forest)])[colSums(df_forest[,52:ncol(df_forest)]) < 5]

list.integratedlessthan5 <- c("SMARCB1", "CHEK2", "ERCC3", "MEN1", "NF2", "SH2B3", "BCL6", "PDGFRB", "ARAF", "AURKA", "ERCC4", "MUTYH", "NFE2L2", "BCL2", "CCNE1")

## Check overlap between these genes and Taibo's intersect list
IntersectGeneList[(IntersectGeneList %in% colnames(df_forest))]
IntersectGeneList[!(IntersectGeneList %in% colnames(df_forest))]
colnames(df_forest)[!(colnames(df_forest) %in% IntersectGeneList)]
temp_df <- nmf_alts[nmf_alts$BEST_EFF_GENE %in% IntersectGeneList[!(IntersectGeneList %in% colnames(df_forest))],]

"PARK2" %in% IntersectGeneList

temp_df <- df_forest[,c("SAMPLE_ACCESSION_NBR", "TERTp_Mut", "TERT_Promoter")]
#temp_df <- df.wide.nmf_alts_input[match(df_forest$SAMPLE_ACCESSION_NBR, df.wide.nmf_alts_input$V1),]


## Quick check that in proper order
#temp_df <- df_forest[,c("SAMPLE_ACCESSION_NBR", "V1", "PANEL_VERSION", "TACC3")]

## Subsets
df_forest.adult <- df_forest[!(df_forest$age_int == "<18"),]
df_forest.adult_gbm <- df_forest.adult[df_forest.adult$subtype == "GBM, IDH-wt",]
df_forest.adult_V600 <- df_forest.adult[df_forest.adult$braf_class == "Class I",]
df_forest.PA <- df_forest[df_forest$subtype == "PA",]

## Overview stats incl median survival
df_forest %>%
  analyse_survival(vars(os, deceased))

## Plot formatting
default_args <- list(break.time.by="breakByMonth",
                     xlab=".OS.years",
                     #legend.title="ECOG Status",
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean",
                     #palette = c("green",  "black", "red", "yellow", "violet", "lightsalmon", "skyblue" ),
                     palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray"),
                     #ggtheme=ggplot2::theme_bw(10))
                     ggtheme=ggplot2::theme_classic(10))
#scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "red", "green", "yellow", "lightsalmon", "skyblue", "violet")),
## Simple survival comparisons and OS curve grid plottings
df_forest %>%
  #mutate(BRAF=recode_factor(braf_class,'Other'="Other", 'Class I'= "Class I", 'Class II'="Class II/III", 'Class III'="Class II/III", 'Fusion'="Fusion",'Gain'="Gain",'Hypermutated'="Hypermutated")) %>%
  mutate(BRAF=recode_factor(braf_class,'Class I'= "Class I", 'Class II'="Class II", 'Class III'="Class III", 'Fusion'="Fusion",'Gain'="Gain",'Other'="Other",'Hypermutated'="Other")) %>%
  analyse_survival(vars(os, deceased), by=BRAF) -> result
result

braf_surv.adults %>%
  #mutate(BRAF=recode_factor(braf_class,'Other'="Other", 'Class I'= "Class I", 'Class II'="Class II/III", 'Class III'="Class II/III", 'Fusion'="Fusion",'Gain'="Gain",'Hypermutated'="Hypermutated")) %>%
  mutate(BRAF=recode_factor(BRAF_Class_Groups,'Class I'= "Class I", 'Class II'="Class II", 'Class III'="Class III", 'Fusion'="Fusion",'Gain'="Gain",'Other'="Other",'Hypermutated'="Other")) %>%
  analyse_survival(vars(OS, Death), by=BRAF) -> result

df_forest.adult %>%
  #mutate(BRAF=recode_factor(braf_class,'Other'="Other", 'Class I'= "Class I", 'Class II'="Class II/III", 'Class III'="Class II/III", 'Fusion'="Fusion",'Gain'="Gain",'Hypermutated'="Hypermutated")) %>%
  mutate(BRAF=recode_factor(braf_class, 'Class I'= "Class I", 'Class II'="Class II", 'Class III'="Class III", 'Fusion'="Fusion",'Gain'="Gain", 'Other'="Other",'Hypermutated'="Other")) %>%
  analyse_survival(vars(os, deceased), by=BRAF) -> result

#map(vars(age_int, grade, cIMPACT, c7Gain_10Loss, CDKN2AB_Loss, subtype, Tumor_Location, TERTp_Mut, IDH1.2, MGMT, BLM, WT1, PDGFRB, NFKBIA, PRDM1, FLT4, FLCN, CCND2, DNMT3A, KDR, REL, PTEN, IKZF1)

#ARID1B, TCF3, GATA3, BRD4, IKZF1, EZH2, FAM131B, AGK, ERCC2, SMO
df_forest.adult_gbm %>%
  analyse_survival(vars(os,deceased), by="BRD4")

print(result)
kaplan_meier_plot(result, default_args)

#temp_df <- df_forest.adult[df_forest.adult$IKZF1 == 1,]

list(result,
     df_forest %>%
       analyse_survival(vars(os, deceased), 
                        by=age_int)
     ) %>%
  kaplan_meier_grid(nrow=2,
                    default_args,
                    break.time.by="breakByQuarterYear",
                    mapped_plot_args=list(
                      legend.title=c("BRAF Alteration", "Age Group"),
                      title=c("A", "B")
                    )) %>%
  print

temp_df<- colSums(Filter(is.numeric, df_forest.adult_gbm))
table(df_forest$TP53)
colnames(df.wide.muts[-1])

## Multiple univariate analysis with all
## Features: each BRAF alteration cluster, 1p/19q loss, IDH mut status, H3K27M mut, chr.7gain, chr 10 loss, TERT promoter mutation, EGFR amp, MGMT methylation status, hypermutation status, age, (tumor location, extent of resection, KPS). CDKN2a/b loss, EGFR mut, TP53, NF1, PTEN, etc
#(colnames(df.wide.nmf_alts_input))[3:241]
#map(vars(grade, subtype), function(by)
#map(vars(age_int, Gender, braf_class, IDH1.2, CDKN2AB_Loss, NF1, PTEN, EGFR, ATRX), function(by)
#map(vars(age_num, age_int, Gender, grade, subtype, Tumor_Location, Tumor_Laterality, Extent_Resection, braf_class, cIMPACT, IDH1.2, x1p_19q, H3K27M, MGMT, c7Gain_10Loss, TERTp_Mut, EGFR_Gain, CDKN2AB_Loss), function(by) #basic
#map(vars(age_int, Gender, Tumor_Location, Tumor_Laterality, Extent_Resection, braf_class, cIMPACT, IDH1.2, x1p_19q, H3K27M, MGMT, c7Gain_10Loss, TERTp_Mut, EGFR_Gain, CDKN2AB_Loss), function(by)

#map(vars(age_num, age_int, Gender, Tumor_Location, Tumor_Laterality, Extent_Resection, grade, subtype, braf_class, cIMPACT, IDH1.2, x1p_19q, H3K27M, MGMT, c7Gain_10Loss, TERTp_Mut, EGFR_Gain, CDKN2AB_Loss, 'CDKN2A/B', KIAA1549, TP53, chr7q, ATRX, chr7p, BRCA2, 'PIK3CA/3R1', NF1, chr10q, chr10p, PTEN, CREBBP, APC, ATM, NOTCH1, SETD2, EGFR, 'IDH1/2', TERT_Promoter, MET, SMO, MSH6, PDGFRA, ROS1, ARID1B, EZH2, CNTNAP2, FGFR3, PRKDC, EP300, AGK, FAM131B, TET2, TSC2, chr19q, DICER1, NOTCH2, ARID2, ASXL1, 'CDK4/6', CUX1, KDR, PTCH1, FLT1, UBE2H, ALK, FANCA, MTOR, ABL1, BCOR, BRCA1, DNMT3A, PTPRZ1, ARID1A, KIT, MYC, PALB2, RB1, SF3B1, BRD4, BRIP1, ERBB3, ERBB4, FLT4, TSC1, AXL, BLM, IGF1R, MAP3K1, 'MDM2/4', PMS2, CARD11, ERBB2, FGFR1, RET, STAG2, DDR2, FANCD2, FLT3, SMARCA4, ZNF217, FGFR4, GLI1, NTRK1, PRDM1, PTPN11, TERT, WRN, CCND2, chr1p, CYLD, MLH1, MSH2, PMS1, BAP1, CBL, DIS3, ERCC5, JAK2, MITF, PIK3C2B, TCF7L2, FBXW7, FGFR2, IKZF1, JAK3, KDM6A, MPL, NTRK3, TCF3, AKT1, CDKN1B, CDKN2C, CNTN1, ERCC2, ETV1, EXT1, GLI2, GNAS, NTRK2, PAX5, SUFU, BCORL1, CDH1, FANCG, MLLT3, RAD21, SMARCB1, STAT6, XPO1, BMPR1A, CD274, CHEK2, ERCC3, ESR1, ETV6, FLCN, GNAI1, HNF1A, HRAS, MCL1, MEN1, MYCL, NF2, PARK2, SH2B3, SMAD2, WT1, BCL6, CBLB, CCDC6, GATA3, MAPK1, PDGFRB, SDHC, SETBP1, SMC3, SOX2, SOX9, STK11, TNFAIP3, AKT2, AR, ARAF, AURKA, ERCC4, GNA11, MAP2K4, MUTYH, MYB, MYBL1, MYD88, NFE2L2, NPM1, PHOX2B, RARA, RUNX1, BCL2, C17orf70, CCNE1, CRLF1, CRTC1, 'ERCC6-PGBD3', KRAS, MEF2B, NBN, PDE7A, SUZ12, ARL11, AURKB, CDC73, CNTN3, ETV4, FANCC, FANCE, FAS, FH, NFKBIA, PDCD1LG2, RALA, SMAD4, AKT3, BCL2L12, BUB1B, CCND1, CDKN1A, CEBPA, CRLF3, CTNNB1, DDB2, FAM46C, FANCF, GATA4, GATA6, GIT2, H3K27M, INHBC, MGMT, PHF6, PIM1, PXMP2, REL, SDHAF2, SDHB, SOCS1, STAT3, TACC3, UTP23, B2M, BCL2L1, GPC3, H3F3A, MAP4K5, 'NKX2-1', RAF1, U2AF1, XPA, ZRSR2), function(by) #all

#map(vars(age_num, age_int, Gender, Tumor_Location, Tumor_Laterality, braf_class, cIMPACT, H3K27M, MGMT, c7Gain_10Loss, TERTp_Mut, EGFR_Gain, CDKN2AB_Loss, 'CDKN2A/B', KIAA1549, TP53, chr7q, ATRX, chr7p, BRCA2, 'PIK3CA/3R1', NF1, chr10q, chr10p, PTEN, CREBBP, APC, ATM, NOTCH1, SETD2, EGFR, TERT_Promoter, MET, SMO, MSH6, PDGFRA, ROS1, ARID1B, EZH2, CNTNAP2, FGFR3, PRKDC, EP300, AGK, FAM131B, TET2, TSC2, chr19q, DICER1, NOTCH2, ARID2, ASXL1, 'CDK4/6', CUX1, KDR, PTCH1, FLT1, UBE2H, ALK, FANCA, MTOR, ABL1, BCOR, BRCA1, DNMT3A, PTPRZ1, ARID1A, KIT, MYC, PALB2, RB1, SF3B1, BRD4, BRIP1, ERBB3, ERBB4, AXL, BLM, IGF1R, MAP3K1, 'MDM2/4', PMS2, CARD11, ERBB2, FGFR1, RET, STAG2, DDR2, FANCD2, FLT3, SMARCA4, ZNF217, FGFR4, GLI1, NTRK1, PRDM1, PTPN11, TERT, WRN, CCND2, chr1p, CYLD, MLH1, MSH2, PMS1, BAP1, CBL, DIS3, ERCC5, JAK2, MITF, PIK3C2B, TCF7L2, FBXW7, FGFR2, IKZF1, JAK3, KDM6A, MPL, NTRK3, TCF3, AKT1, CDKN1B, CDKN2C, CNTN1, ERCC2, EXT1, GLI2, GNAS, NTRK2, PAX5, SUFU, BCORL1, FANCG, MLLT3, RAD21, SMARCB1, STAT6, XPO1, BMPR1A, CD274, CHEK2, ERCC3, ESR1, GNAI1, HNF1A, MCL1, MYCL, NF2, PARK2, SH2B3, SMAD2, WT1, BCL6, CBLB, CCDC6, GATA3, MAPK1, PDGFRB, SDHC, SMC3, SOX2, SOX9, STK11, TNFAIP3, AKT2, AR, AURKA, ERCC4, GNA11, MAP2K4, MUTYH, MYB, MYBL1, MYD88, NFE2L2, NPM1, PHOX2B, RARA, RUNX1, BCL2, CCNE1, CRLF1, CRTC1, 'ERCC6-PGBD3', KRAS), function(by) ##GBM only

#map(vars(age_int, grade, cIMPACT, c7Gain_10Loss, CDKN2AB_Loss, subtype, Tumor_Location, TERTp_Mut, IDH1.2, MGMT, BLM, WT1, PDGFRB, NFKBIA, PRDM1, FLT4, FLCN, CCND2, DNMT3A, KDR, PTEN, IKZF1), function(by) # Gender, Tumor_Location, grade, subtype, braf_class, cIMPACT, IDH1.2, MGMT, c7Gain_10Loss, TERTp_Mut, EGFR_Gain, CDKN2AB_Loss, 'CDKN2A/B', KIAA1549, TP53, chr7q, ATRX, chr7p, BRCA2, 'PIK3CA/3R1', NF1, chr10q, chr10p, PTEN, CREBBP, APC, ATM, NOTCH1, SETD2, EGFR, 'IDH1/2', TERT_Promoter, MET, SMO, MSH6, PDGFRA, ROS1, ARID1B, EZH2, CNTNAP2, FGFR3, PRKDC, EP300, AGK, FAM131B, TET2, TSC2, chr19q, DICER1, NOTCH2, ARID2, ASXL1, 'CDK4/6', CUX1, KDR, PTCH1, FLT1, UBE2H, ALK, FANCA, MTOR, ABL1, BCOR, BRCA1, DNMT3A, PTPRZ1, ARID1A, KIT, MYC, PALB2, RB1, SF3B1, BRD4, BRIP1, ERBB3, ERBB4, FLT4, TSC1, AXL, BLM, IGF1R, MAP3K1, 'MDM2/4', PMS2, CARD11, ERBB2, FGFR1, RET, STAG2, DDR2, FANCD2, FLT3, SMARCA4, ZNF217, FGFR4, GLI1, NTRK1, PRDM1, PTPN11, TERT, WRN, CCND2, chr1p, CYLD, MLH1, MSH2, PMS1, BAP1, CBL, DIS3, ERCC5, JAK2, MITF, PIK3C2B, TCF7L2, FBXW7, FGFR2, IKZF1, JAK3, KDM6A, MPL, NTRK3, TCF3, AKT1, CDKN1B, CDKN2C, CNTN1, ERCC2, ETV1, EXT1, GLI2, GNAS, NTRK2, PAX5, SUFU, BCORL1, CDH1, FANCG, MLLT3, RAD21, SMARCB1, STAT6, XPO1, BMPR1A, CD274, CHEK2, ERCC3, ESR1, ETV6, FLCN, GNAI1, HNF1A, HRAS, MCL1, MEN1, MYCL, NF2, PARK2, SH2B3, SMAD2, WT1, BCL6, CBLB, CCDC6, GATA3, MAPK1, PDGFRB, SDHC, SETBP1, SMC3, SOX2, SOX9, STK11, TNFAIP3, AKT2, AR, ARAF, AURKA, ERCC4, GNA11, MAP2K4, MUTYH, MYB, MYBL1, MYD88, NFE2L2, NPM1, PHOX2B, RARA, RUNX1, BCL2, C17orf70, CCNE1, CRLF1, CRTC1, 'ERCC6-PGBD3', KRAS, MEF2B, NBN, PDE7A, SUZ12, ARL11, AURKB, CDC73, CNTN3, ETV4, FANCC, FANCE, FAS, FH, NFKBIA, PDCD1LG2, RALA, SMAD4, AKT3, BCL2L12, BUB1B, CCND1, CDKN1A, CEBPA, CRLF3, CTNNB1, DDB2, FAM46C, FANCF, GATA4, GATA6, GIT2, H3K27M, INHBC, MGMT, PHF6, PIM1, PXMP2, REL, SDHAF2, SDHB, SOCS1, STAT3, TACC3, UTP23, B2M, BCL2L1, GPC3, H3F3A, MAP4K5, 'NKX2-1', RAF1, U2AF1, XPA, ZRSR2), function(by)
map(vars(age_int, subtype, grade, braf_class, cIMPACT, MGMT, IDH1.2, TERTp_Mut, PTEN, CDKN2AB_Loss,  c7Gain_10Loss)
#map(vars(age_int, grade, cIMPACT, c7Gain_10Loss, CDKN2AB_Loss, subtype, Tumor_Location, TERTp_Mut, IDH1.2, MGMT, BLM, WT1, PDGFRB, NFKBIA, PRDM1, FLT4, FLCN, CCND2, DNMT3A, KDR, REL, PTEN, IKZF1) #Adult signif
#map(vars(age_int, ARID1B, TCF3, GATA3, BRD4, IKZF1, EZH2, FAM131B, AGK, SMO), function(by) #GBM signif
   {
  analyse_multivariate(df_forest.adult_gbm,
                       vars(os, deceased),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names,
                       reference_level_dict = c(age_int="35-50", subtype="Astro, IDH-mt", braf_class="Other"))
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10))

## ***Multiple univariate analysis on all adults with updated features
(colnames(df_forest[,52:ncol(df_forest)]))
(colSums(df_forest[,52:ncol(df_forest)]))
df_forest$age_int
df.wide.nmf_alts_input$V1
#map(vars(age_int, subtype, grade, braf_class, cIMPACT, MGMT, IDH1.2, TERTp_Mut, PTEN, CDKN2AB_Loss, c7Gain_10Loss, NF1, ATRX), function(by)
#map(vars(age_int, subtype, grade, braf_class, cIMPACT, MGMT, IDH1.2, TERTp_Mut, PTEN, CDKN2AB_Loss, c7Gain_10Loss, NF1, ATRX, "BRAF", "CDKN2A/B", "KIAA1549", "TP53", "chr7q", "ATRX", "chr7p", "BRCA2", "PIK3CA/3R1", "NF1", "chr10q", "chr10p", "PTEN", "CREBBP", "APC", "ATM", "NOTCH1", "SETD2", "EGFR", "IDH1/2", "TERT_Promoter", "MET", "SMO", "MSH6", "PDGFRA", "ROS1", "ARID1B", "EZH2", "CNTNAP2", "FGFR3", "PRKDC", "EP300", "AGK", "FAM131B", "TET2", "TSC2", "chr19q", "DICER1", "NOTCH2", "ARID2", "ASXL1", "CDK4/6", "CUX1", "KDR", "PTCH1", "FLT1", "UBE2H", "ALK", "FANCA", "MTOR", "ABL1", "BCOR", "BRCA1", "DNMT3A", "PTPRZ1", "ARID1A", "KIT", "MYC", "PALB2", "RB1", "SF3B1", "BRD4", "BRIP1", "ERBB3", "ERBB4", "FLT4", "TSC1", "AXL", "BLM", "IGF1R", "MAP3K1", "MDM2/4", "PMS2", "CARD11", "ERBB2", "FGFR1", "RET", "STAG2", "DDR2", "FANCD2", "FLT3", "SMARCA4", "ZNF217", "FGFR4", "GLI1", "NTRK1", "PRDM1", "PTPN11", "TERT", "WRN", "CCND2", "chr1p", "CYLD", "MLH1", "MSH2", "PMS1", "BAP1", "CBL", "DIS3", "ERCC5", "JAK2", "MITF", "PIK3C2B", "TCF7L2", "FBXW7", "FGFR2", "IKZF1", "JAK3", "KDM6A", "MPL", "NTRK3", "TCF3", "AKT1", "CDKN1B", "CDKN2C", "CNTN1", "ERCC2", "ETV1", "EXT1", "GLI2", "GNAS", "NTRK2", "PAX5", "SUFU", "BCORL1", "CDH1", "FANCG", "MLLT3", "RAD21", "SMARCB1", "STAT6", "XPO1", "BMPR1A", "CD274", "CHEK2", "ERCC3", "ESR1", "ETV6", "FLCN", "GNAI1", "HNF1A", "HRAS", "MCL1", "MEN1", "MYCL", "NF2", "PARK2", "SH2B3", "SMAD2", "WT1", "BCL6", "CBLB", "GATA3", "MAPK1", "PDGFRB", "SDHC", "SETBP1", "SOX2", "SOX9", "STK11", "TNFAIP3", "AKT2", "AR", "ARAF", "AURKA", "ERCC4", "GNA11", "MAP2K4", "MUTYH", "MYBL1", "NFE2L2", "RUNX1", "BCL2", "C17orf70", "CCNE1", "CRLF1", "CRTC1", "KRAS", "SUZ12"), function(by)
#map(vars("BRAF", "CDKN2A/B", "KIAA1549", "TP53", "chr7q", "ATRX", "chr7p", "BRCA2", "PIK3CA/3R1", "NF1", "chr10q", "chr10p", "PTEN", "CREBBP", "APC", "ATM", "NOTCH1", "SETD2", "EGFR", "IDH1/2", "TERT_Promoter", "MET", "SMO", "MSH6", "PDGFRA", "ROS1", "ARID1B", "EZH2", "CNTNAP2", "FGFR3", "PRKDC", "EP300", "AGK", "FAM131B", "TET2", "TSC2", "chr19q", "DICER1", "NOTCH2", "ARID2", "ASXL1", "CDK4/6", "CUX1", "KDR", "PTCH1", "FLT1", "UBE2H", "ALK", "FANCA", "MTOR", "ABL1", "BCOR", "BRCA1", "DNMT3A", "PTPRZ1", "ARID1A", "KIT", "MYC", "PALB2", "RB1", "SF3B1", "BRD4", "BRIP1", "ERBB3", "ERBB4", "FLT4", "TSC1", "AXL", "BLM", "IGF1R", "MAP3K1", "MDM2/4", "PMS2", "CARD11", "ERBB2", "FGFR1", "RET", "STAG2", "DDR2", "FANCD2", "FLT3", "SMARCA4", "ZNF217", "FGFR4", "GLI1", "NTRK1", "PRDM1", "PTPN11", "TERT", "WRN", "CCND2", "chr1p", "CYLD", "MLH1", "MSH2", "PMS1", "BAP1", "CBL", "DIS3", "ERCC5", "JAK2", "MITF", "PIK3C2B", "TCF7L2", "FBXW7", "FGFR2", "IKZF1", "JAK3", "KDM6A", "MPL", "NTRK3", "TCF3", "AKT1", "CDKN1B", "CDKN2C", "CNTN1", "ERCC2", "ETV1", "EXT1", "GLI2", "GNAS", "NTRK2", "PAX5", "SUFU", "BCORL1", "CDH1", "FANCG", "MLLT3", "RAD21", "SMARCB1", "STAT6", "XPO1", "BMPR1A", "CD274", "CHEK2", "ERCC3", "ESR1", "ETV6", "FLCN", "GNAI1", "HNF1A", "HRAS", "MCL1", "MEN1", "MYCL", "NF2", "PARK2", "SH2B3", "SMAD2", "WT1", "BCL6", "CBLB", "GATA3", "MAPK1", "PDGFRB", "SDHC", "SETBP1", "SOX2", "SOX9", "STK11", "TNFAIP3", "AKT2", "AR", "ARAF", "AURKA", "ERCC4", "GNA11", "MAP2K4", "MUTYH", "MYBL1", "NFE2L2", "RUNX1", "BCL2", "C17orf70", "CCNE1", "CRLF1", "CRTC1", "KRAS", "SUZ12"), function(by)

df_forest.adult_gbm %>%
  analyse_survival(vars(os,deceased), by="BRD4")

var_list.adult = vars("IDH1.2", "CDKN2AB_Loss", "c7Gain_10Loss", "NF1", "ATRX", "PTEN", "IKZF1", "TERT_Promoter", "PDGFRB", "PRDM1")

df_forest.adult %>%
  analyse_survival(vars(os,deceased), by=(var_list.adult), p_adjust_method = "BH")
is.numeric(df_forest$age_num)
## Starting univariate analysis here
map(vars(age_int, subtype, grade, braf_class), function(by) #clinicopathologic
#map(vars(c7Gain_10Loss, "CDKN2A/B", "IDH1.2", "PTEN", "TERT_Promoter", "EGFR", "NF1", "ATRX", "IKZF1", "PRDM1", "DNMT3A", "PARK2", "CCND2", "KIAA1549", "CREBBP", "AGK", "STAG2", "MET"), function(by) # just the genomic data for all, removed CDKN2AB_Loss and PDGFRB
#map(vars(IDH1.2, CDKN2AB_Loss, c7Gain_10Loss, NF1, ATRX, "PTEN", "IKZF1", "TERT_Promoter", "PDGFRB", "PRDM1"), function(by) # just the genomic data for adults only
#map(vars(age_int, subtype, grade, braf_class, cIMPACT, MGMT, IDH1.2, TERTp_Mut, PTEN, CDKN2AB_Loss, c7Gain_10Loss, NF1, ATRX, "CDKN2A/B", "PTEN", "TERT_Promoter", "IKZF1", "PRDM1", "DNMT3A", "PARK2", "CCND2", "KIAA1549", "PDGFRB", "CREBBP", "AGK", "STAG2", "MET"), function(by) # all
#map(vars(age_int, subtype, grade, braf_class, cIMPACT, MGMT, IDH1.2, TERTp_Mut, PTEN, CDKN2AB_Loss, c7Gain_10Loss, NF1, ATRX, "CDKN2A/B", "PTEN", "IKZF1", "TERT_Promoter", "PDGFRB", "PRDM1"), function(by) #adults only
{
  #analyse_multivariate(gbm_os_data,
  #analyse_multivariate(df_forest.adult,
  analyse_multivariate(df_forest,
                       vars(os, deceased),
                       covariates = list(by),
                       #p_adjust_method = "none",
                       covariate_name_dict = covariate_names,
                       reference_level_dict = c(age_int="35-50", braf_class="Other", subtype="HGG"))
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS (months)"),
              #orderer = ~order(HR),
              #orderer = ~order(p),
              labels_displayed = c("endpoint", "factor"),
              values_displayed = c("HR", "p", "n"),
              value_headers = c(HR = "HR", p = "p", n = "n"),
              HR_x_limits = c(0.02, 35),
              HR_x_breaks = c(0.03, 0.06, 0.12, 0.25, 0.5, 1, 2, 4, 8, 16, 32),
              ggtheme = ggplot2::theme_bw(base_size = 10))

temp_df <- df.wide.muts[,c("CDKN2A/B", "PTEN", "TERT_Promoter", "IKZF1", "PRDM1", "DNMT3A", "PARK2", "CCND2", "KIAA1549", "PDGFRB", "CREBBP", "AGK", "STAG2", "MET")]
table(df_forest$`CDKN2A/B`)
table(df_forest.adult$PDGFRB)
table(df.wide.nmf_alts_input$PDGFRB)

## Check surve curve
ggsurvplot(
  fit = survfit(Surv(os, deceased) ~DNMT3A, data = df_forest),
  #title = "Adults with BRAF-altered Astrocytoma, IDH-mt",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("gray50", "red")
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)

#c("CDKN2A/B", "PTEN", "TERT_Promoter", "IKZF1", "PRDM1", "DNMT3A", "PARK2", "CCND2", "KIAA1549", "PDGFRB", "CREBBP", "AGK", "STAG2", "MET")[!(c("CDKN2A/B", "PTEN", "TERT_Promoter", "IKZF1", "PRDM1", "DNMT3A", "PARK2", "CCND2", "KIAA1549", "PDGFRB", "CREBBP", "AGK", "STAG2", "MET") %in% IntersectGeneList)]

colSums(temp_df)
df.wide.mutations$mutations.ep
## Multiple univariate analysis on Adult GBM only
df_forest$membership[df_forest$subtype=="GBM, IDH-wt"]
table(df_forest$subtype)
map(vars(age_int, gender, braf_class, IDH_canonical, CDKN2A.B, NF1, PTEN, EGFR, ATRX), function(by)
{
  #analyse_multivariate(gbm_os_data,
  analyse_multivariate(df_forest[df_forest$subtype=="GBM, IDH-wt",],
                       vars(os, deceased),
                       covariates = list(by), # covariates expects a list
                       covariate_name_dict = covariate_names)
}) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS (months)"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10))

table(df_forest$WHO_Diagnosis[df_forest$subtype == "LGG"])
df_forest$IDH1.2[df_forest$WHO_Diagnosis == "Oligodendroglioma, NEC (1p/19q intact)"]

table(df_forest.adult$braf_class == "Class I")

## Multivariate analysis
table(df_forest.adult$TERTp_Mut)
res.cox <- coxph(Surv(os, deceased) ~ age_int + grade + cIMPACT + c7Gain_10Loss + CDKN2AB_Loss + subtype + Tumor_Location + TERTp_Mut + IDH1.2 + MGMT + BLM + WT1 + PDGFRB + NFKBIA + PRDM1 + FLT4 + FLCN + CCND2 + DNMT3A + KDR + PTEN + IKZF1, data= df_forest.adult_gbm)
res.cox <- coxph(Surv(os, deceased) ~ age_int + braf_class ARID1B + TCF3 + GATA3 + BRD4 + IKZF1 + EZH2 + FAM131B + AGK + SMO, data = df_forest.adult_gbm)
summary(res.cox)
ggforest(res.cox, fontsize=1)

df_forest.adult %>%
  analyse_multivariate(vars(os, deceased),
                       covariates = vars(age_int, grade, cIMPACT, c7Gain_10Loss, CDKN2AB_Loss, subtype, Tumor_Location, TERTp_Mut, IDH1.2, MGMT, BLM, WT1, PDGFRB, NFKBIA, PRDM1, FLT4, FLCN, CCND2, DNMT3A, KDR, REL, PTEN, IKZF1),
                       covariate_name_dict = covariate_names) %>%
  forest_plot(factor_labeller = covariate_names,
              endpoint_labeller = c(time="OS (months)"),
              orderer = ~order(HR),
              labels_displayed = c("endpoint", "factor", "n"),
              ggtheme = ggplot2::theme_bw(base_size = 10))

## Testing out with combined GBM only multivariate
## Covariate names
covariate_names <- c(age_num="Age", age_int="Age Group", gender="Gender", Gender="Gender", subtype_path="Histology", subtype="Pathology", grade="Grade", braf_class="BRAF Mutation Class", cIMPACT="cIMPACT3 Defined GBM", IDH1.2="IDH1/2", IDH_canonical="Canonical IDH-mt", CDKN2A.B="CDKN2A/B Loss",CDKN2AB_Loss="CDKN2A/B Loss", NF1="NF1", PTEN="PTEN", EGFR="EGFR", ATRX="ATRX")

gbm_os_data %>%
  filter(!is.na(os_days)) %>%
  mutate(age_int = as.factor(age_int), braf_class = as.factor(braf_class)) %>%
  analyse_multivariate(vars(os_days, death),
                       covariates = vars(age_num, age_int, gender, braf_class),
                       covariate_name_dict = covariate_names,
                       reference_level_dict = c(age_int="Young Adult", braf_class="Other")) ->
  result
print(result)
forest_plot(result)

### End HR Forest Plot ###

temp_df <- df_forest[df_forest$braf_class == "Class II",]
temp_df$sample[temp_df$EGFR == 1]
temp_df <- df_forest[df_forest$PRDM1 == 1,]


temp_df <- df_forest.adult[df_forest.adult$braf_class == "Class I",]
temp_df <- nmf_alts[nmf_alts$BEST_EFF_GENE == "CDKN2A/B" & nmf_alts$SAMPLE_ACCESSION_NBR %in% adult_classI_samps,]
temp_df <- unique(temp_df)
table(temp_df$variant_classification)
