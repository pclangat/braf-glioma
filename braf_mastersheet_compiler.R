##### Clinical Data Table Parsing #####
## This is for parsing data for master sheet for all BRAF glioma cases (adult, peds, all cohorts)
## Compiling mastersheet df for file to write

## Prepare JHH df
jhh_os_data$SAMPLE_ACCESSION_NBR <- paste("JHH-", jhh_os_data$redcap_patient_id, sep = "")
jhh_os_data$PATIENT_ID <- jhh_os_data$redcap_patient_id

jhh_os_data$age_int2[jhh_os_data$age_at_dx <18] <- "Pediatric"
jhh_os_data$age_int2[jhh_os_data$age_at_dx >= 18 & jhh_os_data$age_at_dx < 35] <- "Young Adult"
jhh_os_data$age_int2[jhh_os_data$age_at_dx > 34] <- "Middle"
jhh_os_data$age_int2[jhh_os_data$age_at_dx > 50] <- "Old"
jhh_os_data$Age_Interval <- jhh_os_data$age_int

jhh_os_data$Age_Interval[jhh_os_data$Age_Interval== 'Pediatric'] <- "<18"
jhh_os_data$Age_Interval[jhh_os_data$Age_Interval == 'Young Adult'] <- "18-34"
jhh_os_data$Age_Interval[jhh_os_data$Age_Interval == 'Middle'] <- "35-50"
jhh_os_data$Age_Interval[jhh_os_data$Age_Interval == 'Old'] <- ">50"

### Adult BRAF Mastersheet ###
adult.braf.mastersheet <- classifier[!(classifier$age_int == "<18") & !(grepl("JHH", classifier$SAMPLE_ACCESSION_NBR)),c("SAMPLE_ACCESSION_NBR", "PATIENT_ID")]
adult.braf.mastersheet <- rbind(adult.braf.mastersheet, (jhh_os_data[jhh_os_data$age_at_dx >18.0, c("SAMPLE_ACCESSION_NBR", "PATIENT_ID")]))

adult.braf.mastersheet$PATIENT_ID_19 <- master.sheet.dfci$pt_id_19[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]

adult.braf.mastersheet$Cohort <- classifier$cohort[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$Cohort[grepl("JHH", adult.braf.mastersheet$SAMPLE_ACCESSION_NBR)] <- "JHH"
                              
adult.braf.mastersheet$Panel_Version <- classifier$PANEL_VERSION[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$Panel_Version <- ifelse(is.na(
  adult.braf.mastersheet$Panel_Version), jhh_panel_data$PANEL_VERSION[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_panel_data$SAMPLE_ACCESSION_NBR
  )], adult.braf.mastersheet$Panel_Version)

adult.braf.mastersheet$Gender <- classifier$Gender[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$Gender <- ifelse(is.na(
  adult.braf.mastersheet$Gender), jhh_os_data$patient_gender[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR
  )], adult.braf.mastersheet$Gender)

adult.braf.mastersheet$Age <- classifier$age_num[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$Age <- ifelse(is.na(
  adult.braf.mastersheet$Age), jhh_os_data$age_at_dx[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$Age)

adult.braf.mastersheet$Age_Interval <- classifier$age_int[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$Age_Interval <- as.character(adult.braf.mastersheet$Age_Interval)
adult.braf.mastersheet$Age_Interval <- ifelse(is.na(
  adult.braf.mastersheet$Age_Interval), jhh_os_data$Age_Interval[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$Age_Interval)

adult.braf.mastersheet$MRN <- master.sheet.dfci$MRN[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]

adult.braf.mastersheet$MRN_BCH <- master.sheet.dfci$CHILDRENS_MRN[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$MRN_BCH[adult.braf.mastersheet$MRN_BCH == 'null'] <- NA 

adult.braf.mastersheet$DOB <- master.sheet.dfci$DOB[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]

adult.braf.mastersheet$DODiagnosis <- adult.outcomes.verified$Date.of.Radiology.Diagnosis[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number)]

adult.braf.mastersheet$DOSurgery <- ifelse(!(is.na(
  adult.braf.mastersheet$PATIENT_ID_19)), master.sheet.dfci$DoS[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR
  )], NA)

adult.braf.mastersheet$DOLastFU <- master.sheet.dfci$DERIVED_LAST_CONTACT_DT[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$DOLastFU[grepl("-", adult.braf.mastersheet$DOLastFU)] <- format(strptime(as.character(
  adult.braf.mastersheet$DOLastFU[grepl("-", adult.braf.mastersheet$DOLastFU)]), "%d-%b-%y"), "%Y-%m-%d")
adult.braf.mastersheet$DOLastFU[grepl("/",adult.braf.mastersheet$DOLastFU)] <- format(strptime(as.character(
  adult.braf.mastersheet$DOLastFU[grepl("/", adult.braf.mastersheet$DOLastFU)]), "%m/%d/%Y"), "%Y-%m-%d")
# correct one sample based on second surgery final report
adult.braf.mastersheet$DOLastFU[adult.braf.mastersheet$SAMPLE_ACCESSION_NBR== "BL-16-W04153"] <- "2015-12-07"
adult.braf.mastersheet$DOLastFU[adult.braf.mastersheet$DOLastFU=='null'] <- NA

adult.braf.mastersheet$DOD <- adult.outcomes.verified$Date.of.Death[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number
)]

adult.braf.mastersheet$Alive <- adult.outcomes.verified$Alive..1.0.[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number)]
adult.braf.mastersheet$Alive <- ifelse(is.na(
  adult.braf.mastersheet$Alive), jhh_os_data$alive[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR
  )], adult.braf.mastersheet$Alive)

#adult.braf.mastersheet$Death <- classifier$death[match(
#  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
#adult.braf.mastersheet$Death <- ifelse(is.na(adult.braf.mastersheet$Death), (jhh_os_data$death[match(
#  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), adult.braf.mastersheet$Death)
#adult.braf.mastersheet$Death <- as.numeric(adult.braf.mastersheet$Death)

## OS 
adult.braf.mastersheet$OS <- (as.Date(as.character(adult.braf.mastersheet$DOD), format="%Y-%m-%d")-as.Date(as.character(adult.braf.mastersheet$DODiagnosis), format="%Y-%m-%d"))/30.4167

adult.braf.mastersheet$OS <- ifelse(is.na(adult.braf.mastersheet$OS), (as.Date(as.character(
  adult.braf.mastersheet$DOD), format="%Y-%m-%d")- as.Date(as.character(
    adult.braf.mastersheet$DOSurgery), format="%Y-%m-%d"))/30.4167, adult.braf.mastersheet$OS)

adult.braf.mastersheet$OS <- ifelse(is.na(adult.braf.mastersheet$OS) & adult.braf.mastersheet$Alive == 1, (as.Date(as.character(adult.braf.mastersheet$DOLastFU), format="%Y-%m-%d")- as.Date(as.character(
    adult.braf.mastersheet$DODiagnosis), format="%Y-%m-%d"))/30.4167, adult.braf.mastersheet$OS)

adult.braf.mastersheet$OS <- ifelse(is.na(adult.braf.mastersheet$OS) & adult.braf.mastersheet$Alive == 1, (as.Date(as.character(adult.braf.mastersheet$DOLastFU), format="%Y-%m-%d")- as.Date(as.character(
  adult.braf.mastersheet$DOSurgery), format="%Y-%m-%d"))/30.4167, adult.braf.mastersheet$OS)

adult.braf.mastersheet$OS <- ifelse(is.na(adult.braf.mastersheet$OS), (tcga.lgg.outcomes$OS.time[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.lgg.outcomes$bcr_patient_barcode
)])/30.4167, adult.braf.mastersheet$OS)

adult.braf.mastersheet$OS <- ifelse(is.na(adult.braf.mastersheet$OS), (tcga.gbm.outcomes$OS..days.[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.gbm.outcomes$Case.ID
)])/30.4167, adult.braf.mastersheet$OS)

adult.braf.mastersheet$OS <- ifelse(is.na(adult.braf.mastersheet$OS), as.double((
  genie.outcomes$Overall.Survival..Months.[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, genie.outcomes$SAMPLE_ACCESSION_NBR)])), adult.braf.mastersheet$OS)

adult.braf.mastersheet$OS <- ifelse(is.na(adult.braf.mastersheet$OS), (jhh_os_data$overall_survival[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), adult.braf.mastersheet$OS)

adult.braf.mastersheet$OS <- ifelse((adult.braf.mastersheet$Alive == 1 & is.na(adult.braf.mastersheet$OS)), 0, adult.braf.mastersheet$OS)

adult.braf.mastersheet$OS <- round(adult.braf.mastersheet$OS, 1)

## OS days
adult.braf.mastersheet$OS_days <- (as.Date(as.character(adult.braf.mastersheet$DOD), format="%Y-%m-%d")-
  as.Date(as.character(adult.braf.mastersheet$DODiagnosis), format="%Y-%m-%d"))
adult.braf.mastersheet$OS_days <- ifelse(is.na(adult.braf.mastersheet$OS_days), as.Date(as.character(
  adult.braf.mastersheet$DOD), format="%Y-%m-%d")- as.Date(as.character(
    adult.braf.mastersheet$DOSurgery), format="%Y-%m-%d"), adult.braf.mastersheet$OS_days)
adult.braf.mastersheet$OS_days <- ifelse(is.na(adult.braf.mastersheet$OS_days) & adult.braf.mastersheet$Alive == 1, (as.Date(as.character(adult.braf.mastersheet$DOLastFU), format="%Y-%m-%d")- as.Date(as.character(
  adult.braf.mastersheet$DODiagnosis), format="%Y-%m-%d")), adult.braf.mastersheet$OS_days)
adult.braf.mastersheet$OS_days <- ifelse(is.na(adult.braf.mastersheet$OS_days) & adult.braf.mastersheet$Alive == 1, (as.Date(as.character(adult.braf.mastersheet$DOLastFU), format="%Y-%m-%d")- as.Date(as.character(
  adult.braf.mastersheet$DOSurgery), format="%Y-%m-%d")), adult.braf.mastersheet$OS_days)
adult.braf.mastersheet$OS_days <- ifelse(is.na(adult.braf.mastersheet$OS_days), tcga.lgg.outcomes$OS.time[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.lgg.outcomes$bcr_patient_barcode
)], adult.braf.mastersheet$OS_days)
adult.braf.mastersheet$OS_days <- ifelse(is.na(adult.braf.mastersheet$OS_days), tcga.gbm.outcomes$OS..days.[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.gbm.outcomes$Case.ID
)], adult.braf.mastersheet$OS_days)
adult.braf.mastersheet$OS_days <- ifelse(is.na(adult.braf.mastersheet$OS_days), as.integer((
  genie.outcomes$Overall.Survival..Months.[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, genie.outcomes$SAMPLE_ACCESSION_NBR)])*30.4167), adult.braf.mastersheet$OS_days)
adult.braf.mastersheet$OS_days <- ifelse(is.na(adult.braf.mastersheet$OS_days), (jhh_os_data$os_days[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), adult.braf.mastersheet$OS_days)
adult.braf.mastersheet$OS_days <- round(adult.braf.mastersheet$OS_days, 0)

## KPS and other clinical
adult.braf.mastersheet$KPS <- adult.outcomes.verified$KPS.at.near.Diagnosis[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number)]

adult.braf.mastersheet$KPS <- ifelse(is.na(adult.braf.mastersheet$KPS), jhh_os_data$KPS[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$KPS)

adult.braf.mastersheet$KPS_Date <- adult.outcomes.verified$Date.of.KPS.Record[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number)]
adult.braf.mastersheet$KPS_Date[adult.braf.mastersheet$KPS_Date == ""] <- NA

adult.braf.mastersheet$Primary <- classifier$Primary[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$Primary[adult.braf.mastersheet$Primary == "TRUE"] <- 1
adult.braf.mastersheet$Primary[adult.braf.mastersheet$Primary == "FALSE"] <- 0
adult.braf.mastersheet$Primary <- ifelse(is.na(adult.braf.mastersheet$Primary), jhh_os_data$Primary[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$Primary)

adult.braf.mastersheet$Tumor_Location <- adult.outcomes.verified$Tumor.Location[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number)]
adult.braf.mastersheet$Tumor_Location <- ifelse(is.na(adult.braf.mastersheet$Tumor_Location), jhh_os_data$Tumor_Location[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$Tumor_Location)
adult.braf.mastersheet$Tumor_Location[adult.braf.mastersheet$Tumor_Location == ""] <- NA

adult.braf.mastersheet$Tumor_Laterality <- adult.outcomes.verified$Tumor.Laterality[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number)]
adult.braf.mastersheet$Tumor_Laterality <- ifelse(is.na(adult.braf.mastersheet$Tumor_Laterality), jhh_os_data$Tumor_Laterality[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$Tumor_Laterality)
adult.braf.mastersheet$Tumor_Laterality[adult.braf.mastersheet$Tumor_Laterality == "Midline"] <- "M"
adult.braf.mastersheet$Tumor_Laterality[adult.braf.mastersheet$Tumor_Laterality == "MIDLINE"] <- "M"
adult.braf.mastersheet$Tumor_Laterality[adult.braf.mastersheet$Tumor_Laterality == "BILATERAL"] <- "R;L"
adult.braf.mastersheet$Tumor_Laterality[adult.braf.mastersheet$Tumor_Laterality == ""] <- NA

adult.braf.mastersheet$Extent_Resection <- adult.outcomes.verified$Extent.of.Resection[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number)]
adult.braf.mastersheet$Extent_Resection <- ifelse(is.na(adult.braf.mastersheet$Extent_Resection), jhh_os_data$Extent_Resection[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$Extent_Resection)
adult.braf.mastersheet$Extent_Resection[adult.braf.mastersheet$Extent_Resection == ""] <- NA
adult.braf.mastersheet$Extent_Resection[adult.braf.mastersheet$Extent_Resection == "BR"] <- 4
adult.braf.mastersheet$Extent_Resection[adult.braf.mastersheet$Extent_Resection == "Biopsy"] <- 1
adult.braf.mastersheet$Extent_Resection[adult.braf.mastersheet$Extent_Resection == "STR"] <- 2
adult.braf.mastersheet$Extent_Resection[adult.braf.mastersheet$Extent_Resection == "GTR"] <- 3

adult.braf.mastersheet$BRAF_Inhibitor <- adult.outcomes.verified$BRAF.inhibitor.ever.[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, adult.outcomes.verified$Sample.Accession.Number)]
adult.braf.mastersheet$BRAF_Inhibitor <- ifelse(is.na(adult.braf.mastersheet$BRAF_Inhibitor), jhh_os_data$BRAF_Inhibitor[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$BRAF_Inhibitor)
adult.braf.mastersheet$BRAF_Inhibitor[adult.braf.mastersheet$BRAF_Inhibitor == ""] <- NA

adult.braf.mastersheet$MGMT <- master.sheet.dfci$MGMT_Mehdi[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$MGMT <- ifelse(is.na(
  adult.braf.mastersheet$MGMT), tcga.gbm.outcomes$MGMT.Status[match(
    adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.gbm.outcomes$Case.ID)], adult.braf.mastersheet$MGMT)
adult.braf.mastersheet$MGMT[adult.braf.mastersheet$MGMT == "UNMETHYLATED"] <- 0
adult.braf.mastersheet$MGMT[adult.braf.mastersheet$MGMT == "METHYLATED"] <- 1
adult.braf.mastersheet$MGMT[adult.braf.mastersheet$MGMT == "PARTIAL"] <- 0.5
adult.braf.mastersheet$MGMT[!(adult.braf.mastersheet$MGMT == "0") & !(adult.braf.mastersheet$MGMT == "1") & !(adult.braf.mastersheet$MGMT == "0.5")] <- 2
adult.braf.mastersheet$MGMT[is.na(adult.braf.mastersheet$MGMT)] <- 2
adult.braf.mastersheet$MGMT[adult.braf.mastersheet$MGMT == "2"] <- NA

adult.braf.mastersheet$IDH <- classifier$IDH1.2[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$IDH[adult.braf.mastersheet$IDH == "TRUE"] <- 1
adult.braf.mastersheet$IDH[adult.braf.mastersheet$IDH == "FALSE"] <- 0

adult.braf.mastersheet$IDH <- ifelse(is.na(adult.braf.mastersheet$IDH), IDH_muts_df$MUT_TYPE[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, IDH_muts_df$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$IDH)
adult.braf.mastersheet$IDH[adult.braf.mastersheet$IDH == "canonical"] <- 1
adult.braf.mastersheet$IDH[adult.braf.mastersheet$IDH == "other"] <- 1

adult.braf.mastersheet$IDH <- ifelse(is.na(adult.braf.mastersheet$IDH), jhh_os_data$IDH[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$IDH)

#Ki67?
#ATRx?
#p53?
#GFAP?
adult.braf.mastersheet$x1p_19q <- classifier$x1p_19q[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$x1p_19q[adult.braf.mastersheet$x1p_19q == "FALSE"] <- 0
adult.braf.mastersheet$x1p_19q[adult.braf.mastersheet$x1p_19q == "TRUE"] <- 1
adult.braf.mastersheet$x1p_19q <- ifelse(is.na(adult.braf.mastersheet$x1p_19q), jhh_os_data$x1p_19q[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$x1p_19q)

#adult.braf.mastersheet$x1p_19q_gen <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% c1p19qloss, 1, NA)

adult.braf.mastersheet$H3K27M <- classifier$H3K27M[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$H3K27M[adult.braf.mastersheet$H3K27M == "FALSE"] <- 0
adult.braf.mastersheet$H3K27M[adult.braf.mastersheet$H3K27M == "TRUE"] <- 1
adult.braf.mastersheet$H3K27M <- ifelse(is.na(adult.braf.mastersheet$H3K27M), jhh_os_data$H3K27M[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], adult.braf.mastersheet$H3K27M)

adult.braf.mastersheet$WHO_Grade <- classifier$grade[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$WHO_Grade <- ifelse(is.na(adult.braf.mastersheet$WHO_Grade), (jhh_os_data$who_grade[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), adult.braf.mastersheet$WHO_Grade)

adult.braf.mastersheet$WHO_Diagnosis <- classifier$WHO_Diagnosis[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$WHO_Diagnosis <- ifelse(is.na(adult.braf.mastersheet$WHO_Diagnosis), (jhh_os_data$who_dx[match(
  adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), adult.braf.mastersheet$WHO_Diagnosis)

#adult.braf.mastersheet$Subtype_Path <- path_subtypes_who$Subtype_Path[match(adult.braf.mastersheet$WHO_Diagnosis, path_subtypes_who$WHO_Diagnosis)]
adult.braf.mastersheet$Subtype_Path <- path_subtypes_who$Subtype[match(adult.braf.mastersheet$WHO_Diagnosis, path_subtypes_who$WHO_Diagnosis)]
path_subtypes_who$Subtype
adult.braf.mastersheet$WHO_Diagnosis[adult.braf.mastersheet$SAMPLE_ACCESSION_NBR == "JHH-52"]

adult.braf.mastersheet$Subtype <- classifier$subtype[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$Subtype <- ifelse(is.na(adult.braf.mastersheet$Subtype), adult.braf.mastersheet$Subtype_Path, adult.braf.mastersheet$Subtype)

adult.braf.mastersheet$c7Gain_10Loss <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% classifier$SAMPLE_ACCESSION_NBR, 0, NA)
adult.braf.mastersheet$EGFR_Gain <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% classifier$SAMPLE_ACCESSION_NBR, 0, NA)
adult.braf.mastersheet$TERTP_Mut <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% classifier$SAMPLE_ACCESSION_NBR, 0, NA)
adult.braf.mastersheet$c7Gain_10Loss <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% c7gain10loss, 1, adult.braf.mastersheet$c7Gain_10Loss)
adult.braf.mastersheet$EGFR_Gain <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% egfrgain, 1, adult.braf.mastersheet$EGFR_Gain)
adult.braf.mastersheet$TERTP_Mut <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% tertpmut, 1, adult.braf.mastersheet$TERTP_Mut)

adult.braf.mastersheet$cIMPACT <- classifier$cIMPACT[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]

## Additional gene mutations

adult.braf.mastersheet$IDH_Mutation <- ifelse(adult.braf.mastersheet$IDH == 0, "wt", NA)
adult.braf.mastersheet$IDH_Mutation <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% IDH_muts_canonical, "canonical", adult.braf.mastersheet$IDH_Mutation)
adult.braf.mastersheet$IDH_Mutation <- ifelse(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% IDH_muts_other, "other", adult.braf.mastersheet$IDH_Mutation)

## Other specific mutations
# adult.braf.mastersheet$ATRX <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "ATRX",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "ATRX",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$TP53 <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "TP53",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "TP53",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$CDKN2A.B <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "CDKN2A/B",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "CDKN2A/B",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$NF1 <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "NF1",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "NF1",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$PTEN <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "PTEN",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "PTEN",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$KIAA1549 <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "KIAA1549",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "KIAA1549",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$CIC <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "CIC",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "CIC",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$TERT <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "TERT",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "TERT",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$EGFR <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "EGFR",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "EGFR",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$FUBP1 <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "FUBP1",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "FUBP1",]$SAMPLE_ACCESSION_NBR)]
# 
# adult.braf.mastersheet$PEBP <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "PEBP",]$variant_classification)[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "PEBP",]$SAMPLE_ACCESSION_NBR)]

adult.braf.mastersheet$TMB <- classifier$mut_num[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]

adult.braf.mastersheet$Hypermutated <- classifier$hypermutated[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]

#adult.braf.mastersheet$TMB <- 'TBD'
#adult.braf.mastersheet$Hypermutated <- 'TBD'

adult.braf.mastersheet$BRAF_Class <- classifier$braf_class[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
adult.braf.mastersheet$BRAF_Class <- ifelse(is.na(adult.braf.mastersheet$BRAF_Class), (jhh_os_data$braf_mutation_class[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), adult.braf.mastersheet$BRAF_Class)
adult.braf.mastersheet$BRAF_Class[!(is.na(adult.braf.mastersheet$Hypermutated)) & adult.braf.mastersheet$Hypermutated == TRUE] <- "Hypermutated"

table((adult.braf.mastersheet$BRAF_Class))
#adult.braf.mastersheet$BRAF_Class <- ifelse(is.na(adult.braf.mastersheet$BRAF_Class), (jhh_os_data$braf_mutation_class[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), levels(adult.braf.mastersheet$BRAF_Class)[adult.braf.mastersheet$BRAF_Class])

adult.braf.mastersheet$BRAF_Class
#levels(adult.braf.mastersheet$BRAF_Class)[adult.braf.mastersheet$BRAF_Class]

adult.braf.mastersheet$BRAF_Mutation_Specific <- classifier$braf_mutation[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]

adult.braf.mastersheet$BRAF_Mutation_Specific <- ifelse(is.na(adult.braf.mastersheet$BRAF_Mutation_Specific), (jhh_os_data$braf_mutation_specify[match(adult.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), adult.braf.mastersheet$BRAF_Mutation_Specific)

## Pull from the old code
# braf.master.sheet$Therapy <- tcga.gbm.outcomes$Therapy.Class[match(
#   braf.master.sheet$SAMPLE_ACCESSION_NBR, tcga.gbm.outcomes$Case.ID)]
# 
# braf.master.sheet$Therapy <- ifelse(is.na(
#   braf.master.sheet$Therapy), dfci.req.treatment.plan$STD_CHEMO_PLAN[match(
#     braf.master.sheet$PATIENT_ID_19, dfci.req.treatment.plan$PATIENT_ID
#   )], braf.master.sheet$Therapy)
# 
# braf.master.sheet$DOTherapyStart <- format(strptime(as.character(dfci.req.treatment.plan$TPLAN_START_DT[match(braf.master.sheet$PATIENT_ID_19, dfci.req.treatment.plan$PATIENT_ID)]), "%d-%b-%Y"), "%Y-%m-%d")
# 
# braf.master.sheet$CycleDates <- paste(format(strptime(as.character(dfci.req.treatment.plan$FIRST_CYCLE_DT[match(braf.master.sheet$PATIENT_ID_19, dfci.req.treatment.plan$PATIENT_ID)]), "%d-%b-%Y"), "%Y-%m-%d"), format(strptime(as.character(dfci.req.treatment.plan$LAST_RCVD_CYCLE_DT[match(braf.master.sheet$PATIENT_ID_19, dfci.req.treatment.plan$PATIENT_ID)]), "%d-%b-%Y"), "%Y-%m-%d"), sep = ";")
# 
# braf.master.sheet$CycleDates[braf.master.sheet$CycleDates == "NA;NA"] <- NA
# 
# braf.master.sheet$DCReason <- dfci.req.treatment.plan$TPLAN_DC_REASON[match(braf.master.sheet$PATIENT_ID_19, dfci.req.treatment.plan$PATIENT_ID)]
# 
# braf.master.sheet$DCReason[braf.master.sheet$DCReason == "AUTOMATICALLY DISCONTINUED - PATIENT DECEASED"] <- "DECEASED"
# braf.master.sheet$DCReason[braf.master.sheet$DCReason == "PATIENT EXPIRED"] <- "DECEASED"
# 
# unique(braf.master.sheet$Histo_Dx)
# table(braf.master.sheet$BRAF_class[!is.na(braf.master.sheet$BRAF_class) & braf.master.sheet$BRAF_class == "Class I"])
# table(braf.master.sheet$BRAF_class[!is.na(braf.master.sheet$BRAF_class) & braf.master.sheet$BRAF_class == "Class II"])
# table(braf.master.sheet$BRAF_class[!is.na(braf.master.sheet$BRAF_class) & braf.master.sheet$BRAF_class == "Class III"])
# braf.master.sheet$SAMPLE_ACCESSION_NBR[braf.master.sheet$Histo_Dx == "Glioblastoma"]
# braf.master.sheet$SAMPLE_ACCESSION_NBR[braf.master.sheet$Histo_Dx == "Glioblastoma"]
# unique(braf.master.sheet$Histo_Dx)

rownames(adult.braf.mastersheet) <- adult.braf.mastersheet$SAMPLE_ACCESSION_NBR
rownames(adult.braf.mastersheet)

#write.table(pathology_classes, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/Data/classifications//braf_all_pathology_classifications.2021-04-02.csv", sep = ",")

### Peds BRAF data sheet ###
peds.braf.mastersheet <- classifier[classifier$age_int == "<18" ,c("SAMPLE_ACCESSION_NBR", "PATIENT_ID")]
peds.braf.mastersheet <- rbind(peds.braf.mastersheet, (jhh_os_data[jhh_os_data$age_at_dx <18.0, c("SAMPLE_ACCESSION_NBR", "redcap_patient_id")]))
peds.braf.mastersheet$PATIENT_ID_19 <- master.sheet.dfci$pt_id_19[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$Cohort <- classifier$cohort[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$Cohort[grepl("JHH", peds.braf.mastersheet$SAMPLE_ACCESSION_NBR)] <- "JHH"

peds.braf.mastersheet$Panel_Version <- classifier$PANEL_VERSION[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$Panel_Version[grepl("JHH", peds.braf.mastersheet$SAMPLE_ACCESSION_NBR)] <- "4.0"

peds.braf.mastersheet$Gender <- classifier$Gender[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$Gender <- ifelse(is.na(
  peds.braf.mastersheet$Gender), jhh_os_data$patient_gender[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR
  )], peds.braf.mastersheet$Gender)

peds.braf.mastersheet$Age <- classifier$age_num[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$Age <- ifelse(is.na(
  peds.braf.mastersheet$Age), jhh_os_data$age_at_dx[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], peds.braf.mastersheet$Age)

peds.braf.mastersheet$Age_Interval <- classifier$age_int[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$Age_Interval <- as.character(peds.braf.mastersheet$Age_Interval)
peds.braf.mastersheet$Age_Interval <- ifelse(is.na(
  peds.braf.mastersheet$Age_Interval), jhh_os_data$age_int2[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)], peds.braf.mastersheet$Age_Interval)

peds.braf.mastersheet$MRN <- master.sheet.dfci$MRN[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$MRN_BCH <- master.sheet.dfci$CHILDRENS_MRN[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$MRN_BCH[peds.braf.mastersheet$MRN_BCH == 'null'] <- NA 

peds.braf.mastersheet$DOB <- master.sheet.dfci$DOB[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$DODiagnosis <- ifelse(!(is.na(
 peds.braf.mastersheet$PATIENT_ID_19)), format(strptime(as.character(dfci.req.cancer.diagnosis.careg$DIAGNOSIS_DT[match(peds.braf.mastersheet$PATIENT_ID_19, dfci.req.cancer.diagnosis.careg$PATIENT_ID
  )]), "%d-%b-%Y"), "%Y-%m-%d"), tcga.lgg.outcomes$initial_pathologic_dx_year[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.lgg.outcomes$bcr_patient_barcode)])

peds.braf.mastersheet$DOSurgery <- ifelse(!(is.na(
  peds.braf.mastersheet$PATIENT_ID_19)), master.sheet.dfci$DoS[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR
  )], NA)

peds.braf.mastersheet$DOLastFU <- peds.outcomes.verified$DoLC[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, peds.outcomes.verified$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$DOD <- peds.outcomes.verified$DoD[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, peds.outcomes.verified$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$Alive <- peds.outcomes.verified$HYBRID_DEATH_IND[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, peds.outcomes.verified$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$Alive <- ifelse(is.na(
  peds.braf.mastersheet$Alive), tcga.lgg.outcomes$vital_status[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.lgg.outcomes$bcr_patient_barcode
  )], peds.braf.mastersheet$Alive)
peds.braf.mastersheet$Alive <- ifelse(is.na(
  peds.braf.mastersheet$Alive), tcga.gbm.outcomes$Vital.Status[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.gbm.outcomes$Case.ID)], peds.braf.mastersheet$Alive)
peds.braf.mastersheet$Alive <- ifelse(is.na(
  peds.braf.mastersheet$Alive), genie.outcomes$Patient.s.Vital.Status[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, genie.outcomes$SAMPLE_ACCESSION_NBR
  )], peds.braf.mastersheet$Alive)

peds.braf.mastersheet$Alive[peds.braf.mastersheet$Alive == "Y"] <- 0
peds.braf.mastersheet$Alive[peds.braf.mastersheet$Alive == "N"] <- 1
peds.braf.mastersheet$Alive[peds.braf.mastersheet$Alive == "DECEASED"] <- 0
peds.braf.mastersheet$Alive[peds.braf.mastersheet$Alive == "LIVING"] <- 1
peds.braf.mastersheet$Alive[peds.braf.mastersheet$Alive == "Dead"] <- 0
peds.braf.mastersheet$Alive[peds.braf.mastersheet$Alive == "Alive"] <- 1
peds.braf.mastersheet$Alive[peds.braf.mastersheet$Alive == "DEAD"] <- 0
peds.braf.mastersheet$Alive[peds.braf.mastersheet$Alive == "ALIVE"] <- 1

#peds.braf.mastersheet$Death <- peds.braf.mastersheet$Alive
#peds.braf.mastersheet$Death[peds.braf.mastersheet$Alive == "0"] <- 1
#peds.braf.mastersheet$Death[peds.braf.mastersheet$Alive == "1"] <- 0
#peds.braf.mastersheet$Death <- as.numeric(peds.braf.mastersheet$Death)

## OS 
peds.braf.mastersheet$OS <- (as.Date(as.character(peds.braf.mastersheet$DOD), format="%Y-%m-%d")-as.Date(as.character(peds.braf.mastersheet$DODiagnosis), format="%Y-%m-%d"))/30.4167

peds.braf.mastersheet$OS <- ifelse(is.na(peds.braf.mastersheet$OS), (as.Date(as.character(
  peds.braf.mastersheet$DOD), format="%Y-%m-%d")-as.Date(as.character(
    peds.braf.mastersheet$DOSurgery), format="%Y-%m-%d"))/30.4167, peds.braf.mastersheet$OS)

peds.braf.mastersheet$OS <- ifelse(is.na(peds.braf.mastersheet$OS) & peds.braf.mastersheet$Alive == 1, (as.Date(as.character(peds.braf.mastersheet$DOLastFU), format="%Y-%m-%d")-as.Date(as.character(
  peds.braf.mastersheet$DODiagnosis), format="%Y-%m-%d"))/30.4167, peds.braf.mastersheet$OS)

peds.braf.mastersheet$OS <- ifelse(is.na(peds.braf.mastersheet$OS) & peds.braf.mastersheet$Alive == 1, (as.Date(as.character(peds.braf.mastersheet$DOLastFU), format="%Y-%m-%d")- as.Date(as.character(
  peds.braf.mastersheet$DOSurgery), format="%Y-%m-%d"))/30.4167, peds.braf.mastersheet$OS)

peds.braf.mastersheet$OS <- ifelse(is.na(peds.braf.mastersheet$OS), (tcga.lgg.outcomes$OS.time[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.lgg.outcomes$bcr_patient_barcode
)])/30.4167, peds.braf.mastersheet$OS)

peds.braf.mastersheet$OS <- ifelse(is.na(peds.braf.mastersheet$OS), (tcga.gbm.outcomes$OS..days.[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.gbm.outcomes$Case.ID
)])/30.4167, peds.braf.mastersheet$OS)

peds.braf.mastersheet$OS[peds.braf.mastersheet$OS < 0] <- NA

peds.braf.mastersheet$OS <- ifelse((peds.braf.mastersheet$Alive == 1 & is.na(peds.braf.mastersheet$OS)), 0, peds.braf.mastersheet$OS)

peds.braf.mastersheet$OS <- ifelse(is.na(peds.braf.mastersheet$OS), as.double((
  genie.outcomes$Overall.Survival..Months.[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, genie.outcomes$SAMPLE_ACCESSION_NBR)])), peds.braf.mastersheet$OS)

#peds.braf.mastersheet$OS <- ifelse(is.na(peds.braf.mastersheet$OS), (jhh_os_data$overall_survival[match(
#  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), peds.braf.mastersheet$OS)

peds.braf.mastersheet$OS[peds.braf.mastersheet$OS < 0] <- NA
peds.braf.mastersheet$OS <- round(peds.braf.mastersheet$OS, 1)
table(is.na(adult.braf.mastersheet$OS))

## OS days
peds.braf.mastersheet$OS_days <- (as.Date(as.character(peds.braf.mastersheet$DOD), format="%Y-%m-%d")-
                                     as.Date(as.character(peds.braf.mastersheet$DODiagnosis), format="%Y-%m-%d"))
peds.braf.mastersheet$OS_days <- ifelse(is.na(peds.braf.mastersheet$OS_days), as.Date(as.character(
  peds.braf.mastersheet$DOD), format="%Y-%m-%d")- as.Date(as.character(
    peds.braf.mastersheet$DOSurgery), format="%Y-%m-%d"), peds.braf.mastersheet$OS_days)
peds.braf.mastersheet$OS_days <- ifelse(is.na(peds.braf.mastersheet$OS_days) & peds.braf.mastersheet$Alive == 1, (as.Date(as.character(peds.braf.mastersheet$DOLastFU), format="%Y-%m-%d")- as.Date(as.character(
  peds.braf.mastersheet$DODiagnosis), format="%Y-%m-%d")), peds.braf.mastersheet$OS_days)
peds.braf.mastersheet$OS_days <- ifelse(is.na(peds.braf.mastersheet$OS_days) & peds.braf.mastersheet$Alive == 1, (as.Date(as.character(peds.braf.mastersheet$DOLastFU), format="%Y-%m-%d")- as.Date(as.character(
  peds.braf.mastersheet$DOSurgery), format="%Y-%m-%d")), peds.braf.mastersheet$OS_days)
peds.braf.mastersheet$OS_days <- ifelse(is.na(peds.braf.mastersheet$OS_days), tcga.lgg.outcomes$OS.time[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.lgg.outcomes$bcr_patient_barcode
)], peds.braf.mastersheet$OS_days)
peds.braf.mastersheet$OS_days <- ifelse(is.na(peds.braf.mastersheet$OS_days), tcga.gbm.outcomes$OS..days.[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.gbm.outcomes$Case.ID
)], peds.braf.mastersheet$OS_days)
peds.braf.mastersheet$OS_days <- ifelse(is.na(peds.braf.mastersheet$OS_days), as.integer((
  genie.outcomes$Overall.Survival..Months.[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, genie.outcomes$SAMPLE_ACCESSION_NBR)])*30.4167), peds.braf.mastersheet$OS_days)
peds.braf.mastersheet$OS_days <- ifelse(is.na(peds.braf.mastersheet$OS_days), (jhh_os_data$os_days[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), peds.braf.mastersheet$OS_days)
peds.braf.mastersheet$OS_days[peds.braf.mastersheet$OS_days < 0] <- NA
peds.braf.mastersheet$OS_days <- round(peds.braf.mastersheet$OS_days, 0)

peds.braf.mastersheet$KPS <- NA

peds.braf.mastersheet$KPS_Date <- NA

peds.braf.mastersheet$Primary <- classifier$Primary[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$Primary[peds.braf.mastersheet$Primary == "TRUE"] <- 1
peds.braf.mastersheet$Primary[peds.braf.mastersheet$Primary == "FALSE"] <- 0

peds.braf.mastersheet$Tumor_Location <- dfci.req.cancer.diagnosis.careg$SITE_DESCR[match(peds.braf.mastersheet$PATIENT_ID_19, dfci.req.cancer.diagnosis.careg$PATIENT_ID)]
peds.braf.mastersheet$Tumor_Location <- sub(" LOBE", "", peds.braf.mastersheet$Tumor_Location)
peds.braf.mastersheet$Tumor_Location <-  sub(", NOS", "", peds.braf.mastersheet$Tumor_Location)
peds.braf.mastersheet$Tumor_Location <-  sub(", OVERLAPPING LESION", "", peds.braf.mastersheet$Tumor_Location)

peds.braf.mastersheet$Tumor_Laterality <- dfci.req.cancer.diagnosis.careg$LATERALITY_DESCR[match(peds.braf.mastersheet$PATIENT_ID_19, dfci.req.cancer.diagnosis.careg$PATIENT_ID)]
peds.braf.mastersheet$Tumor_Laterality[grepl("MIDLINE", peds.braf.mastersheet$Tumor_Laterality)] <- "MIDLINE"
peds.braf.mastersheet$Tumor_Laterality <- sub(" .*", "", peds.braf.mastersheet$Tumor_Laterality)
peds.braf.mastersheet$Tumor_Laterality[peds.braf.mastersheet$Tumor_Laterality == "NOT"] <- NA
peds.braf.mastersheet$Tumor_Laterality[peds.braf.mastersheet$Tumor_Laterality == "PAIRED"] <- NA
peds.braf.mastersheet$Tumor_Laterality[peds.braf.mastersheet$Tumor_Laterality == "BILATERAL"] <- "R;L"
peds.braf.mastersheet$Tumor_Laterality[grepl("MIDLINE", peds.braf.mastersheet$Tumor_Laterality)] <- "M"
peds.braf.mastersheet$Tumor_Laterality[peds.braf.mastersheet$Tumor_Laterality == "LEFT"] <- "L"
peds.braf.mastersheet$Tumor_Laterality[peds.braf.mastersheet$Tumor_Laterality == "RIGHT"] <- "R"

peds.braf.mastersheet$Extent_Resection <- NA

peds.braf.mastersheet$BRAF_Inhibitor <- peds.outcomes.verified$braf_inhibitor[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, peds.outcomes.verified$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$MGMT <- master.sheet.dfci$MGMT_Mehdi[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, master.sheet.dfci$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$MGMT <- ifelse(is.na(
  peds.braf.mastersheet$MGMT), tcga.gbm.outcomes$MGMT.Status[match(
    peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, tcga.gbm.outcomes$Case.ID)], peds.braf.mastersheet$MGMT)
peds.braf.mastersheet$MGMT[peds.braf.mastersheet$MGMT == "UNMETHYLATED"] <- 0
peds.braf.mastersheet$MGMT[peds.braf.mastersheet$MGMT == "METHYLATED"] <- 1
peds.braf.mastersheet$MGMT[peds.braf.mastersheet$MGMT == "PARTIAL"] <- 0.5
peds.braf.mastersheet$MGMT[!(peds.braf.mastersheet$MGMT == "0") & !(peds.braf.mastersheet$MGMT == "1") & !(peds.braf.mastersheet$MGMT == "0.5")] <- 2
peds.braf.mastersheet$MGMT[is.na(peds.braf.mastersheet$MGMT)] <- 2
peds.braf.mastersheet$MGMT[peds.braf.mastersheet$MGMT == "2"] <- NA

peds.braf.mastersheet$IDH <- classifier$IDH1.2[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$IDH[peds.braf.mastersheet$IDH == "TRUE"] <- 1
peds.braf.mastersheet$IDH[peds.braf.mastersheet$IDH == "FALSE"] <- 0

peds.braf.mastersheet$IDH <- ifelse(is.na(peds.braf.mastersheet$IDH), IDH_muts_df$MUT_TYPE[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, IDH_muts_df$SAMPLE_ACCESSION_NBR)], peds.braf.mastersheet$IDH)
peds.braf.mastersheet$IDH[peds.braf.mastersheet$IDH == "canonical"] <- 1
peds.braf.mastersheet$IDH[peds.braf.mastersheet$IDH == "other"] <- 1
peds.braf.mastersheet$IDH[is.na(peds.braf.mastersheet$IDH)] <- 0
peds.braf.mastersheet$IDH


#Ki67?
#ATRx?
#p53?
#GFAP?
peds.braf.mastersheet$x1p_19q <- classifier$x1p_19q[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$x1p_19q[peds.braf.mastersheet$x1p_19q == "FALSE"] <- 0
peds.braf.mastersheet$x1p_19q[peds.braf.mastersheet$x1p_19q == "TRUE"] <- 1
#peds.braf.mastersheet$x1p_19q_gen <- ifelse(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% c1p19qloss, 1, 0)

peds.braf.mastersheet$H3K27M <- classifier$H3K27M[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$H3K27M[peds.braf.mastersheet$H3K27M == "FALSE"] <- 0
peds.braf.mastersheet$H3K27M[peds.braf.mastersheet$H3K27M == "TRUE"] <- 1

peds.braf.mastersheet$WHO_Grade <- classifier$grade[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$WHO_Grade <- ifelse(is.na(peds.braf.mastersheet$WHO_Grade), (jhh_os_data$who_grade[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), peds.braf.mastersheet$WHO_Grade)

peds.braf.mastersheet$WHO_Diagnosis <- classifier$WHO_Diagnosis[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$WHO_Diagnosis <- ifelse(is.na(peds.braf.mastersheet$WHO_Diagnosis), (jhh_os_data$who_dx[match(
  peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, jhh_os_data$SAMPLE_ACCESSION_NBR)]), peds.braf.mastersheet$WHO_Diagnosis)

#peds.braf.mastersheet$Subtype_Path <- path_subtypes_who$Subtype_Path[match(peds.braf.mastersheet$WHO_Diagnosis, path_subtypes_who$WHO_Diagnosis)]
peds.braf.mastersheet$Subtype_Path <- path_subtypes_who$Subtype[match(peds.braf.mastersheet$WHO_Diagnosis, path_subtypes_who$WHO_Diagnosis)]

peds.braf.mastersheet$Subtype <- classifier$subtype[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$Subtype <- ifelse(is.na(peds.braf.mastersheet$Subtype), peds.braf.mastersheet$Subtype_Path, peds.braf.mastersheet$Subtype)

peds.braf.mastersheet$c7Gain_10Loss <- ifelse(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% c7gain10loss, 1, 0)

peds.braf.mastersheet$EGFR_Gain <- ifelse(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% egfrgain, 1, 0)

peds.braf.mastersheet$TERTP_Mut <- ifelse(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% tertpmut, 1, 0)

peds.braf.mastersheet$cIMPACT <- classifier$cIMPACT[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
table(peds.braf.mastersheet$cIMPACT)
peds.braf.mastersheet$IDH_Mutation <- ifelse(peds.braf.mastersheet$IDH == 0, "wt", NA)
peds.braf.mastersheet$IDH_Mutation <- ifelse(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% IDH_muts_canonical, "canonical", peds.braf.mastersheet$IDH_Mutation)
peds.braf.mastersheet$IDH_Mutation <- ifelse(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% IDH_muts_other, "other", peds.braf.mastersheet$IDH_Mutation)

# ## Other genes
# peds.braf.mastersheet$ATRX <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "ATRX",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "ATRX",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$TP53 <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "TP53",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "TP53",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$CDKN2A.B <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "CDKN2A/B",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "CDKN2A/B",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$NF1 <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "NF1",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "NF1",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$PTEN <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "PTEN",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "PTEN",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$KIAA1549 <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "KIAA1549",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "KIAA1549",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$CIC <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "CIC",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "CIC",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$TERT <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "TERT",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "TERT",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$EGFR <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "EGFR",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "EGFR",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$FUBP1 <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "FUBP1",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "FUBP1",]$SAMPLE_ACCESSION_NBR)]
# 
# peds.braf.mastersheet$PEBP <- (nmf_alts[nmf_alts$BEST_EFF_GENE == "PEBP",]$variant_classification)[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, nmf_alts[nmf_alts$BEST_EFF_GENE == "PEBP",]$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$TMB <- classifier$mut_num[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]

peds.braf.mastersheet$Hypermutated <- classifier$hypermutated[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]

#peds.braf.mastersheet$TMB <- 'TBD'
#peds.braf.mastersheet$Hypermutated <- 'TBD'

peds.braf.mastersheet$BRAF_Class <- classifier$braf_class[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]
peds.braf.mastersheet$BRAF_Class[!(is.na(peds.braf.mastersheet$Hypermutated)) & peds.braf.mastersheet$Hypermutated == TRUE] <- "Hypermutated"
table(peds.braf.mastersheet$BRAF_Class)

peds.braf.mastersheet$BRAF_Mutation_Specific <- classifier$braf_mutation[match(peds.braf.mastersheet$SAMPLE_ACCESSION_NBR, classifier$SAMPLE_ACCESSION_NBR)]

rownames(peds.braf.mastersheet) <- peds.braf.mastersheet$SAMPLE_ACCESSION_NBR

## Combine all
all.braf.mastersheet <- rbind(adult.braf.mastersheet, peds.braf.mastersheet)

##### end Clinical Table Parsing #####
colnames(peds.braf.mastersheet)[!(colnames(peds.braf.mastersheet) %in% colnames(adult.braf.mastersheet))]
colnames(adult.braf.mastersheet)[!(colnames(adult.braf.mastersheet) %in% colnames(peds.braf.mastersheet))]

## Write characteristics summary table by alteration type
table(all.braf.mastersheet$BRAF_Class)
table(adult.braf.mastersheet$BRAF_Class)
table(peds.braf.mastersheet$BRAF_Class)

for (bc in c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other")){
  print(bc)
  print(table(all.braf.mastersheet$Gender[all.braf.mastersheet$BRAF_Class == bc]))
}

for (bc in c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other")){
  print(bc)
  print(table(all.braf.mastersheet$Age_Interval[all.braf.mastersheet$BRAF_Class == bc]))
}

for (bc in c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other")){
  print(bc)
  print(table(all.braf.mastersheet$Subtype[all.braf.mastersheet$BRAF_Class == bc]))
}

for (bc in c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other")){
  print(bc)
  print(table(all.braf.mastersheet$WHO_Grade[all.braf.mastersheet$BRAF_Class == bc]))
}

for (bc in c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other")){
  print(bc)
  print(table(all.braf.mastersheet$Alive[all.braf.mastersheet$BRAF_Class == bc]))
}

survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.all)
survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.adults)
survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.peds)

for (bc in c("Class I", "Class II", "Class III", "Fusion", "Gain", "Hypermutated", "Other")){
  print(bc)
  print(table(peds.braf.mastersheet$Alive[peds.braf.mastersheet$BRAF_Class == bc]))
}

## Double checking WHO grade numbers discrepancy
temp_df <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/clinical-master-sheets/braf_all_glioma_mastersheet.2021-07-09.csv", stringsAsFactors = F)
table(temp_df$WHO_Grade)
table(all.braf.mastersheet$WHO_Grade)

table(temp_df$WHO_Grade[temp_df$BRAF_Class == "Class I"])
table(all.braf.mastersheet$WHO_Grade[all.braf.mastersheet$BRAF_Class == "Class I"])

table(censored.adult.braf.mastersheet$WHO_Grade[censored.adult.braf.mastersheet$BRAF_Class == "Class I"])
table(adult.braf.mastersheet$WHO_Grade[adult.braf.mastersheet$BRAF_Class == "Class I"])

table(censored.peds.braf.mastersheet$WHO_Grade[censored.peds.braf.mastersheet$BRAF_Class == "Class I"])
table(peds.braf.mastersheet$WHO_Grade[peds.braf.mastersheet$BRAF_Class == "Class I"])

table(all.braf.mastersheet$WHO_Diagnosis[all.braf.mastersheet$WHO_Grade == "4"])
## It appears to be related to JHH-52: Astro IDH-mt, considered to be Grade 4 vs Grade 2!!!

## Incorrect OS method (should use KM median and 95% CI)
#table(peds.braf.mastersheet$Alive[peds.braf.mastersheet$BRAF_Class == "Other"])
#sort(peds.braf.mastersheet$OS[peds.braf.mastersheet$BRAF_Class == "Fusion"])

#median(peds.braf.mastersheet$OS[peds.braf.mastersheet$BRAF_Class == "Other"], na.rm=TRUE)
#range(peds.braf.mastersheet$OS[peds.braf.mastersheet$BRAF_Class == "Other"], na.rm=TRUE)

#median(adult.braf.mastersheet$OS[adult.braf.mastersheet$BRAF_Class == "Class I"], na.rm=TRUE)
#range(adult.braf.mastersheet$OS[adult.braf.mastersheet$BRAF_Class == "Other"], na.rm=TRUE)

## Pathology classification testing
#sort(table(all.braf.mastersheet$WHO_Diagnosis))
#table(adult.braf.mastersheet$WHO_Diagnosis)
#pathology_classes <- all.braf.mastersheet[,c("WHO_Grade", "WHO_Diagnosis", "Subtype_Path", "Subtype_Molec", "IDH", "x1p_19q_gen", "cIMPACT", "c7Gain_10Loss", "EGFR_Gain", "TERTP_Mut")]
#pathology_classes <- unique(pathology_classes)
#pathology_classes_simple <- unique(all.braf.mastersheet[,c("WHO_Diagnosis", "Subtype_Path")])
#adult_path_classes_simple <- unique(adult.braf.mastersheet[,c("WHO_Diagnosis", "Subtype_Path")])

## Clarifying the additional 15 cases
#path_clarifying_list <- c("GENIE-JHU-03601-04216", "GENIE-MDA-6677-18354", "BL-15-T48446", "BL-16-A40851", "TCGA-DU-7309", "GENIE-JHU-02493-02971", "BL-14-D14982", "BL-15-E03366", "BL-15-G14321", "BL-15-X17861", "BL-16-W42664", "BL-17-M25111", "GENIE-MSK-P-0005609-T01-IM5", "TCGA-DU-7008", "TCGA-FG-8182")
#path_clarifying_list
#path_clarifying.mastersheet <- all.braf.mastersheet[all.braf.mastersheet$SAMPLE_ACCESSION_NBR %in% path_clarifying_list,]

## Write table
#write.table(adult.braf.mastersheet, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/Data/clinical-master-sheets/braf_adult_glioma_mastersheet.2021-08-11.csv", sep = ",")
#write.table(peds.braf.mastersheet, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/Data/clinical-master-sheets/braf_peds_glioma_mastersheet.2021-08-11.csv", sep = ",")
#write.table(all.braf.mastersheet, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/Data/clinical-master-sheets/braf_all_glioma_mastersheet.2021-08-11.csv", sep = ",")


#write.table(path_clarifying.mastersheet, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/Data/clinical-master-sheets/path_clarifying_info_mastersheet.2021-06-03.csv", sep = ",")

## Checking things randomly

###

table(all.braf.mastersheet$BRAF_Class)
table(adult.braf.mastersheet$BRAF_Class)

all.braf.mastersheet$SAMPLE_ACCESSION_NBR[all.braf.mastersheet$BRAF_Class == "Hypermutated"]
all.braf.mastersheet$SAMPLE_ACCESSION_NBR[all.braf.mastersheet$Hypermutated == TRUE]

#write.table(all.braf.mastersheet[all.braf.mastersheet$BRAF_Class == "Fusion",], "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/clinical-master-sheets/braf_all_glioma_fusions.2022-07-30.csv", sep = ",")

