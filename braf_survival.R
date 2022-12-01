## Pinky Langat 2020
## This code
## This code is built from tutorial https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html

#install.packages("survminer")

library(survival)
library(survminer)
library(lubridate)
library(RColorBrewer)

braf_surv.all <- all.braf.mastersheet[!is.na(all.braf.mastersheet$OS),]
braf_surv.adults <- adult.braf.mastersheet[!is.na(adult.braf.mastersheet$OS),]
braf_surv.peds <- peds.braf.mastersheet[!is.na(peds.braf.mastersheet$OS),]

## Preparing survial dataframes
braf_surv.all$Death <- ifelse(braf_surv.all$Alive == 0, 1, 0)
braf_surv.all$OS_yrs <- round(braf_surv.all$OS/12.0, 2)

braf_surv.adults$Death <- ifelse(braf_surv.adults$Alive == 0, 1, 0)
braf_surv.adults$OS_yrs <- round(braf_surv.adults$OS/12.0, 2)

braf_surv.peds$Death <- ifelse(braf_surv.peds$Alive == 0, 1, 0)
braf_surv.peds$OS_yrs <- round(braf_surv.peds$OS/12.0, 2)

braf_surv.adult_gbm <- braf_surv.adults[braf_surv.adults$Subtype == "GBM, IDH-wt",]

## Testing code
lung
Surv(lung$time, lung$status)[1:10]

braf_surv
Surv(braf_surv$time, braf_surv$status)
f1 <- survfit(Surv(time,status) ~1, data = braf_surv)
names(f1)

braf_survival
Surv(braf_survival$OS, braf_survival$Death)
table(braf_survival)

table(braf_surv$class)
braf_surv$class[braf_surv$class == "Class II"] <- "Class II/III"
braf_surv$class[braf_surv$class == "Class III"] <- "Class II/III"
braf_surv$class[braf_surv$class == "Rearrangement"] <- "Other"
table(braf_surv$class)

plot(survfit(Surv(time,status) ~class, data = braf_surv), xlab = "Days", ylab = "Overall survival probability")
ggsurvplot(
  fit = survfit(Surv(time, status) ~class, data = braf_surv),
  xlab = "Days",
  ylab = "Overall survival probability"#,
  #conf.int = TRUE
)

survfit(formula = Surv(time, status) ~class, data = braf_surv)
survdiff(formula = Surv(time, status) ~class, data = braf_surv)
sd <- survdiff(Surv(time, status) ~class, data = braf_surv)
1 - pchisq(sd$chisq, length(sd$n) - 1)

## cIMPACT vs non-cIMPACT GBM
plot(survfit(Surv(OS,Death) ~Hypermutated, data = braf_surv.adult_gbm), xlab = "Months", ylab = "Overall survival probability")

## Testing complete

#### Starts here ####
## Median overal survival statistics
survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.all)
survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.adults)
survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.peds)
survfit(Surv(OS, Death) ~1, data = braf_surv.adults[braf_surv.adults$BRAF_Class == "Fusion",])
survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults)
survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.adult_gbm)
survfit(Surv(OS, Death) ~1, data = braf_surv.adult_gbm)
survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adults)
survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adult_gbm)

## All pts by BRAF
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.all),
  title = "All patients with BRAF-altered glioma",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray50", "gray")#+
  #scale_fill_brewer(type = "qual", palette = "Spectral")
)

## All pts by BRAF with simplified groups
braf_surv.all$BRAF_Class_Groups <- braf_surv.all$BRAF_Class
braf_surv.all$BRAF_Class_Groups[braf_surv.all$BRAF_Class_Groups == "Hypermutated"] <- "Other"

#braf_surv.all$BRAF_Class_Groups[braf_surv.all$BRAF_Class_Groups == "Class II"] <- "Class II/III"
#braf_surv.all$BRAF_Class_Groups[braf_surv.all$BRAF_Class_Groups == "Class III"] <- "Class II/III"

ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.all),
  title = "All patients with BRAF-altered glioma",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)
pairwise_survdiff(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.all)

## All GBM by Class Groups
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.all[grepl("GBM", braf_surv.all$Subtype),]),
  title = "All patients with BRAF-altered GBM",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)

## All GBM by Class I vs not
braf_surv.all$BRAF_Class_Groups <- "Other"
braf_surv.all$BRAF_Class_Groups[braf_surv.all$BRAF_Class == "Class I"] <- "Class I"

ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.all[grepl("GBM", braf_surv.all$Subtype),]),
  title = "All patients BRAF-altered GBM",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  palette = c(brewer.pal(6, "Spectral")[c(6)], "gray")#+
)

## All Adults by BRAF
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.adults),
  title = "Adults",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray50", "gray")#+
  #scale_fill_brewer(type = "qual", palette = "Spectral")
)

surv_median(fit = survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.adults), combine = FALSE)
survfit(Surv(OS, Death) ~BRAF_Class, data = braf_surv.adults)

## All Adults by BRAF with simplified groups
braf_surv.adults$BRAF_Class_Groups <- braf_surv.adults$BRAF_Class
braf_surv.adults$BRAF_Class_Groups[braf_surv.adults$BRAF_Class_Groups == "Hypermutated"] <- "Other"

#braf_surv.adults$BRAF_Class_Groups[braf_surv.adults$BRAF_Class_Groups == "Class II"] <- "Class II/III"
#braf_surv.adults$BRAF_Class_Groups[braf_surv.adults$BRAF_Class_Groups == "Class III"] <- "Class II/III"

ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults),
  title = "Adults with BRAF-altered glioma",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)

df_forest.adult$braf_class_group <- df_forest.adult$braf_class
df_forest.adult$braf_class_group[df_forest.adult$braf_class_group == 'Hypermutated'] <- "Other"
ggsurvplot(
  fit = survfit(Surv(os, deceased) ~braf_class_group, data = df_forest.adult),
  title = "Adults with BRAF-altered glioma with molecular data",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)
pairwise_survdiff(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults, p.adjust.method = "BH")

## All Adults by BRAF Class I vs not
braf_surv.adults$BRAF_Class_Groups <- "Non-Class I"
braf_surv.adults$BRAF_Class_Groups[braf_surv.adults$BRAF_Class == "Class I"] <- "Class I"

ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults),
  title = "Adults with BRAF-altered glioma",
  xlab = "Time (months)",
  ylab = "Overall survival probability (%)",
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  fun = "pct",
  xlim = c(0,280),
  break.time.by = 48,
  #xscale = "m_y",
  #legend = "bottom",
  risk.table.height = 0.3,
  legend = "none",
  legend.labs = c("Class I (V600E)", "Non-Class I"),
  #surv.median.line = "hv",
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")
  #palette = c(brewer.pal(6, "Spectral")[c(6)], "darkgray")
  palette = c(brewer.pal(6, "Spectral")[c(6)], "gray50")
)

#p4a <-
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults),
  title = "Adult glioma",
  xlab = "Time (months)",
  ylab = "Overall survival probability (%)",
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  fun = "pct",
  #xlim = c(0,280),
  #break.time.by = 48,
  #xscale = "m_y",
  #legend = "bottom",
  #risk.table.height = 0.25,
  #surv.scale = "percent",
  base_size = 20,
  font.main = 20,
  font.legend = 20,
  font.x = 20,
  font.y = 20,
  font.tickslab = 20,
  pval.size = 7,
  #risk.table.col = "strata",
  risk.table.fontsize = 7,
  risk.table.title.fontsize = 7,
  risk.table.y.text = FALSE,
  risk.table.y.col = TRUE,
  legend = "right",
  legend.title = "BRAF Alteration",
  legend.labs = c("Class I (V600E)", "Non-Class I"),
  #surv.median.line = "hv",
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")
  #palette = c(brewer.pal(6, "Spectral")[c(6)], "darkgray")
  palette = c(brewer.pal(6, "Spectral")[c(6)], "gray50")
)  

pairwise_survdiff(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults)
surv_median(fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults), combine = FALSE)

brewer.pal("grey90")

## Adult GBM by BRAF class groups
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults[grepl("GBM", braf_surv.adults$Subtype),]),
  title = "Adults with BRAF-altered GBM",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)

## Adult GBM by Class I vs not
braf_surv.adults$BRAF_Class_Groups <- "Other"
braf_surv.adults$BRAF_Class_Groups[braf_surv.adults$BRAF_Class == "Class I"] <- "Class I"

#p5b <- 
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults[grepl("GBM", braf_surv.adults$Subtype),]),
  title = "BRAF-altered Adult GBM",
  xlab = "Time (months)",
  ylab = "Overall survival probability (%)",
  #xlim = c(0,10),
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  #risk.table.height = 0.25,
  surv.scale = "percent",
  base_size = 20,
  font.main = 20,
  font.legend = 20,
  font.x = 20,
  font.y = 20,
  font.tickslab = 20,
  pval.size = 7,
  #risk.table.col = "strata",
  risk.table.fontsize = 7,
  risk.table.title.fontsize = 7,
  risk.table.y.text = FALSE,
  risk.table.y.col = TRUE,
  legend = "right",
  legend.title = "BRAF Alteration",
  legend.labs = c("Class I", "Non-Class I"),
  #palette = c(brewer.pal(3, "Set1")[c(3,2,1)])
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  palette = c(brewer.pal(6, "Spectral")[c(6)], "gray")#+
  #palette = c("#2b83ba", "gray")
) 
pairwise_survdiff(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults)


## Adult Astro by BRAF class groups
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.adults[grepl("Astro", braf_surv.adults$Subtype),]),
  title = "Adults with BRAF-altered Astrocytoma, IDH-mt",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)

## Adult Fusions by Subtype
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~Subtype, data = braf_surv.adults[grepl("Fusion", braf_surv.adults$BRAF_Class),]),
  title = "Adult gliomas with BRAF Fusions",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#a6d854", "black", "#ffed6f", "skyblue")
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)

## Adult Fusions by Grade
#p5a <- 
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~WHO_Grade, data = braf_surv.adults[grepl("Fusion", braf_surv.adults$BRAF_Class),]),
  title = "BRAF-rearranged Adult Gliomas",
  xlab = "Time (months)",
  ylab = "Overall survival probability (%)",
  #xlim = c(0,10),
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  #risk.table.height = 0.25,
  surv.scale = "percent",
  base_size = 20,
  font.main = 20,
  font.legend = 20,
  font.x = 20,
  font.y = 20,
  font.tickslab = 20,
  pval.size = 7,
  #risk.table.col = "strata",
  risk.table.fontsize = 7,
  risk.table.title.fontsize = 7,
  risk.table.y.text = FALSE,
  risk.table.y.col = TRUE,
  legend = "right",
  legend.title = "WHO Grade",
  legend.labs = c("Grade 1", "Grade 2", "Grade 4"),  
  palette = c(brewer.pal(3, "Set1")[c(3,2,1)])
  #palette = c("#a6d854", "black", "#ffed6f", "skyblue")
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)
# scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "red", "#a6d854", "#ffed6f", "lightsalmon", "skyblue", "violet"))

# WHO Grade Groups 1 vs 2-4
braf_surv.adults$WHO_Grade_Groups[braf_surv.adults$WHO_Grade == 1] <- "1"
braf_surv.adults$WHO_Grade_Groups[braf_surv.adults$WHO_Grade %in% c(2,3,4)] <- "2-4"
braf_surv.adults$WHO_Grade_Groups[braf_surv.adults$WHO_Grade == "4" & braf_surv.adults$BRAF_Class == "Fusion"]
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~WHO_Grade_Groups, data = braf_surv.adults[grepl("Fusion", braf_surv.adults$BRAF_Class),]),
  title = "BRAF-rearranged Adult Gliomas",
  xlab = "Time (months)",
  ylab = "Overall survival probability (%)",
  #xlim = c(0,10),
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  #risk.table.height = 0.25,
  surv.scale = "percent",
  base_size = 20,
  font.main = 20,
  font.legend = 20,
  font.x = 20,
  font.y = 20,
  font.tickslab = 20,
  pval.size = 7,
  #risk.table.col = "strata",
  risk.table.fontsize = 7,
  risk.table.title.fontsize = 7,
  risk.table.y.text = FALSE,
  risk.table.y.col = TRUE,
  legend = "right",
  legend.title = "WHO Grade",
  legend.labs = c("Grade 1", "Grade 2-4"),  
  palette = c(brewer.pal(3, "Set1")[c(2,1)])
  #palette = c("#a6d854", "black", "#ffed6f", "skyblue")
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)

## All peds by BRAF class
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.peds),
  title = "Pediatric BRAF-altered glioma",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,4,2,1)], "gray50", "gray")#+
  #scale_fill_brewer(type = "qual", palette = "Spectral")
)
pairwise_survdiff(Surv(OS, Death) ~BRAF_Class, data = braf_surv.peds)

## All peds by BRAF simplified class
braf_surv.peds$BRAF_Class_Groups <- braf_surv.peds$BRAF_Class
braf_surv.peds$BRAF_Class_Groups[braf_surv.peds$BRAF_Class_Groups == "Hypermutated"] <- "Other"

ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.peds),
  title = "Pediatric BRAF-altered glioma",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,4,2,1)], "gray")#+
  #scale_fill_brewer(type = "qual", palette = "Spectral")
)
pairwise_survdiff(Surv(OS, Death) ~BRAF_Class_Groups, data = braf_surv.peds)

## Peds Fusions by Subtype
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~Subtype, data = braf_surv.peds[grepl("Fusion", braf_surv.peds$BRAF_Class),]),
  title = "Pediatric gliomas with BRAF Fusions",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#a6d854", "black", "#ffed6f", "skyblue")
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray")#+
  #palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
)

## Adults & Peds all gliomas by certain BRAF alterations
#braf_surv.all$Age_BRAF <- 
braf_surv.all$BRAF_Class_Groups <- braf_surv.all$BRAF_Class
#braf_surv.all$BRAF_Class_Groups[braf_surv.all$BRAF_Class_Groups == "Hypermutated"] <- "Other"
braf_surv.all$BRAF_Class_Groups[braf_surv.all$BRAF_Class_Groups %in% c("Hypermutated", "Class II", "Class III", "Gain")] <- "Other"

braf_surv.all$Age_Group <- "Adult"
braf_surv.all$Age_Group[braf_surv.all$Age_Interval == "<18"] <- "Pediatric"

braf_surv.all$Age_BRAF_Group <- paste(braf_surv.all$Age_Group, braf_surv.all$BRAF_Class_Groups)
braf_surv.all$Age_BRAF_Group[!(braf_surv.all$BRAF_Class_Groups %in% c("Class I", "Fusion", "Other"))] <- NA

braf_surv.all[!is.na(braf_surv.all$Age_BRAF_Group),]
table(braf_surv.all$Age_BRAF_Group)

ggsurvplot(
  fit = survfit(Surv(OS, Death) ~Age_BRAF_Group, data = braf_surv.all[!is.na(braf_surv.all$Age_BRAF_Group),]),
  title = "Adult vs pediatric gliomas by BRAF alteration type",
  xlab = "Time (months)",
  ylab = "Overall survival probability (%)",
  #xlim = c(0,10),
  pval = TRUE,
  #pval.method = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  fun = "pct",
  xlim = c(0,280),
  break.time.by = 48,
  #xscale = "m_y",
  linetype = c(1,1,1,3,3,3),
  palette = c(brewer.pal(12, "Paired")[c(1,7)], "gray", brewer.pal(6, "Spectral")[c(6,2)], "darkgray"),
  legend.labs = c("Adult Class I", "Adult Fusion", "Adult Other", "Ped Class I", "Ped Fusion", "Ped Other"),
  risk.table.height = 0.3,
  legend = "none"
  #aspect.ratio = 1
  #ggtheme = theme(aspect.ratio = 1) #+ theme_classic()
  #ggtheme = theme(aspect.ratio = 1)
)
# + ggtheme = theme(aspect.ratio = 1)
p
p$plot + theme_void()
#p$plot + scale_linetype_discrete(name = 'Age_Group', breaks=c(1,3), labels = c('Adult', 'Pediatric'))

brewer.pal(n=6,"Spectral")

braf_surv.all$Age_BRAF_Group[!(braf_surv.all$BRAF_Class_Groups %in% c("Class I", "Fusion"))] <- NA
braf_surv.all[!is.na(braf_surv.all$Age_BRAF_Group),]
table(braf_surv.all$Age_BRAF_Group)

#p4f <-
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~Age_BRAF_Group, data = braf_surv.all[!is.na(braf_surv.all$Age_BRAF_Group),]),
  title = "Adult and pediatric glioma",
  xlab = "Time (months)",
  ylab = "Overall survival probability (%)",
  #xlim = c(0,10),
  pval = TRUE,
  #pval.method = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  fun = "pct",
  #xlim = c(0,280),
  #break.time.by = 48,
  #xscale = "m_y",
  base_size = 20,
  font.main = 20,
  font.legend = 20,
  font.x = 20,
  font.y = 20,
  font.tickslab = 20,
  pval.size = 7,
  #risk.table.col = "strata",
  risk.table.fontsize = 7,
  risk.table.title.fontsize = 7,
  risk.table.y.text = FALSE,
  risk.table.y.col = TRUE,
  legend = "right",
  legend.title = "BRAF Alteration",
  linetype = c(1,1,3,3),
  #palette = c(brewer.pal(12, "Paired")[c(1,7)], "gray", brewer.pal(6, "Spectral")[c(6,2)], "darkgray"),
  palette = c(brewer.pal(12, "Paired")[c(1,7)], brewer.pal(6, "Spectral")[c(6,2)]),
  legend.labs = c("Adult Class I", "Adult Rearrangement", "Ped Class I", "Ped Rearrangement")
  #aspect.ratio = 1
  #ggtheme = theme(aspect.ratio = 1) #+ theme_classic()
  #ggtheme = theme(aspect.ratio = 1)
)

## Adults by BRAFi use
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adults),
  title = "All adults with BRAF-altered glioma",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#377eb8", "#e41a1c")
  #palette = c(brewer.pal(4, "RdBu")[c(4,1)])
)

## Adult Grade3-4 by BRAFi
braf_surv.adults$WHO_Grade[braf_surv.adults$WHO_Grade >2]
table(braf_surv.adults$WHO_Grade)
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adults[braf_surv.adults$WHO_Grade >2,]),
  title = "Adult high-grade glioma (Grade 3-4)",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#377eb8", "#e41a1c")
  #palette = c(brewer.pal(4, "RdBu")[c(4,1)])
)

# p4f <- 
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adults[braf_surv.adults$WHO_Grade >2,]),
  title = "Adult high-grade glioma (Grade 3-4)",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  fun = "pct",
  base_size = 20,
  font.main = 20,
  font.legend = 20,
  font.x = 20,
  font.y = 20,
  font.tickslab = 20,
  pval.size = 7,
  risk.table.fontsize = 7,
  risk.table.title.fontsize = 7,
  risk.table.y.text = FALSE,
  risk.table.y.col = TRUE,
  legend = "right",
  legend.title = "Targeted therapy",
  legend.labs = c("No BRAF/MEKi", "BRAF/MEKi"),
  palette = c("#e41a1c", "#377eb8")
  #palette = c(brewer.pal(4, "RdBu")[c(4,1)])
)

## Adult Class I by BRAFi
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adults[braf_surv.adults$BRAF_Class=="Class I",]),
  title = "Adults with Class I BRAF-altered glioma",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#377eb8", "#e41a1c")
  #palette = c(brewer.pal(4, "RdBu")[c(4,1)])
)

## Adult GBM BRAFi vs not
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adult_gbm),
  title = "Adults with BRAF-altered GBM",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#377eb8", "#e41a1c")
  #palette = c(brewer.pal(4, "RdBu")[c(4,1)])
)

ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adults[grepl("GBM", braf_surv.adults$Subtype),]),
  title = "Adults with BRAF-altered GBM",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#377eb8", "#e41a1c")
  #palette = c(brewer.pal(4, "RdBu")[c(4,1)])
)

## Adult GBM Class I BRAFi vs not
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adults[grepl("GBM", braf_surv.adults$Subtype) & braf_surv.adults$BRAF_Class=="Class I",]),
  title = "Adults with Class I BRAF-altered GBM",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#377eb8", "#e41a1c")
  #palette = c(brewer.pal(4, "RdBu")[c(4,1)])
)

## Adult PXA BRAFi vs not
ggsurvplot(
  fit = survfit(Surv(OS, Death) ~BRAF_Inhibitor, data = braf_surv.adults[grepl("PXA", braf_surv.adults$Subtype),]),
  title = "Adults with BRAF-altered PXA",
  xlab = "Survival (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#377eb8", "#e41a1c")
  #palette = c(brewer.pal(4, "RdBu")[c(4,1)])
)

## Adult GBM BRAF Class Groups
ggsurvplot(
  fit = survfit(Surv(OS,Death) ~BRAF_Class_Groups, data = braf_surv.adult_gbm),
  xlab = "Months",
  ylab = "Overall survival probability",
  pval = TRUE,
  risk.table = TRUE#,
  #conf.int = TRUE
)

## BRAF Classes
ggsurvplot(
  fit = survfit(Surv(OS,Death) ~BRAF_Class, data = braf_surv.adult_gbm),
  #title = "GBM Class I vs Class II/III/Fusions",
  xlab = "Months",
  ylab = "Overall survival probability",
  #pval = TRUE,
  risk.table = TRUE#,
  #conf.int = TRUE
)

## BRAF Class Groups
ggsurvplot(
  fit = survfit(Surv(OS,Death) ~BRAF_Class_Groups, data = braf_surv.adult_gbm),
  #title = "GBM Class I vs Class II/III/Fusions",
  xlab = "Months",
  ylab = "Overall survival probability",
  pval = TRUE,
  palette = c("#00BFC4","#7CAE00","#F8766D"),
  risk.table = TRUE#,
  #conf.int = TRUE
)

## Age intervals
ggsurvplot(
  fit = survfit(Surv(OS,Death) ~Age_Interval, data = braf_surv.adult_gbm),
  title = "GBM",
  xlab = "Months",
  ylab = "Overall survival probability",
  pval = TRUE,
  risk.table = TRUE#,
  #conf.int = TRUE
)

## cIMPACT
ggsurvplot(
  fit = survfit(Surv(OS,Death) ~cIMPACT, data = braf_surv.adult_gbm),
  title = "GBM",
  xlab = "Months",
  ylab = "Overall survival probability",
  pval = TRUE,
  risk.table = TRUE#,
  #conf.int = TRUE
)

## By Class
class_survival_data <- braf_surv.adults[braf_surv.adults$BRAF_Class == "Class III" & braf_surv.adults$Subtype == "GBM, IDH-wt",]
ggsurvplot(
  fit = survfit(Surv(OS_yrs,Death) ~BRAF_Class, data = class_survival_data),
  title = "Class III",
  xlab = "Survival (years)",
  ylab = "Overall survival probability",
  xlim = c(0,10),
  #pval = TRUE,
  #risk.table = TRUE#,
  conf.int = FALSE
)

## All GBM classes
braf_surv.adult_gbm$BRAF_Class_Groups <- braf_surv.adult_gbm$BRAF_Class
braf_surv.adult_gbm$BRAF_Class_Groups[braf_surv.adult_gbm$BRAF_Class_Groups == "Class II"] <- "Class II/III/Fusion"
braf_surv.adult_gbm$BRAF_Class_Groups[braf_surv.adult_gbm$BRAF_Class_Groups == "Class III"] <- "Class II/III/Fusion"
braf_surv.adult_gbm$BRAF_Class_Groups[braf_surv.adult_gbm$BRAF_Class_Groups == "Fusion"] <- "Class II/III/Fusion"
braf_surv.adult_gbm$BRAF_Class_Groups[braf_surv.adult_gbm$BRAF_Class_Groups == "Gain"] <- "Other"
braf_surv.adult_gbm$BRAF_Class_Groups[braf_surv.adult_gbm$BRAF_Class_Groups == "Hypermutated"] <- "Other"

ggsurvplot(
  fit = survfit(Surv(OS_yrs, Death) ~BRAF_Class, data = braf_surv.adult_gbm_simplified),
  title = "Adult GBM",
  xlab = "Survival (years)",
  ylab = "Overall survival probability",
  xlim = c(0,10),
  #pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
  #scale_fill_brewer(type = "qual", palette = "Spectral")
)

## Adult GBM by Age
braf_surv.adult_gbm[c("Age", "Age_Interval")]
factor(braf_surv.adult_gbm$Age_Interval)
braf_surv.adult_gbm$Age_Interval <- factor(braf_surv.adult_gbm$Age_Interval, levels = c("18-34", "35-50", ">50"))
#braf_surv.adult_gbm_simplified$Age_Interval
braf_surv.adult_gbm$OS
#p5c <- 
ggsurvplot(
  #fit = survfit(Surv(OS_yrs, Death) ~Age_Interval, data = braf_surv.adult_gbm_simplified),
  fit = survfit(Surv(OS, Death) ~Age_Interval, data = braf_surv.adult_gbm),
  title = "BRAF-altered Adult GBM",
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  #xlim = c(0,10),
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  #risk.table.height = 0.25,
  surv.scale = "percent",
  base_size = 20,
  font.main = 20,
  font.legend = 20,
  font.x = 20,
  font.y = 20,
  font.tickslab = 20,
  pval.size = 7,
  #risk.table.col = "strata",
  risk.table.fontsize = 7,
  risk.table.title.fontsize = 7,
  risk.table.y.text = FALSE,
  risk.table.y.col = TRUE,
  legend = "right",
  legend.title = "Patient Age (years)",
  #legend.labs = c("18-34", "35-50", ">50"),
  palette = c(brewer.pal(3, "Set1")[c(3,2,1)])
  #palette = c(brewer.pal(8, "OrRd")[c(7,3,5)])#+
  #scale_fill_brewer(type = "qual", palette = "Spectral")
)
c(brewer.pal(8, "OrRd")[c(7,3,5)])
c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")
survfit(formula = Surv(OS, Death) ~cIMPACT, data = braf_surv.adult_gbm)
survdiff(formula = Surv(OS, Death) ~cIMPACT, data = braf_surv.adult_gbm)
sd <- survdiff(Surv(OS, Death) ~cIMPACT, data = braf_surv.adult_gbm)
1 - pchisq(sd$chisq, length(sd$n) - 1)

## Adult GBM BRAF Class Groups
ggsurvplot(
  fit = survfit(Surv(OS_yrs, Death) ~BRAF_Class_Groups, data = braf_surv.adult_gbm),
  title = "Adult GBM",
  xlab = "Survival (years)",
  ylab = "Overall survival probability",
  xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c('#619FFF', '#00BA38', '#F8766D')
)

## GBM cIMPACT
ggsurvplot(
  fit = survfit(Surv(OS_yrs, Death) ~cIMPACT, data = braf_surv.adult_gbm),
  title = "Adult GBM",
  xlab = "Survival (years)",
  ylab = "Overall survival probability",
  xlim = c(0,10),
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  #palette = c('#619FFF', '#00BA38', '#F8766D')
)

## Adult GBM BRAF class 
braf_surv.adult_gbm
braf_surv.adult_gbm_simplified <- braf_surv.adult_gbm
braf_surv.adult_gbm_simplified$BRAF_Class[braf_surv.adult_gbm_simplified$BRAF_Class == "Class II"] <- "Class II/III"
braf_surv.adult_gbm_simplified$BRAF_Class[braf_surv.adult_gbm_simplified$BRAF_Class == "Class III"] <- "Class II/III"
ggsurvplot(
  fit = survfit(Surv(OS_yrs, Death) ~BRAF_Class, data = braf_surv.adult_gbm_simplified),
  title = "Adult GBM",
  xlab = "Survival (years)",
  ylab = "Overall survival probability",
  xlim = c(0,10),
  #pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c(brewer.pal(6, "Spectral")[c(6,5,2,1)], "gray")#+
  #scale_fill_brewer(type = "qual", palette = "Spectral")
)

## Make a layout page
require("survminer")
splots <- list()
splots[[1]] <- p5a
splots[[2]] <- p5b
splots[[3]] <- p5c

# Arrange multiple ggsurvplots and print the output
layout <- arrange_ggsurvplots(splots, print = TRUE, ncol = 2, nrow = 3)
ggsave("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/5-Figures/drafting-figures/fig5/Fig5_R_draft.pdf", layout, width = 8, height = 10, units = "in")


#fit.coxph <- coxph(surv_object ~ rx + resid.ds + age_group + ecog.ps, 
 #                  data = ovarian)
#ggforest(fit.coxph, data = ovarian)
braf_surv.adult_gbm 
