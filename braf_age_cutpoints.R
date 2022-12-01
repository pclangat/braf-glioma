## Pinky Langat 2020
## This code is designed to run after braf_oncoprint.R
## Has supporting analyses
#install.packages("logisticPCA")
library(logisticPCA)

library(tidyverse)
library(tidytidbits)
library(survivalAnalysis) 


### OS data and cutpoint analysis ###
# Standardizing and combining inst'l OS data
names(jhh_os_data)[names(jhh_os_data) == "age_at_dx"] <- "age_num"
names(jhh_os_data)[names(jhh_os_data) == "redcap_patient_id"] <- "membership"
names(jhh_os_data)[names(jhh_os_data) == "braf_mutation_class"] <- "braf_class"
names(jhh_os_data)[names(jhh_os_data) == "who_grade"] <- "grade"
names(jhh_os_data)[names(jhh_os_data) == "patient_gender"] <- "Gender"

gbm_os_data <- rbind((classifier[classifier$subtype_path == "GBM" & classifier$membership %in% sub_samps, c("membership", "age_num", "death", "os_days", "age_int", "Gender", "subtype_path", "grade", "braf_class")]), 
                     (jhh_os_data[jhh_os_data$subtype_path == "GBM", c("membership","age_num","death","os_days", "age_int", "Gender", "subtype_path", "grade", "braf_class")]))

names(gbm_os_data)[names(gbm_os_data) == "Gender"] <- "gender"
gbm_os_data$os_days <- as.numeric(gbm_os_data$os_days)
#gbm_os_data <- gbm_os_data[!is.na(gbm_os_data$os_days),]
count(gbm_os_data[gbm_os_data$age_num > 64 & !is.na(gbm_os_data$os_days),])

# Cutpoint analysis
gbm_age_cutpt <- surv_cutpoint(gbm_os_data,
                               time = "os_days",
                               event = "death",
                               variables = c("age_num")
)

summary(gbm_age_cutpt)
plot(gbm_age_cutpt, "age_num", palette = "npg")

gbm_age_cutpt_cat <- surv_categorize(gbm_age_cutpt)

fit <- surv_fit(Surv(os_days, death) ~ age_num,
                data = gbm_age_cutpt_cat)
ggsurvplot(
  fit,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for
  # point estimaes of survival curves.
  #xlim = c(18,90),        # present narrower X axis, but not affect
  # survival estimates.
  #break.time.by = 10,    # break X axis in time intervals by 500.
  ggtheme = theme_classic2(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)

# Cutpoint 2 analysis
gbm_age_cutpt2 <- surv_cutpoint(gbm_os_data[gbm_os_data$age_num < 52,],
                                time = "os_days",
                                event = "death",
                                variables = c("age_num")
)

summary(gbm_age_cutpt2)
plot(gbm_age_cutpt2, "age_num", palette = "npg")

gbm_age_cutpt_cat2 <- surv_categorize(gbm_age_cutpt2)

fit2 <- surv_fit(Surv(os_days, death) ~ age_num,
                 data = gbm_age_cutpt_cat2)
ggsurvplot(
  fit2,                     # survfit object with calculated statistics.
  risk.table = TRUE,       # show risk table.
  pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE,         # show confidence intervals for
  # point estimaes of survival curves.
  #xlim = c(18,90),        # present narrower X axis, but not affect
  # survival estimates.
  #break.time.by = 10,    # break X axis in time intervals by 500.
  ggtheme = theme_classic2(), # customize plot and risk table with a theme.
  risk.table.y.text.col = T, # colour risk table text annotations.
  risk.table.y.text = FALSE # show bars instead of names in text annotations
  # in legend of risk table
)
### End OS data and cutpoint analysis ###

