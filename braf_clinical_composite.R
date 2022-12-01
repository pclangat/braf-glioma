## Pinky Langat 2020
## This code is designed to run after braf_oncoprint.R
## Has supporting analyses
#install.packages("logisticPCA")
library(logisticPCA)

library(tidyverse)
library(tidytidbits)
library(survivalAnalysis) 


### Composite figure  ###
#combined_classifier <- rbind((classifier_subsamps[, c("membership", "braf_class", "age_num", "subtype_path", "grade")]), jhh_os_data[, c("membership", "braf_class", "age_num", "subtype_path", "grade")])
#combined_classifier$subtype_path <- factor(combined_classifier$subtype_path, levels = c("GBM", "Astro", "Oligo", "PA", "PXA", "Other"))
#combined_classifier$grade <- factor(combined_classifier$grade, levels = c("1", "Low", "2", "3", "4"))

composite_df <- adult.braf.mastersheet

#composite_df$Age[combined_classifier$Age == "17.9"] <- 9.0

## Age by BRAF class density plots
#ggplot(data=combined_classifier, aes(x=age_num, group=braf_class)) +
ggplot(data=composite_df, aes(x=Age, group=BRAF_Class)) +
  #geom_density(fill="#69b3a2") +
  geom_histogram(fill="#69b3a2", binwidth = 2) +
  scale_x_continuous(breaks = seq(10, 90, by = 10))
ggplot(data=composite_df, aes(x=Age, group=BRAF_Class)) +
  geom_density(fill="#69b3a2") +
  #geom_histogram(fill="#69b3a2", binwidth = 5) +
  #scale_x_continuous(breaks = seq(10, 90, by = 5)) +
  #theme_ipsum() +
  theme_minimal() +
  facet_wrap(~BRAF_Class) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks.x=element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
table(composite_df$BRAF_Class)
plot((composite_df$Age[composite_df$BRAF_Class == "Gain"]))

ggplot(data=composite_df, aes(x=Age, group=BRAF_Class)) +
  #geom_density(fill="#69b3a2") +
  geom_histogram(fill="#69b3a2", binwidth = 2) +
  scale_x_continuous(breaks = seq(10, 90, by = 10))

## Histology pie chart
histopath_df <- data.frame(table(composite_df$Subtype, composite_df$BRAF_Class))
names(histopath_df) <- c("Subtype", "BRAF_Class", "value")
histopath_df
histopath_df %>%
  ggplot(aes(x="", y=value, fill=Subtype, group=BRAF_Class)) +
  geom_bar(stat="identity", width=1, position = 'fill') +
  coord_polar("y", start=0) +
  facet_wrap(~BRAF_Class) +
  theme_void() + 
  #scale_fill_manual(breaks = c("GBM-IDHwt" = "black", "Oligo" = "lightsalmon", "HGG" = "red", "Astro-IDHmut" = "green", "LGG" = "yellow", "Other" = "gray", "PA" = "skyblue", "Other High Grade" = "red", "PXA" = "violet", "Other glioma"= "purple")
  #scale_fill_manual(values=c("black", "green", "yellow", "skyblue", "violet", "gray"))
  scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "red", "green", "yellow", "lightsalmon", "skyblue", "violet"))
ano_col
###

## Grades stacked bar plot 
grades_df <- data.frame(table(composite_df$WHO_Grade, composite_df$BRAF_Class))
names(grades_df) <- c("Grade", "BRAF_Class", "value")
grades_df
grades_df %>%
  ggplot(aes(x="", y=value, fill=Grade)) +
  geom_bar(stat = 'identity', position = position_fill(reverse = TRUE)) +
  coord_flip() +
  facet_wrap(~BRAF_Class) +
  theme_minimal() +
  scale_fill_manual(values=c(brewer.pal(9, "YlGnBu")[c(3,5,7,9)])) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

## Summary BRAF mutations pie chart
table(composite_df$BRAF_Class)
data.frame(table(composite_df$BRAF_Class)) %>%
  ggplot(aes(x="", y=Freq, fill=Var1)) +
  geom_bar(stat="identity", width=1, position = 'fill') +
  coord_polar("y", start=0) +
  theme_void() +
  scale_fill_manual(values=c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray50", "gray"))
#scale_fill_manual(values=c(brewer.pal(6, "Spectral")[c(6,5,4,2)],"#fee0d2", "gray"))
#scale_fill_manual(values=c(brewer.pal(6, "PRGn")[c(6,5,4,2,3)], "gray"))
#scale_fill_manual(values=c(brewer.pal(6, "RdYlBu")[c(6,5,4,2,3)], "gray"))
#scale_fill_manual(values=c(brewer.pal(6, "RdBu")[c(6,5,4,2,3)], "gray"))
#scale_fill_manual(values=c(brewer.pal(6, "BrBG")[c(6,5,4,2,3)], "gray"))
#scale_fill_manual(values=c("#799B85","#8BBE9B", "#ACDACB", "#F7B9B2", "#B9D1E2", "#F9F9D2"))
#scale_fill_manual(values=c("#00441b", "#238b45", "#66c2a4", "#fb8072", "#80b1d3", "#ffffb3"))
#scale_fill_manual(values=c("#FF6131", "#FFB341", "#D7D452", "#4D8CE8","#9BD478", "gray"))
#scale_fill_manual(values=c("#4D8CE8", "#6ABFB6", "#9BD478", "#D7D452","#FF6131", "#D9D9D9"))

#braf_col <- c("Class I" = "#4D8CE8", "Class II" = "#6ABFB6", "Class III" = "#9BD478", "Fusion" = "#D7D452", "Rearrangement" = "#FFB341", "Gain" = "#FF6131", "Other" = "#D9D9D9")

brewer.pal(6, "Spectral")[c(6,5,4,2,1)]

### Stacked bar plot summaries
## BRAF Class
stack_df <- all.braf.mastersheet
stack_df$Patient_Cohort <- ifelse(stack_df$Age_Interval == "<18", "Pediatric", "Adult")
stack_df <- data.frame(table(stack_df$Patient_Cohort, stack_df$BRAF_Class))
names(stack_df) <- c("Patient_Cohort", "BRAF_Class", "frequency")
stack_df
ggplot(stack_df, aes(fill=BRAF_Class, y=frequency, x=Patient_Cohort)) + 
  #geom_bar(position="stack", stat="identity", color="black") +
  geom_bar(position=position_fill(reverse = TRUE), stat="identity", color="white") +
  theme_void() +
  scale_fill_manual(values=c(brewer.pal(6, "Spectral")[c(6,5,4,2,1)], "gray50", "gray"))

## Subtype stack
stack_df <- all.braf.mastersheet
stack_df$Patient_Cohort <- ifelse(stack_df$Age_Interval == "<18", "Pediatric", "Adult")
stack_df <- data.frame(table(stack_df$Patient_Cohort, stack_df$Subtype))
names(stack_df) <- c("Patient_Cohort", "Subtype", "frequency")
stack_df$Subtype <- factor(stack_df$Subtype, levels = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA"))
stack_df$Subtype
ggplot(stack_df, aes(fill=Subtype, y=frequency, x=Patient_Cohort)) + 
  #geom_bar(position="stack", stat="identity", color="black") +
  geom_bar(position=position_fill(reverse = TRUE), stat="identity", color="white") +
  theme_void() +
  #scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "red", "green", "yellow", "lightsalmon", "skyblue", "violet"))
  #scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "red", "chartreuse3", "gold", "lightsalmon", "skyblue", "violet"))
  #scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "#fc8d62", "#a6d854", "#ffd92f", "#e5c494", "skyblue", "violet"))
  #scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "red", "#a6d854", "#ffd92f", "lightsalmon", "skyblue", "violet"))
  scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "red", "#a6d854", "#ffed6f", "lightsalmon", "skyblue", "violet"))


## Grades stack
stack_df <- all.braf.mastersheet
stack_df$Patient_Cohort <- ifelse(stack_df$Age_Interval == "<18", "Pediatric", "Adult")
stack_df <- data.frame(table(stack_df$Patient_Cohort, stack_df$WHO_Grade))
names(stack_df) <- c("Patient_Cohort", "WHO_Grade", "frequency")
stack_df
ggplot(stack_df, aes(fill=WHO_Grade, y=frequency, x=Patient_Cohort)) + 
  #geom_bar(position="stack", stat="identity", color="black") +
  geom_bar(position="fill", stat="identity", color="white") +
  theme_void() +
  scale_fill_manual(values=c(brewer.pal(9, "YlGnBu")[c(3,5,7,9)]))

## Gender stack
stack_df <- all.braf.mastersheet
stack_df$Patient_Cohort <- ifelse(stack_df$Age_Interval == "<18", "Pediatric", "Adult")
stack_df <- data.frame(table(stack_df$Patient_Cohort, stack_df$Gender))
names(stack_df) <- c("Patient_Cohort", "Gender", "frequency")
stack_df
ggplot(stack_df, aes(fill=Gender, y=frequency, x=Patient_Cohort)) + 
  #geom_bar(position="stack", stat="identity", color="black") +
  geom_bar(position="fill", stat="identity", color="white") +
  theme_void() +
  scale_fill_manual(values=c(brewer.pal(3, "PiYG")[c(1,3)]))

## Fusion Partners stack
fusion_sankey_df
stack_df <- fusion_sankey_df
stack_df$Patient_Cohort <- ifelse(stack_df$BRAF == "BRAF Pediatric", "Pediatric", "Adult")
stack_df
#stack_df <- data.frame(table(stack_df$Patient_Cohort, stack_df$WHO_Grade))
#names(stack_df) <- c("Patient_Cohort", "WHO_Grade", "frequency")
stack_df$Fusion_Partner <- factor(stack_df$Fusion_Partner, levels = c("KIAA1549", "FAM131B", "AGK", "GNAI1", "UBE2H", "BCAS1", "CCDC6", "GIT2", "PTPRZ11"))
stack_df
ggplot(stack_df, aes(fill=Fusion_Partner, y=Freq, x=Patient_Cohort)) + 
  #geom_bar(position="stack", stat="identity", color="black") +
  #geom_bar(position="fill", stat="identity", color="white") +
  geom_bar(position=position_fill(reverse = TRUE), stat="identity", color="white") +
  theme_void() +
  #theme(legend.position="bottom") +
  scale_fill_brewer(palette="Set3")
#scale_fill_manual(values=c(brewer.pal(9, "YlGnBu")[c(3,5,7,9)]))

# Fusion Partner Pie Chart
fusion_df <- data.frame(table(fusion_sankey_df$Fusion_Partner, fusion_sankey_df$BRAF))
fusion_sankey_df
fusion_df <- fusion_sankey_df[c("Fusion_Partner", "BRAF", "Freq")]
#names(fusion_df) <- c("Fusion_Partner", "BRAF", "value")
fusion_df
fusion_df %>%
  ggplot(aes(x="", y=Freq, fill=Fusion_Partner, group=BRAF)) +
  geom_bar(stat="identity", width=1, position = 'fill', color='white') +
  coord_polar("y", start=0) +
  facet_wrap(~BRAF) +
  theme_void() +
  theme(legend.position="bottom") +
  #geom_text(aes(y =fill, label = Fusion_Partner), color = "white", size=6) #+
  scale_fill_brewer(palette="Set3")
#scale_fill_manual((values=c(brewer.pal(9, "Set3"))))
#scale_fill_manual(breaks = c("GBM-IDHwt" = "black", "Oligo" = "lightsalmon", "HGG" = "red", "Astro-IDHmut" = "green", "LGG" = "yellow", "Other" = "gray", "PA" = "skyblue", "Other High Grade" = "red", "PXA" = "violet", "Other glioma"= "purple")
#scale_fill_manual(values=c("black", "green", "yellow", "skyblue", "violet", "gray"))
#scale_fill_manual(breaks = c("GBM, IDH-wt", "HGG", "Astro, IDH-mt", "LGG", "Oligo", "PA", "PXA" ), values=c("black", "red", "green", "yellow", "lightsalmon", "skyblue", "violet"))
ano_col
