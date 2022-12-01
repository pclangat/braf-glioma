library(ggpubr)
library(tidyverse)
library(data.table)

adult.data <- fread("censored_braf_adult_glioma_mastersheet.2021-07-09.csv")
peds.data <- fread("censored_braf_peds_glioma_mastersheet.2021-07-09.csv")


metadata <- read.csv("braf_adult_glioma_classifier_datasheet.braf_mut_classes_revised.2020-10-30.csv",
                     header = T, stringsAsFactors = F)

gene.list <- fread("geneList.csv")

not.covered.1 <- fread("NotCovered-1.csv")

all.dfci.genes <- setdiff(gene.list$Gene, not.covered.1$x)
all.dfci.genes <- all.dfci.genes[all.dfci.genes!=""]
## 287 genes 

panels.excluded <- c("JHU-50GP-V2", "MDA-46-V1", "MDA-50-V1")
individuals.excluded <- unique(adult.data$SAMPLE_ACCESSION_NBR[which(adult.data$PANEL_VERSION %in% panels.excluded)])
peds.excluded <- unique(peds.data$SAMPLE_ACCESSION_NBR[which(peds.data$PANEL_VERSION %in% panels.excluded)])
## none excluded 


genie.genes.trimmed <- fread("GENIE.geneList.csv")
genie.genes <- unique(genie.genes.trimmed$Hugo_Symbol)
genie.panels <- unique(genie.genes.trimmed$SEQ_ASSAY_ID)
genie.panels <- intersect(genie.panels, unique(metadata$PANEL_VERSION))
genie.panels <- setdiff(genie.panels, panels.excluded)
for(g in genie.panels){
  genie.genes <- intersect(genie.genes, genie.genes.trimmed$Hugo_Symbol[which(genie.genes.trimmed$SEQ_ASSAY_ID == g)])
}

## 178 genes 



############ Load mutation data

data <- read.table("alterations_braf_samples.txt", header = T, quote = "", stringsAsFactors = F, sep = "\t")
names(data) <- c("sample", "gene", "classification")

data$sample <- gsub('\\"', '', data$sample)
data$gene <- gsub('\\"', '', data$gene)
data$classification <- gsub('\\"', '', data$classification)

data <- data[!duplicated(data),]

table(data$classification)

data <- subset(data, sample %in% union(adult.data$SAMPLE_ACCESSION_NBR,
                                       peds.data$SAMPLE_ACCESSION_NBR))





# neuropathpanel <- scan("NeuropathNGS-panel.txt", what = character())
jhuv2 <- scan("SolidTumorPanel-II_GeneList_v2.0.txt", what = character())

genes.used <- intersect(all.dfci.genes, genie.genes)
genes.used <- intersect(genes.used, jhuv2)
## 136 genes 
# saveRDS(genes.used, file = "IntersectGeneList.RDS"

jhu.data <- read.csv("JHH-genomics.csv", header = T, stringsAsFactors = F)
names(jhu.data) <- c("sample", "panel", "gene")

jhu.individuals.excluded <- jhu.data$sample[which(grepl("NGS", jhu.data$panel))]
jhu.data <- jhu.data[!(jhu.data$sample %in% jhu.individuals.excluded),]

jhu.data <- jhu.data[which(jhu.data$sample != 13),]



combined.data <- subset(data[,c("sample", "gene")], gene %in% genes.used)
combined.data$cohort <- unlist(lapply(combined.data$sample, function(x) strsplit(x,"-")[[1]][[1]]))
tmp <- subset(jhu.data[,c("sample", "gene")], gene %in% genes.used)
tmp$sample <- paste0("JHH-", tmp$sample)
tmp$cohort <- "JHH"
combined.data <- rbind(combined.data, tmp)

all.patients <- unique(combined.data$sample)
all.patients.cohort <- unlist(lapply(all.patients, function(x) strsplit(x, "-")[[1]][[1]]))
## 261 patients, 13 from JHU



######## Merging metadata 

adult.metadata <- adult.data[,c("SAMPLE_ACCESSION_NBR",
                                "Cohort",
                                "Gender",
                                "Age_Interval",
                                "KPS",
                                "Tumor_Laterality",
                                "WHO_Grade",
                                "BRAF_Class",
                                "Subtype_Path",
                                "Subtype")]

peds.metadata <- peds.data[,c("SAMPLE_ACCESSION_NBR",
                                "Cohort",
                                "Gender",
                                "Age_Interval",
                                "KPS",
                                "Tumor_Laterality",
                                "WHO_Grade",
                                "BRAF_Class",
                                "Subtype_Path",
                                "Subtype")]






combined.metadata <- rbind(adult.metadata, peds.metadata)
combined.metadata$Patient <- combined.metadata$SAMPLE_ACCESSION_NBR

all.patients.gender <- combined.metadata$Gender[match(all.patients, combined.metadata$Patient)]
all.patients.gender[which(is.na(all.patients.gender))] <- "Unknown"

all.patients.subtypepath <- combined.metadata$Subtype_Path[match(all.patients, combined.metadata$Patient)]
all.patients.subtypepath[which(is.na(all.patients.subtypepath))] <- "Unknown"

all.patients.subtype <- combined.metadata$Subtype[match(all.patients, combined.metadata$Patient)]
all.patients.subtype[which(is.na(all.patients.subtype))] <- "Unknown"


all.patients.grade <- combined.metadata$WHO_Grade[match(all.patients, combined.metadata$Patient)]
all.patients.grade[which(is.na(all.patients.grade))] <- "Unknown"

all.patients.brafclass <- combined.metadata$BRAF_Class[match(all.patients, combined.metadata$Patient)]
all.patients.brafclass[which(is.na(all.patients.brafclass))] <- "Unknown"

all.patients.laterality <- combined.metadata$Tumor_Laterality[match(all.patients, combined.metadata$Patient)]
all.patients.laterality[which(is.na(all.patients.laterality))] <- "Unknown"

all.patients.agegroup <- combined.metadata$Age_Interval[match(all.patients, combined.metadata$Patient)]
all.patients.agegroup[which(is.na(all.patients.agegroup))] <- "Unknown"






df <- data.frame(matrix(0, ncol = length(genes.used), nrow = length(all.patients)))
colnames(df) <- genes.used
rownames(df) <- all.patients

for(i in 1:nrow(combined.data)){
  df[combined.data$sample[[i]], combined.data$gene[[i]]] <- 1
}

saveRDS(df, file = "CombinedMutationData.RDS")

df <- readRDS("CombinedMutationData.RDS")


hypermutated <- c("TCGA-DU-6392","TCGA-14-1396", "TCGA-06-5416", "TCGA-19-5956",
                  "GENIE-MSK-P-0003529-T01-IM5", "GENIE-MSK-P-0009499-T01-IM5",
                  "BL-15-N26622", "GENIE-MSK-P-0023123-T01-IM6", "GENIE-MSK-P-0013293-T02-IM6",
                  "BL-15-N38303", "BL-16-A40851", "TCGA-16-1460", "GENIE-MSK-P-0019556-T01-IM6 BL-16-R13109",
                  "BL-15-T48446", "GENIE-MSK-P-0000500-T01-IM3 BL-17-J18600", "GENIE-MSK-P-0017675-T01-IM5",
                  "BL-14-J36353",  "BL-14-A45755", "GENIE-UCSF-1612-1771T")
all.patients.hypermutated <- ifelse(all.patients %in% hypermutated, "Yes", "No")


#### Sept 11, 2021: remove hypermutated
to.remove <- union(which(rownames(df) %in% hypermutated), which(all.patients.subtype == "Oligo"))
df.removehyper <- df[-to.remove,]
all.patients.gender <- all.patients.gender[-to.remove]
all.patients.subtype <- all.patients.subtype[-to.remove]
all.patients.grade <- all.patients.grade[-to.remove]
all.patients.brafclass <- all.patients.brafclass[-to.remove]
all.patients.agegroup <- all.patients.agegroup[-to.remove]






###########################################
## Chisq analysis: correlation gene ~ features
###########################################

library(vcd)

stat.data <- as.data.frame(df.removehyper)

stat.data <- stat.data %>%
  mutate_all(as.factor)


gene.names <- names(stat.data)

stat.data$Gender <- all.patients.gender
stat.data$Grade <- all.patients.grade
stat.data$BRAFClass <- all.patients.brafclass
stat.data$Age <- all.patients.agegroup
stat.data$Subtype <- all.patients.subtype
# stat.data$Hypermutated <- all.patients.hypermutated


variable.names <- c("Gender", "Grade", "BRAFClass", "Age", "Subtype")


results.data <- data.frame()

for(i in 1:length(gene.names)){
  for(j in 1:length(variable.names)){
    tab <- xtabs(as.formula(paste0("~ ", gene.names[[i]], " + ", variable.names[[j]])),
                 data = stat.data)
    tmp <- summary(assocstats(tab))
    
    
    this.result <- data.frame(Gene = gene.names[[i]],
                              Variable = variable.names[[j]],
                              Chisq = tmp$summary$statistic,
                              Pvalue = tmp$summary$p.value)
    
    results.data <- rbind(results.data, this.result)
    
  }
  if(i %% 10 == 0) print(i)
}

fwrite(results.data, file = "Correlation-CremerVStat.csv",
       col.names = T, row.names = F, quote = F)


############################################
## post hoc test
############################################

library(chisq.posthoc.test)


results.data <- data.frame()

for(i in 1:length(gene.names)){
  for(j in 1:length(variable.names)){
    
    tab <- xtabs(as.formula(paste0("~ ", gene.names[[i]], " + ", variable.names[[j]])),
                 data = stat.data)
    if(min(dim(tab)) == 1) next
    
    tmp <- chisq.posthoc.test(as.matrix(tab))
    
    tmp <- tmp[which(tmp$Dimension == 1),]
    
    # tmp <- tmp[grepl("p values", tmp[,2]),]
    rownames(tmp) <- tmp[,c("Value")]
    tmp <- tmp[,-c(1, 2)]
    tmp <- as.data.frame(t(tmp))
    tmp$Residuals <- as.numeric(tmp$Residuals)
    
    tmp$Direction <- sign(tmp$Residuals)
    
    
    this.result <- data.frame(Gene = gene.names[[i]],
                              Variable = variable.names[[j]],
                              Category = rownames(tmp),
                              Residuals = tmp$Residuals,
                              Pvalue = tmp$`p values`)
    
    results.data <- rbind(results.data, this.result)
    
  }
  if(i %% 10 == 0) print(i)
}


results.data$p_numeric <- as.numeric(gsub("\\*", "", unlist(lapply(results.data$Pvalue, toString))))

fwrite(results.data, file = "Correlation-Chisq-posthoc-removehypermutated.txt",
       col.names = T, row.names = F, quote = F, sep = "\t")


### plotting:
library(pheatmap)
results.data <- fread("Correlation-Chisq-posthoc-removehypermutated.txt")


plot.data <- results.data %>%
  filter(Variable == "Subtype") %>%
  select(c("Gene", "Category", "p_numeric")) 

plot.data$p_adj <- p.adjust(plot.data$p_numeric, method = "fdr")
plot.data$p_adj_masked <- ifelse(plot.data$p_adj < 0.2, plot.data$p_adj, NA)

plot.data <- plot.data %>%
  select(c("Gene", "Category", "p_adj_masked"))  %>%
  pivot_wider(names_from = Category, values_from = p_adj_masked)



genes <- plot.data$Gene
plot.data <- as.matrix(-log10(plot.data[,2:ncol(plot.data)]+ 1e-4))
rownames(plot.data) <- genes

plot.data <- plot.data[rowSums(is.na(plot.data)) < 6,]



p <- pheatmap((plot.data),
              fontsize = 15,
              # annotation_col = patient.annotation,
              # clustering_distance_cols = "euclidean",
              cluster_rows = F,
              # clustering_distance_rows = "correlation",
              cluster_cols = F,
              show_colnames = T) 

ggsave(p, filename = "Heatmap_Subtype_posthoc.png", width = 10, height = 15)


#############################

plot.data <- results.data %>%
  filter(Variable == "BRAFClass") %>%
  select(c("Gene", "Category", "p_numeric")) 

plot.data$p_adj <- p.adjust(plot.data$p_numeric, method = "fdr")
plot.data$p_adj_masked <- ifelse(plot.data$p_adj < 0.2, plot.data$p_adj, NA)

plot.data <- plot.data %>%
  select(c("Gene", "Category", "p_adj_masked"))  %>%
  pivot_wider(names_from = Category, values_from = p_adj_masked)



genes <- plot.data$Gene
plot.data <- as.matrix(-log10(plot.data[,2:ncol(plot.data)]+ 1e-4))
rownames(plot.data) <- genes

plot.data <- plot.data[rowSums(is.na(plot.data)) < ncol(plot.data),]



p <- pheatmap((plot.data),
              fontsize = 15,
              # annotation_col = patient.annotation,
              # clustering_distance_cols = "euclidean",
              cluster_rows = F,
              # clustering_distance_rows = "correlation",
              cluster_cols = F,
              show_colnames = T) 

ggsave(p, filename = "Heatmap_BRAFClass_posthoc.png", width = 10, height = 15)



#############################

plot.data <- results.data %>%
  filter(Variable == "Grade") %>%
  select(c("Gene", "Category", "p_numeric")) 

plot.data$p_adj <- p.adjust(plot.data$p_numeric, method = "fdr")
plot.data$p_adj_masked <- ifelse(plot.data$p_adj < 0.2, plot.data$p_adj, NA)

plot.data <- plot.data %>%
  select(c("Gene", "Category", "p_adj_masked"))  %>%
  pivot_wider(names_from = Category, values_from = p_adj_masked)



genes <- plot.data$Gene
plot.data <- as.matrix(-log10(plot.data[,2:ncol(plot.data)]+ 1e-4))
rownames(plot.data) <- genes

plot.data <- plot.data[rowSums(is.na(plot.data)) < ncol(plot.data),]



p <- pheatmap((plot.data),
              fontsize = 15,
              # annotation_col = patient.annotation,
              # clustering_distance_cols = "euclidean",
              cluster_rows = F,
              # clustering_distance_rows = "correlation",
              cluster_cols = F,
              show_colnames = T,
              angle_col = 0) 

ggsave(p, filename = "Heatmap_Grade_posthoc.png", width = 6, height = 8)



### We stopped here 














# 
# 
# M <- as.table(rbind(c(762, 327, 468), c(484, 239, 477)))
# dimnames(M) <- list(gender = c("F", "M"),
#                     party = c("Democrat","Independent", "Republican"))
# chisq.test(M)
# #> 
# #>  Pearson's Chi-squared test
# #> 
# #> data:  M
# #> X-squared = 30.07, df = 2, p-value = 2.954e-07
# # Standarized residuals
# # As a form of post hoc analysis the standarized residuals can be analysed. A rule of thumb is that standarized residuals of above two show significance.
# 
# chisq.results <- chisq.test(M)
# chisq.results$stdres
# #>       party
# #> gender   Democrat Independent Republican
# #>      F  4.5020535   0.6994517 -5.3159455
# #>      M -4.5020535  -0.6994517  5.3159455
# # Post Hoc Analysis
# # However, the above two rule is a rule of thumb. These standarized residuals can be used to calculate p-values, which is what this package is designed for as shown in the following example.
# 
# chisq.posthoc.test(M,
#                    method = "bonferroni")
# 
# 
# 
































# 
# 
# data("Arthritis")
# tab <- xtabs(~Improved + Treatment, data = Arthritis)
# tmp <- summary(assocstats(tab))





###############################
## Addendum: logistic PCA
###############################

library(logisticPCA)

logpca_cv = cv.lpca(df, ks = 2, ms = 1:10)
plot(logpca_cv)

# optimal m = 3

## remove columns with all 0 or 1

calSD <- df %>%
  summarise_all(sd)
df.cleaned <- df[,which(calSD > 0)]


clogpca_model = convexLogisticPCA(df.cleaned, k = 2, m = which.min(logpca_cv))



### plotting

# party = rownames(house_votes84)
# plot(logsvd_model, type = "scores") + geom_point(aes(colour = party)) + 
#   ggtitle("Exponential Family PCA") + scale_colour_manual(values = c("blue", "red"))


plot(clogpca_model, type = "scores") + geom_point(aes(colour = party)) + 
  ggtitle("Convex Logistic PCA") +
  scale_colour_manual(values = c("blue", "red"))



logPCA.data <- data.frame(PC1 = clogpca_model$PCs[,1],
                          PC2 = clogpca_model$PCs[,2],
                       Cohort = all.patients.cohort,
                       Gender = all.patients.gender,
                       Grade = all.patients.grade,
                       SubtypePath = all.patients.subtypepath,
                       SubtypeMolecular = all.patients.subtypemole,
                       Hypermutated = all.patients.hypermutated,
                       Laterality = all.patients.laterality,
                       BRAF_Class = all.patients.brafclass,
                       Age = all.patients.agegroup)



plist <- lapply(names(logPCA.data)[3:ncol(logPCA.data)], function(var.name){
  ggscatter(logPCA.data, x = 'PC1', y = "PC2",
            #alpha = 0.3,
            #size = 0.6,
            color = var.name, #shape = "Category",
            #title = plot.title,
            palette = "jco",
            # palette = c("#e41a1c", "#377eb8", "#4daf4a"),
            title = paste0(var.name)) +
    # geom_label_repel(aes(label = label,
    #                      fill = label), color = 'white',
    #                  size = 3.5) +
    theme_pubclean() +
    theme(axis.title = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.text = element_text(size = 12))
})





p <- ggarrange(plotlist = plist, ncol = 3, nrow = 3)


ggsave(p, filename = "Genomics-logPCA-WithPeds.png", 
       width = 15, height = 15)

















##########################
# Analysis 1: PCA
##########################

pcaData <- prcomp((as.matrix(df)))


pca.data <- data.frame(PC1 = pcaData$x[,"PC1"],
                       PC2 = pcaData$x[,"PC2"],
                       Cohort = all.patients.cohort,
                       Gender = all.patients.gender,
                       Grade = all.patients.grade,
                       SubtypePath = all.patients.subtypepath,
                       SubtypeMolecular = all.patients.subtypemole,
                       Hypermutated = all.patients.hypermutated,
                       Laterality = all.patients.laterality,
                       BRAF_Class = all.patients.brafclass,
                       Age = all.patients.agegroup)


# saveRDS(pca.data, file = paste0("PCAData-", length(genes.used), "Genes-", length(all.patients), "Patients.RDS"))


percentVar <- round(100 * (pcaData$sdev^2/sum(pcaData$sdev^2)))



plist <- lapply(names(pca.data)[3:ncol(pca.data)], function(var.name){
  ggscatter(pca.data, x = 'PC1', y = "PC2",
            #alpha = 0.3,
            #size = 0.6,
            color = var.name, #shape = "Category",
            #title = plot.title,
            palette = "jco",
            # palette = c("#e41a1c", "#377eb8", "#4daf4a"),
            title = paste0(var.name),
            xlab = paste0("PC1: ", percentVar[1], "% variance"),
            ylab = paste0("PC2: ", percentVar[2], "% variance")) +
    # geom_label_repel(aes(label = label,
    #                      fill = label), color = 'white',
    #                  size = 3.5) +
    theme_pubclean() +
    theme(axis.title = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.text = element_text(size = 12))
})


p <- ggarrange(plotlist = plist, ncol = 3, nrow = 3)
ggsave(p, filename = "Genomics-PCA-WithPeds.png", 
       width = 15, height = 15)












############################
# Analysis 2: UMAP
############################

library(umap)

tmp <- umap(df)



plot.data <- data.frame(UMAP1 = tmp$layout[,1],
                       UMAP2 = tmp$layout[,2],
                       Cohort = all.patients.cohort,
                       Gender = all.patients.gender,
                       Grade = all.patients.grade,
                       SubtypePath = all.patients.subtypepath,
                       SubtypeMolecular = all.patients.subtypemole,
                       Hypermutated = all.patients.hypermutated,
                       Laterality = all.patients.laterality,
                       BRAF_Class = all.patients.brafclass,
                       Age = all.patients.agegroup)




plist <- lapply(names(plot.data)[3:ncol(plot.data)], function(var.name){
  ggscatter(plot.data, x = 'UMAP1', y = "UMAP2",
            #alpha = 0.3,
            #size = 0.6,
            color = var.name, #shape = "Category",
            #title = plot.title,
            palette = "jco",
            # palette = c("#e41a1c", "#377eb8", "#4daf4a"),
            title = paste0(var.name)) +
    # geom_label_repel(aes(label = label,
    #                      fill = label), color = 'white',
    #                  size = 3.5) +
    theme_pubclean() +
    theme(axis.title = element_text(size = 13, face = "bold"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          legend.position = "right",
          axis.text = element_text(size = 12))
})


p <- ggarrange(plotlist = plist, ncol = 3, nrow = 3)
ggsave(p, filename = "Genomics-UMAP-WithPeds.png", 
       width = 15, height = 15)



