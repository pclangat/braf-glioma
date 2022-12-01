## Original code:
## Noah Greenwald
## 9/17/15
## Revisions by:
## Pinky Langat 2022

## Helpful R functions for manipulating MAF output files for data analysis. 
#source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering/ExacFilter.R")
#source("C:/p/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//ESPFilter.R")
#source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//PoNFilter.R")
#source("C:/Users/Noah/OneDrive/Work/Coding/R/Scripts/Meningioma/Filtering//CallCompare.R")

## Deprecated
#source("http://www.bioconductor.org/biocLite.R")
#biocLite(c("rtracklayer", "AnnotationHub", "Rsamtools"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#BiocManager::install(c("rtracklayer", "AnnotationHub", "Rsamtools"))

#install.packages("dplyr")
#install.packages("ggplot2")
#install.packages("devtools")
#install.packages("R.utils")
#install.packages("kohonen")
#install.packages("nycflights13")
#install.packages("randomcoloR")
#install.packages("NMF")
#install.packages("gridExtra")
#install.packages("cowplot")
#install.packages("vegan")
#install.packages("limma")


## Nice packages
#library(grDevices)
library(RColorBrewer)
library(dplyr)
#library(rtracklayer)
library(ggplot2)
#library(curl)
library(R.utils)
#library(Biobase)
#library(kohonen)
#library(nycflights13)
library(reshape2)
library(randomcoloR)
library(NMF)
library(gridExtra)
#library(BiocManager)
library(vegan)


## ggplot settings
rameen_theme <- theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
                      axis.title.y = element_blank(), axis.line.x = element_line(size = 1.25, color = "black"), axis.line.y = element_line(size = 1.25, color = "black"), 
                      axis.text.x = element_text(angle = 90, hjust = 1))
rameen_theme2 <- theme(panel.background = element_blank(), axis.ticks.y = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
                      axis.title.y = element_blank(), axis.line.x = element_line(size = 1.25, color = "black"), axis.line.y = element_line(size = 1.25, color = "black"), 
                      axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = 6))

FilterMaf <- function(maf, list, column.name, keep = TRUE){
  # Takes a maf, and returns a maf with only those rows who have a value in given column that matches an entry in the list. 
  
  # Args:
  #     maf: a maf
  #     list: list of column values to keep
  #     column.name: column in the maf which will be checked
  #     keep: optional, logical value. When false, discards selected rows
  
  # Returns: A maf matching input criteria
  
  ## Error Checking
  if (is.null(maf[[column.name]])){
    stop("Invalid column name")
  }
  
  attendance <- list %in% maf[[column.name]]
  if (sum(attendance) < length(attendance)){
    warning("List element ", list[which(!attendance)], " not found")
  }
  
  if (keep == TRUE){
    idx <- which(maf[[column.name]] %in% list)
    new.maf <- maf[idx, ]
    return (new.maf)
  }else{
    idx <- which(!(maf[[column.name]] %in% list))
    new.maf <- maf[idx, ]
    return (new.maf)
  }
}

PlotMaf <- function(maf, column.name, percent = 100, lower.margin = 1, title = ""){
    ## Takes a Maf file and produces a bar plot 
    
    # Args: 
    #     Maf: a maf file in array format. 
    #     column.name: a string with the name of the column that the barplot will be produced from. 
    #     lower.margin: optional, a multiplier that adjusts the marign at the bottom of the graph. 
    #     title: optional, gives the graph a title
    #     percent: optional, allows for plotting of only top x% of hits ranked by frequency
    
    ## Error checking
    if (is.null(maf[[column.name]])){
        stop("Invalid column name")
    }
    
    if (percent <= 0 | percent > 100){
        stop("Learn how to calculate a percent")
    }
    data <- maf[, column.name]
    tbl <- table(data)
    tbl <- sort(tbl, decreasing = T)
    tbl <- tbl[1:floor(length(tbl)*percent/100)]
    par(mar= c(5*lower.margin,4,4,2) + .1)
    size = (length(unique(data)))
    barplot(tbl, las = 2, main = title)
    
    #   if (size > 60){
    #       x = ceiling(size / 60)
    #       for (i in 1:x){
    #         end
    #         barplot(tbl, las = 2, main = c(title," ", i, " of ", x))
    #       }
    #   }
    #   else{
    #   barplot(tbl, las = 2, main = title)
    #   }
}

AverageMaf <- function(maf, list, list.column, int.column){
  # Takes a maf and produces a graph with the average value of a number of categories within the maf
  
  # Args:
  #     maf: a maf file in data.table format
  #     list: a list of the groupings the function will use to calcuate averages
  #     list.column: column of the maf which the list is drawn from
  #     int.column: column whose values will be summed for each entry
  
  # Returns: a list of averages corresponding to each item in input list, values of which are found in int.column
  
  ## Error Checking
  if (is.null(maf[[list.column]]) | is.null(maf[[int.column]])){
    stop("Invalid column name")
  }
  
  attendance <- list %in% maf[[list.column]]
  if (sum(attendance) < length(attendance)){
    warning("List element ", list[which(!attendance)], " not found")
  }
  
  ## Creates a vector to store average for each entry
  values <- c()
  
  ## loop through list, computing average for each entry and storing in vector
  for(i in 1:length(list)){
    temp.frame <- maf[which(maf[[list.column]] == list[i]), ]
    x <-  mean(temp.frame[[int.column]])
    values[i] <- x
  }
  return(values)
}


ReccurentMaf <- function(maf, column.name, threshold = 1){
  # Takes a maf file, and returns a filtered maf with only those entries in designated column with
  # a frequency greater than the threshold value
  
  # Args:
  #     maf: a maf file
  #     column.name: the designated column in maf file
  #     threshold: cutoff value for # of times column entry must be seen to be included
  
  # Returns: a filtered maf
  
  ## Error Checking
  if (is.null(maf[[column.name]])){
    stop("Invalid column name")
  }
  
  tbl <- table(maf[[column.name]])
  approved.list <- c()
  for (i in 1:length(tbl)){
    if (tbl[[i]] > threshold){
      approved.list <- c(approved.list, names(tbl[i]))
    }
  }
  FilterMaf(maf, approved.list, column.name)
}

RatioMaf <- function (maf, column.name, list1, list2){
  ## Takes a maf file, and returns the number of occurences of variables found in list 1 and 2
  
  ## Args:
  #       maf: a maf file
  #       column.name: the column in maf file that list elements are drawn from
  #       list1: list of possible values for given column
  #       list2: second list of possible values for given column
  
  # Returns: summary of occurences of selected elements in each list. 
  
  
  temp1 <- FilterMaf(maf, list1, column.name)
  temp2 <- FilterMaf(maf, list2, column.name)
  val1 <- nrow(temp1)
  val2 <- nrow(temp2)
  valtot <- nrow(maf)
  paste("There are a total of ", valtot, "mutations in your maf. Of these, ", val1, 
        " were in the first subset, and ", val2, " were in the second subset.")
}

PerSampleMaf <- function(maf, column.name, identifier.column = "Tumor_Sample_Barcode", list = maf[[column.name]]){
  ## Takes a maf, and returns a filtered maf that limits occurence to one per sample
  
  ## Args:
  #       maf: a maf file
  #       column.name: name of the column to be checked for uniqueness
  #       identifier.column: name of the column with sample IDs            
  #       list: values that will be kept in filtered maf. Default is all values in column
  
  # Returns: a maf with only those rows which contain an element from list in selected column, 
  # limited to once per sample. 
  filtered <- FilterMaf(maf, list, column.name, TRUE)
  cols <- c(identifier.column, column.name)
  dups <- filtered[, cols]
  idx <- duplicated(dups)
  filtered <- filtered[!idx, ]
}            

CompareMutsMaf <- function(maf, cutoff = 0, title = ""){
  ## Takes a maf and a minimum value for # of mutations, and returns a comparison chart of total 
  ## number of mutations in given genes, as well as # of samples with given mutation
  
  ## Args:
  #       maf: a maf file
  #       cutoff: minimum number of times a mutation is seen in dataset in order to be kept
  
  # Returns: a side-by-side barplot comparing frequencies of total vs per sample mutations per gene
  
  filtered <- ReccurentMaf(maf, "Hugo_Symbol", cutoff)
  per.gene <- table(filtered[["Hugo_Symbol"]])
  per.sample.maf <- PerSampleMaf(filtered, "Hugo_Symbol")
  per.sample <- table(per.sample.maf[["Hugo_Symbol"]])
  data <- rbind(per.gene, per.sample)
  par(mar=c(10, 4, 4, 2) + .1)
  barplot(data, beside = TRUE, legend.text = c("Total mutations per gene", "Samples with gene mutated"), 
          las = 2, main = title )
}

## Compares different cutoffs for filtering


FilterCutoffMaf <- function(maf1, maf2, cutoff, cut1 = "Filter 1", cut2 = "Filter2",  title = ""){
  maf1 <- PerSampleMaf(maf1, "Hugo_Symbol")
  maf1 <- ReccurentMaf(maf1, "Hugo_Symbol", cutoff)
  maf2 <- PerSampleMaf(maf2, "Hugo_Symbol")
  maf2 <- ReccurentMaf(maf2, "Hugo_Symbol", cutoff)
  data1 <- table(maf1[["Hugo_Symbol"]])
  data2 <- table(maf2[["Hugo_Symbol"]])
  data <- EqualizeTable(data1, data2)
  barplot(data, beside = TRUE, legend.text = c(cut1, cut2),
          las = 2, main = title)
}


## Takes two tables of values, makes sure all values from each table appear in the other, and orders the table
## so that all entries are in same position. 

EqualizeTable <- function(table1, table2){
  
  ## Takes all entries not found in table 1, sets them to zero value
  table1.names <- names(table1)
  table2.names <- names(table2)
  free <- !(table2.names %in% table1.names)
  vals <- c(rep(0, sum(free)))
  names(vals) <- table2.names[free]
  table1 <- c(table1, vals)
  table1 <- table1[order(names(table1))]
  
  ## Does the same thing, vice aversa
  free <- !(table1.names %in% table2.names)
  vals <- c(rep(0, sum(free)))
  names(vals) <- table1.names[free]
  table2 <- c(table2, vals)
  table2 <- table2[order(names(table2))]

  table3 <- rbind(table1, table2)
  table3
}


CleanYourRoom <- function(input.data, naughty.list = c()){
## Takes a data frame or matrix, and removes rows with problemtic entries

## Args:
##      input.data: the data to be input
##      naughty.list: additional parameters to exclude

    for(i in 1:ncol(input.data)){
        bad.boys <- input.data[, i] %in% c("", "-", "N/A", "na", NA, naughty.list)
        input.data <- input.data[!bad.boys, ]
    }
    return (input.data)
}

TestStatistic <- function(vector, mu){
## takes a vector of values, determines if average is significantly different from input
    
## Args:
##      vector: vector of numeric values
##      mu: value to compare vector values to    
    deviation <- sd(vector)
    average <- mean(vector)
    denom <- deviation / sqrt(length(vector))
    return((mu - average) / denom)
}

is.wholenumber <- function(x, tol = .Machine$double.eps^0.5){
    ## checks if a number is an integer
    ## taken from Integer help page documentation examples
    abs(x - round(x)) < tol
}

PlotComut <- function(mut.maf1 = NULL, mut.maf2 = NULL, focal.cnv.maf = NULL, broad.cnv.maf = NULL, samples, input.samples, input.genes = "Hugo_Symbol", input.mutations = "Variant_classification", 
                      gene.cutoff = 1, file.path = NULL, unique.muts = NULL, col.vector = NULL, phenotypes = NULL, manual.order = NULL, fixed.order = NULL, return.matrix = FALSE, 
                      return.plot = FALSE, coverage = TRUE, y.axis.font.size = 8, legend.font.size = 8, dimensions = c(7, 7), title = "Comut Plot", pheno_top = FALSE){
## Takes in a maf file, and generates a comut plot of common mutations in each of the samples
## adapted from code at https://benchtobioinformatics.wordpress.com/2015/05/25/how-to-make-a-co-mutation-plot/
    

## Args:
    ## mut.maf1: a maf file with mutation calls
    ## mut.maf2: an optional maf file with mutation calls, all of which will sorted below the first and labeled tier 4
    ## focal.cnv.maf: an optional maf (long format) of CNV calls with sample identifier as first column, gene as second, and type of CNV as third
    ## broad.cnv.maf: an optional maf
    ## samples: a vector with all individuals to be plotted (don't take from MAF, could be absent due to no mutations)
    ## input.samples: a column wthin the maf which specifies samples
    ## input.genes: a column within the maf which specifies genes
    ## input.mutations: a column within the maf describing the type of mutation
    ## gene.cutoff: a threshold for deciding how many genes to include in comut. if beteen 0 and 1, treated as percentage of total. otherwise, number of genes
    ## file.path: optional path specifying where the comut will be saved
    ## unique.muts: an optional vector with all the unique mutation classes present in dataset. Useful for consistent colors across multiple comuts
    ## col.vector
    ## phenotypes: an optional dataframe with samples as the first column, and subsequent columns taken as additional information to plot in given order
    ## manual.order: a vector of variables names which will overide the default ordering of rows by frequency for alternate sorting
    ## fixed.order: a vector of sample identifiers which determines ordering of all samples
    ## return.matrix: if TRUE, returns the wide dataframe used to sort the comut to enable statistical testing
    ## title: a title for the plot
    
    
    ## Error Checking
    print("beginning")
    if(!missing(mut.maf1)){
        if (is.null(mut.maf1[[input.samples]]) |is.null(mut.maf1[[input.genes]]) | is.null(mut.maf1[[input.mutations]])){
            stop("Invalid column name for mut.maf")
        }
     
        if (nrow(mut.maf1) == 0){
            stop("Invalid mut.maf, doesn't contain any entries")
        }
    }
    
    if(!missing(mut.maf2)){
        if (is.null(mut.maf2[[input.samples]]) |is.null(mut.maf2[[input.genes]]) | is.null(mut.maf2[[input.mutations]])){
            stop("Invalid column name for mut.maf2")
        }   
        
        if(nrow(mut.maf2) == 0){
            stop("Invalid mut.maf2, doesn't contain any entries")
        }
    }
    
    
    if (gene.cutoff < 0){
        stop("Invalid value for gene cutoff: must be an integer greater than 1 or a real number between 0 and 1")
    }
    print("check1")
    if (gene.cutoff > 1){
        if(!is.wholenumber(gene.cutoff)){
            stop("Invalid value for gene cutoff: must be an integer greater than 1 or a real number between 0 and 1")
        }
    }
    print("check2")
    if (!is.vector(samples)){
        warning("samples not passed in as vector, converting first column")
        samples <- samples[, 1]
    }
    print("check3")
    if (!missing(fixed.order)){
        if (sum(!(fixed.order %in% samples) != 0)){
            stop("Some samples in fixed ordering variable aren't present in samples list")
        }
    }
    print("starting samples")
    all.samples <- unique(samples)
    ## initailize long dataframe with all possible gene-sample combinations
    if(!missing(mut.maf1)){
        print("mut 1 not missing")
        all.samples <- unique(samples)
        all.genes <- unique(mut.maf1[, input.genes])
        df <- expand.grid(all.samples, all.genes, stringsAsFactors = F)
        colnames(df) <- c("samples", "genes")
        df <- df[!is.na(df$samples), ]
        print(head(df))
        
        ## optimized code: use merge + melt to generate df with all called mutations
        print("Initializing mut.maf1")
        print(head(mut.maf1))
        mut.maf1.long <- melt(mut.maf1[, c(input.samples, input.genes, input.mutations)], id = c(input.samples, input.genes))
        print(head(mut.maf1.long))
        mut.maf1.long <- mut.maf1.long[-3]
        colnames(mut.maf1.long) <- c("samples", "genes", "mutations")
        df <- merge(df, mut.maf1.long, c("samples", "genes"), all.x = TRUE)
        print(head(df))
        
        ## annotate coverage information
        if (coverage){
            print("annotating coverage")
            for (i in 1:length(unique(df$samples))){                #print("annotating coverage 2")
                tmp.sample <- unique(df$samples)[i]
                #tmp.sample.ver <- master.sheet[master.sheet$SAMPLE_ACCESSION_NBR == tmp.sample, "PANEL_VERSION"]
                #tmp.sample.ver <- GENIE.samples[GENIE.samples$Sample.Identifier == tmp.sample, "Sequence.Assay.ID"]
                #tmp.sample.ver <- overlap.samples[overlap.samples$Sample.Identifier == tmp.sample, "Sequence.Assay.ID"]
                tmp.sample.ver <- classifier$PANEL_VERSION[classifier$SAMPLE_ACCESSION_NBR == tmp.sample]
                #tmp.sample.ver <- gsc.samples$Sequence.Assay.ID[gsc.samples$Sample.Identifier == tmp.sample]
                ##oncomap coverage stored separately till we figure out how to display
                ## checks to see how many genes for each sample were not covered in that assay, labels them as such
                #if (!tmp.sample.ver){
                    #ghosts <- df[df$samples == tmp.sample & df$genes %in% not.covered.map, ]
                    #if (nrow(ghosts) > 0){
                        #df[df$samples == tmp.sample & df$genes %in% not.covered.map, ]$mutations <- "nc"
                    #}else{
                        ## do nothing
                    #}
                #print("Testing!")
                #print(tmp.sample)
                #if (FALSE) {
                if (classifier$cohort[classifier$SAMPLE_ACCESSION_NBR == tmp.sample] == "DFCI"){
                  #print(i) 
                  #print("annotating coverage 2a")
                  #print("DFCI!")
                  tmp.sample.ver <- as.numeric(tmp.sample.ver)
                  #print(tmp.sample.ver)
                  ghosts <- df[df$samples == tmp.sample & (df$genes %in% not.covered.updated[[tmp.sample.ver]] | !(df$genes %in% gene.list.extended$Gene)), ]
                  #ghosts <- df[df$samples == tmp.sample & (df$genes %in% not.covered[[tmp.sample.ver]] | !(df$genes %in% gene.list$Gene)), ]
                  print("DFCI Ghosts!")
                  print(ghosts)
                  #print(head(df))
                  if (nrow(ghosts) > 0){
                    df$mutations[df$samples == tmp.sample & (df$genes %in% not.covered.updated[[tmp.sample.ver]] | !(df$genes %in% gene.list.extended$Gene))] <- "nc"
                  }else{
                    ## do nothing
                  }
                #}else if (TRUE){
                 }else if (classifier$cohort[classifier$SAMPLE_ACCESSION_NBR == tmp.sample] == "JHH"){
                  #print(i) 
                  #print("annotating coverage 2a")
                  #print("JHH!")
                  #tmp.sample.ver <- as.numeric(tmp.sample.ver)
                  #print(tmp.sample.ver)
                  ghosts <- df[df$samples == tmp.sample & (df$genes %in% not.covered.updated[[tmp.sample.ver]] | !(df$genes %in% gene.list.extended$Gene)), ]
                  #ghosts <- df[df$samples == tmp.sample & (df$genes %in% not.covered[[tmp.sample.ver]] | !(df$genes %in% gene.list$Gene)), ]
                  print("JHH Ghosts!")
                  print(ghosts)
                  #print(head(df))
                  if (nrow(ghosts) > 0){
                    df$mutations[df$samples == tmp.sample & (df$genes %in% not.covered.updated[[tmp.sample.ver]] | !(df$genes %in% gene.list.extended$Gene))] <- "nc"
                  }else{
                    ## do nothing
                  }
                  #}else if (TRUE){
                }else if (classifier$cohort[classifier$SAMPLE_ACCESSION_NBR == tmp.sample] == "GENIE"){
                  #print("GENIE!")
                    #ghosts <- df[df$samples == tmp.sample & !(df$genes %in% GENIE.gene.list.extended.trimmed$Hugo_Symbol[GENIE.gene.list.extended.trimmed$SEQ_ASSAY_ID == tmp.sample.ver]), ]
                  ghosts <- df[df$samples == tmp.sample & !(df$genes %in% GENIE.gene.list.extended.trimmed$Hugo_Symbol), ]
                    #print("Success!")
                    if (nrow(ghosts) > 0){
                        #df[df$samples == tmp.sample & !(df$genes %in% GENIE.gene.list.extended.trimmed$Hugo_Symbol[GENIE.gene.list.extended.trimmed$SEQ_ASSAY_ID == tmp.sample.ver]), ]$mutations <- "nc"
                        df[df$samples == tmp.sample & !(df$genes %in% GENIE.gene.list.extended.trimmed$Hugo_Symbol), ]$mutations <- "nc"
                        print(ghosts)
                    }else{
                        ## do nothing
                    }
                
              }else{
                #print("Other!")
                # ghosts <- df[df$samples == tmp.sample & df$genes, ]
                # #print("Success!")
                # if (nrow(ghosts) > 0){
                #   df[df$samples == tmp.sample & !(df$genes %in% GENIE.gene.list.extended$Hugo_Symbol[GENIE.gene.list.extended$SEQ_ASSAY_ID == tmp.sample.ver]), ]$mutations <- "nc"
                # }else{
                  ## do nothing
              }
            }
          print("TERTp & MGMT Check!")
          print(table(df$mutations[df$genes == "TERT_Promoter"]))
          #print(df[df$samples %in% TCGA.info$SAMPLE_ACCESSION_NBR & df$genes == "TERT_Promoter", ])
          df$mutations[df$samples %in% tcga.tertp.status$Case[is.na(tcga.tertp.status$TERT.promoter.status)] & df$genes == "TERT_Promoter"] <- "nc"
          print(table(df$mutations[df$genes == "TERT_Promoter"]))
          print(table(df$mutations[df$genes == "MGMT"]))
          #print(df[df$samples %in% TCGA.info$SAMPLE_ACCESSION_NBR & df$genes == "TERT_Promoter", ])
          #print(sum(df$samples %in% master.sheet.glioma$SAMPLE_ACCESSION_NBR[master.sheet.glioma$MGMT_Mehdi == "N/A"]))
          #print(sum(df$samples %in% master.sheet.glioma$SAMPLE_ACCESSION_NBR[master.sheet.glioma$MGMT_Mehdi == "N/A"] & df$genes == "MGMT"))
          #df$mutations[df$samples %in% master.sheet.glioma$SAMPLE_ACCESSION_NBR[master.sheet.glioma$MGMT_Mehdi == "N/A"] & df$genes == "MGMT"] <- "nc"
          print(table(df$mutations[df$genes == "MGMT"]))
        }
        # print("STS Check!")
        # print(table(df$mutations[df$samples %in% sarc.samples$Sample.Identifier]))
        ## orders rows based on gene most frequently altered genes, removing NAs and uncovered genes
        #df$genes[df$genes %in% c("CDKN2A", "CDKN2B")] <- "CDKN2A/B"
        df_sub <- subset(df, !is.na(df$mutations))
        df_sub <- subset(df, df$mutations != "nc")
        ord <- names(sort(table(df_sub$genes), decreasing = T))
        df.combined <- df
    }
    ## checks if second maf of mutations supplied. If so, duplicate above code to produce second df
    if(!missing(mut.maf2)){
        ## iniatialize df.2
        all.samples <- unique(samples)
        all.genes <- unique(mut.maf2[, input.genes])
        df.2 <- expand.grid(all.samples, all.genes, stringsAsFactors = F)
        colnames(df.2) <- c("samples", "genes")
        
        print("Initializing mut.maf2")
        mut.maf2.long <- melt(mut.maf2[, c(input.samples, input.genes, input.mutations)], id = c(input.samples, input.genes))
        mut.maf2.long <- mut.maf2.long[-3]
        colnames(mut.maf2.long) <- c("samples", "genes", "mutations")
        df.2 <- merge(df.2, mut.maf2.long, c("samples", "genes"), all.x = TRUE)
        
        if (coverage){
            ## annotate coverage information
            print("annotating coverage")
            for (i in 1:length(unique(df.2$samples))){
                tmp.sample <- unique(df.2$samples)[i]
                tmp.sample.ver <- classifier[classifier$SAMPLE_ACCESSION_NBR == tmp.sample, "PANEL_VERSION"]
                
                ##oncomap coverage stored separately till we figure out how to display
                ## checks to see how many genes for each sample were not covered in that assay, labels them as such
                if (!tmp.sample.ver){
                    ghosts <- df.2[df.2$samples == tmp.sample & df.2$genes %in% not.covered.map, ]
                    if (nrow(ghosts) > 0){
                        df.2[df.2$samples == tmp.sample & df.2$genes %in% not.covered.map, ]$mutations <- "nc"
                    }else{
                        ## do nothing
                    }
                }else{
                    ghosts <- df.2[df.2$samples == tmp.sample & df.2$genes %in% not.covered.updated[[tmp.sample.ver]], ]
                    if (nrow(ghosts) > 0){
                        df.2[df.2$samples == tmp.sample & df.2$genes %in% not.covered.updated[[tmp.sample.ver]], ]$mutations <- "nc"
                    }else{
                        ## do nothing
                    }
                }
            }
        }
        
        df_sub <- subset(df.2, !is.na(df.2$mutations))
        df_sub <- subset(df.2, df.2$mutations != "nc")
        
        ## generate empty row to separate tier 4 from rest
        df.3 <- expand.grid(all.samples, "Tier4")
        df.3$mutations <- "nc"
        colnames(df.3) <- c("samples", "genes", "mutations")
        
        ord.2 <- names(sort(table(df_sub$genes), decreasing = T))
        
        ## checks if first supplied argument or not, either updates or creates tracking variables as appropriate
        if(exists("df.combined")){
            df.combined <- rbind(df.combined, df.3, df.2)
            ord <- c(ord, "Tier4", ord.2)
        }else{
            df.combined <- df.2
            ord <-ord.2
        }
        
    }
    
    ## keep given % of genes based on cutoff if mutations data supplied
    if(!missing(mut.maf1) | !missing(mut.maf2)){
        if (gene.cutoff == 1){
            ## keep all genes
        }else if (gene.cutoff < 1){
            ## keep given percentage of top hits
            idx <- ceiling(length(ord) * gene.cutoff)
            ord <- ord[1:idx]
        }else{
            ## keep given number of top hits
            idx <- min(gene.cutoff, length(ord))
            ord <- ord[1:idx]
        }
    }
    
    ## checks to see if focal.cnv information supplied to comut
    print ("starting focal")
    if (!missing(focal.cnv.maf)){
        if (!is.na(focal.cnv.maf)){
            print("focal")
            all.samples <- unique(samples)
            ## generates df with all possible cnvs
            all.genes <- unique(focal.cnv.maf$GENE)
            df.cnv <- expand.grid(all.samples, all.genes, stringsAsFactors = F)
            print(head(df.cnv))
            colnames(df.cnv) <- c("samples", "genes")
            
            print("Initializing df.cnv.long")
            colnames(focal.cnv.maf) <- c("samples", "genes", "mutations")
            df.cnv <- merge(df.cnv, focal.cnv.maf, c("samples", "genes"), all.x = TRUE)
            
            # annotes gene coverage
            print("annotating coverage")
            for (i in 1:length(unique(df.cnv$samples))){
                tmp.sample <- unique(df.cnv$samples)[i]
                tmp.sample.ver <- classifier[classifier$SAMPLE_ACCESSION_NBR == tmp.sample, "PANEL_VERSION"]

                ## checks to see how many genes for each sample were not covered in that assay, labels them as such
                ghosts <- df.cnv[df.cnv$samples == tmp.sample & df.cnv$genes %in% not.covered.cnv[[tmp.sample.ver]], ]
                if (nrow(ghosts) > 0){
                    df.cnv[df.cnv$samples == tmp.sample & df.cnv$genes %in% not.covered.cnv[[tmp.sample.ver]], ]$mutations <- "nc"
                }else{
                    ## do nothing
                }

            }
            
            df.cnv$genes <- sapply(df.cnv$genes, paste, "-cnv", sep = "", USE.NAMES = F)
            
            ## generate subset with mutations for accurate counting
            df_sub <- df.cnv[!is.na(df.cnv$mutations), ]
            ord.cnv <- names(sort(table(df_sub$genes), decreasing = T))
            
            if(exists("df.combined")){
                df.combined <- rbind(df.combined, df.cnv)
                ord <- c(ord, ord.cnv)
            }else{
                df.combined <- df.cnv
                ord <- ord.cnv
            }
            
        
        }
    }
    
    ## checks to see if broad.cnv information supplied to comut
    if (!missing(broad.cnv.maf)){
        if (!is.na(broad.cnv.maf)){
            print("broad not missing")
            ## renames df to enable merging with mutations df, adds to previous dataframe, 
            ## then inserts names at end of dataframe for factor ordering
            
            ## first transpose so samples are columns
            missing_samples <- all.samples[!(all.samples %in% colnames(broad.cnv.maf))]
            #missing_samples <- unique(samples)
            cnv_missing <- matrix(nrow = nrow(broad.cnv.maf), ncol = length(missing_samples))
            cnv_missing <- as.data.frame(cnv_missing)
            rownames(cnv_missing) <- rownames(broad.cnv.maf)
            colnames(cnv_missing) <- missing_samples
            #print(dim(cnv_missing))
            cnv_missing[ , ] <- 3
            print("number of samples w/out broad calls")
            print(dim(cnv_missing))
            print(head(colnames(cnv_missing)))
            broad.cnv.maf <- cbind(broad.cnv.maf, cnv_missing)
            print(head(colnames(broad.cnv.maf)))            
            broad.cnv.maf <- t(broad.cnv.maf)
            broad.cnv.maf <- as.data.frame(broad.cnv.maf)
            broad.cnv.maf$samples <- rownames(broad.cnv.maf)
            print(head(broad.cnv.maf$samples))
            broad.cnv.maf$samples <- gsub("\\.", "\\-", broad.cnv.maf$samples)
            print(head(broad.cnv.maf$samples))
            broad.cnv.maf <- melt(broad.cnv.maf, id = "samples")
            
            colnames(broad.cnv.maf) <- c("samples", "genes", "mutations")
            #print(head(broad.cnv.maf))
            ## convert 1 and -1 to gain and loss for clarity
            print("converting numbers to labels")
            broad.cnv.maf$mutations[broad.cnv.maf$mutations %in% c(-1, -1.5)] <- "arm-level loss"
            broad.cnv.maf$mutations[broad.cnv.maf$mutations %in% c(1, 2)] <- "arm-level gain"
            #broad.cnv.maf$mutations[broad.cnv.maf$mutations %in% c(1, 2)] <- "altered"
            #broad.cnv.maf$mutations[is.na(broad.cnv.maf$mutations)] <- "?"
            broad.cnv.maf$mutations[is.na(broad.cnv.maf$mutations)] <- "not covered"
            broad.cnv.maf$mutations[broad.cnv.maf$mutations == 0.46] <- "Mixed"
            broad.cnv.maf$mutations[broad.cnv.maf$mutations == 0] <- NA
            broad.cnv.maf$mutations[broad.cnv.maf$mutations == 3] <- "ncc"
            print("assigning order")
            print(head(broad.cnv.maf))
            ord.broad.cnv <- names(sort(table(broad.cnv.maf$genes), decreasing = T))
            print(ord.broad.cnv)
            print(table(broad.cnv.maf$mutations))
            
            ## add to df
            print("adding to df")
            if(exists("df.combined")){
                df.combined <- rbind(df.combined, broad.cnv.maf)
                print(dim(df.combined))
                ord <- c(ord, ord.broad.cnv)
            }else{
                df.combined <- broad.cnv.maf
                ord <- ord.broad.cnv
            }
        }
    }
    
    
    ## checks to see if additional phenotype information supplied for comut
    print("phenotypes")
    if (!missing(phenotypes)){
        ## transforms to long format so it can be added to df
        df.pheno <- melt(phenotypes, id = colnames(phenotypes)[1])
        
        ## renames df to enable merging with mutations df, adds to previous dataframe, 
        ## then inserts names at end of dataframe for factor ordering
        colnames(df.pheno) <- c("samples", "genes", "mutations")
        df.combined <- rbind(df.combined, df.pheno)
        ord <- c(ord, colnames(phenotypes)[-1])
        print(table(df.combined$mutations[df.combined$genes == "Subtype"]))
        print(head(phenotypes))
    }
    ## if manual ordering supplied, resorts for factor ordering
    print("man order check")
    if (!missing(manual.order)){
        ## make sure manual names are actually in df
        if (sum(!(ord %in% df.combined$genes)) > 0){
            stop(paste("Manual order term contains genes not found in data frame: ", ord[!(ord %in% df.combined$genes)]))
        }else {
            ord <- c(manual.order, ord)
            ord <- ord[!duplicated(ord)]
        }   
    }
    
    ## creats factors based on levels, removed genes not present in ord conditions
    print("df wide")
    df.combined$genes <- factor(df.combined$genes, levels = ord)
    df.combined <- df.combined[order(df.combined$genes), ]
    df.combined <- df.combined[!is.na(df.combined$genes), ]
    print(unique(df.combined$genes))
    
    ## reshape long df to wide format to determine ordering of samples (columns)
    df.wide <- reshape(df.combined, v.names = "mutations", idvar = "samples", timevar = "genes", direction = "wide")
    print(colnames(df.wide))
    for (i in 2:ncol(df.wide)){
        ## change NA's and uncovered to filler which will always order last
        df.wide[, i][is.na(df.wide[, i])] <- "z"
        df.wide[, i][df.wide[, i] == "nc"] <- "z"
        df.wide[, i][df.wide[, i] == "ncc"] <- "z"
        df.wide[, i][df.wide[, i] == "not covered"] <- "z"
        ## change mutations to "mutation" so type doesn't factor into cascading order
        df.wide[, i][df.wide[, i] %in% c("frameshift_indel", "missense", "nonsense", "splice_site", "stop_codon", "in_frame_indel", "other", "TSS", 
                                         "damaging mutation", "focal gain", "focal loss", "rearrangement")] <- "mutation"
    }
    
    ##sort age manually
    print("fillered")
    
    ## sorts by single gene if only one significantly mutated, otherwise all columns
    if (ncol(df.wide) == 2){
        df.wide <- df.wide[order(df.wide[, 2]), ]
        
    }else{
      if ("mutations.Age" %in% colnames(df.wide)){
         df.wide$mutations.Age <- factor(df.wide$mutations.Age, levels = c("Pediatric", "Young Adult", "Middle", "Old", "Unknown"))
         #df.wide$mutations.Subtype <- factor(df.wide$mutations.Subtype, levels = c("Oligo", "Other LGG", "Astro", "Glioblastoma"))
      }
        print("do.call")
        df.wide <- df.wide[do.call("order", df.wide[, -1]), ]
    }
    print(df.wide[1:5, 1:5])
    
    ## takes given order, and sorts samples based on ordering
    ## add on samples with no mutations to factor creation so they don't get NA'd
    print("missing again")
     missing <- df.combined$sample[!(df.combined$samples %in% df.wide$samples)]
    df.combined$samples <- factor(df.combined$samples, levels = c(df.wide$sample, missing))
    
    ## overwrites ordering if fixed order supplied
    if(!missing(fixed.order)){
        df.combined$samples <- factor(df.combined$samples, levels = fixed.order)
    }
    
    ## final cleanup
    print("final cleanup")
    
    ## switches factor order back for plotting so most frequent is highest on y axis
    print(levels(df.combined$genes))
    if (pheno_top){
      print(c(colnames(phenotypes)[colnames(phenotypes) %in% df.combined$genes], levels(df.combined$genes)[!(levels(df.combined$genes) %in% colnames(phenotypes))]))
      #print(rev(levels(df.combined$genes)[!(levels(df.combined$genes) %in% colnames(phenotypes))]))
      df.combined$genes <- factor(df.combined$genes, c(colnames(phenotypes)[colnames(phenotypes) %in% df.combined$genes], levels(df.combined$genes)[!(levels(df.combined$genes) %in% colnames(phenotypes))]))
    }
    df.combined$genes <- factor(df.combined$genes, rev(levels(df.combined$genes)))
    print(levels(df.combined$genes))
    #print(table(df.combined$mutations[df.combined$genes %in% c("1p", "7p", "10q", "19q")]))
    df.combined$mutations[is.na(df.combined$mutations)] <- "wt"
    #print(table(df.combined$mutations[df.combined$genes %in% c("1p", "7p", "10q", "19q")]))
    #print(dim(df.combined))
    #df.combined$mutations[is.na(df.combined$mutations)] <- "Other LGG"
    df.combined$mutations[df.combined$mutations %in% c("nc", "ncc")] <- "not covered"
    #print(dim(df.combined))
    print(head(df.combined[!(df.combined$samples %in% samples), ]), 20)
    df.combined <- df.combined[df.combined$samples %in% samples,]
    #print(dim(df.combined))
    print(table(df.combined$mutations))
    print("plotting")
    
    ## generate color scheme for plotting
    
    if (!missing(col.vector)){
        print("using supplied color scheme")
        
    }else{
        print("using random color scheme")
        
        ## first set defaults
        default.names <- c("arm-level gain", "arm-level loss", "HA", "2DEL", "focal loss", "focal gain")
        default.colors <- c("red", "blue", "red", "blue", "red", "blue")
        names(default.colors) <- default.names
        col.vector <- default.colors
    }

    ## gets arguments that don't have assigned color
    missing.names <- unique(df.combined$mutations)[!(unique(df.combined$mutations) %in% names(col.vector))]
    print(col.vector)
    print(length(missing.names))
    if (length(missing.names > 0)){
      missing.colors <- distinctColorPalette(length(missing.names))
      print(missing.colors)
      names(missing.colors) <- missing.names
      print(names(missing.colors))
      col.vector <- c(col.vector, missing.colors)
    }
    print(col.vector)
    
    ## sets wt to grey, then removes those factors that aren't present
    print("assigning colors")
    col.vector[which(names(col.vector) == "wt")] <- "beige"
    col.vector[which(names(col.vector) %in% c("nc", "ncc", "not covered"))] <- "white"
    col.vector[which(names(col.vector) == "?")] <- "grey"
    col.vector <- col.vector[names(col.vector) %in% unique(df.combined$mutations)]
    print(col.vector)
    colScale <- scale_fill_manual(values = col.vector)
    
    ## generate plot
    print("generating plot")
    mut <- ggplot(df.combined, aes(x=samples, y=genes, height=0.8, width=1)) + 
        geom_tile(aes(fill=mutations)) +
        colScale +
        ggtitle(title) +
        theme(
            legend.key = element_rect(fill='NA'),
            legend.key.size = unit(0.4, 'cm'),
            legend.title = element_blank(),
            legend.position="bottom",
            legend.text = element_text(size=legend.font.size, face="bold"),
            axis.ticks.x=element_blank(),
            axis.ticks.y=element_blank(),
            #axis.text.x=element_text(angle = 90, hjust = 1),
            axis.text.x=element_blank(),
            axis.text.y=element_text(colour="Black", size = y.axis.font.size),
            axis.title.x=element_text(face="bold"),
            axis.title.y=element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.background=element_blank()
        )
    print(mut) 
    
    ## prints to screen or to filepath
    if(missing(file.path)){
       print(mut) 
    }else{
        print(paste("printing Comut to ", file.path))
        pdf(file.path, width = dimensions[1], height = dimensions[2])
        print(mut)
        dev.off()
    }
    
    if(return.matrix){
        return(df.combined)
    }
    
    if(return.plot){
        #return(mut)
        return(as.character(df.combined$samples))
    }
}


PlotHistogramCoverage <- function(maf, gene.column, samples, title, type = "mut"){
## Takes in a maf file, and generates a histogram plot by percent samples affected, controling for uncovered genes

## Args:
    ## maf: a maf file
    ## gene.column: the column in maf with gene mutation information
    ## samples: a vector with all samples
    
    ## error checking

    if (is.null(maf[[gene.column]])){
        stop("Invalid column name")
    }
    
    if(missing(samples)){
        stop("Missing sample vector")
    }
    
    if (type == "focal"){
      print("doing focal analysis")
      print(paste("No. samples input: ", length(samples)))
      samples <- samples[samples %in% c(colnames(GENIE.cna), TCGA.GBM.focal.wide$Tumor_Sample_Barcode, TCGA.LGG.focal.wide$Tumor_Sample_Barcode, classifier$SAMPLE_ACCESSION_NBR)]
      print(paste("samples with cnv coverage: ", length(samples)))
    }
    ## iniatlize dataframe with naive incidence counts, correcting for multiple samples
    print("initializing")
    maf <- maf[maf$SAMPLE_ACCESSION_NBR %in% samples, ]
    genes <- maf[, gene.column]
    print(sum(genes == "TERT_Promoter"))
    tbl <- table(genes)
    df <- data.frame(tbl)
    
    ## Returns the total number of samples that the gene was covered in
    TrueDenominator <- function(gene){
        ## figure out which versions it was in
        if (type == "focal") {
          print("Focal!")
          dfci.idx <- c("1", "2", "3", "3.1")[gene.list.extended[gene, c(3,5,8,8)] != 0] 
          genie.idx <- unique(GENIE.samples$Sequence.Assay.ID[GENIE.samples$Sample.Identifier %in% colnames(GENIE.cna[ ,-1])[!is.na(GENIE.cna[gene,-1])]])
          tcga.idx <- c("TCGA_LGG", "TCGA_GBM")[c(gene %in% colnames(TCGA.LGG.focal.wide), gene %in% colnames(TCGA.GBM.focal.wide))]
          jhh.idx <- c("STPv2", "STPv3", "STPv3.2", "STPv4", "STPL", "NPP")[gene.list.extended[gene, c(12,13,14,15,16,17)] != 0] 
          idx <- c(dfci.idx, genie.idx, tcga.idx, jhh.idx)
          idx <- idx[!is.na(idx)]
        }
      if (type == "sv") {
        print("SV!")
        dfci.idx <- c("1", "2", "3", "3.1")
        genie.idx <- unique(GENIE.samples$Sequence.Assay.ID[GENIE.samples$Sequence.Assay.ID %in% GENIE.data$Sequence.Assay.ID[GENIE.data$Sample.Identifier %in% GENIE.fusions.all$Tumor_Sample_Barcode]])
        tcga.idx <- c("TCGA_LGG", "TCGA_GBM")
        jhh.idx <- c("STPv2", "STPv3", "STPv3.2", "STPv4", "STPL", "NPP")
        idx <- c(dfci.idx, genie.idx, tcga.idx, jhh.idx)
        idx <- idx[!is.na(idx)]
      }
        else {
        print("Other variant section!")
        print(head(idx))
        #idx <- c(c("1", "2", "3", "3.1")[(as.numeric(gene.list.extended[gene, c(3,5,8,8)])+as.numeric(gene.list.extended[gene, c(4,6,9,9)])+c(0,as.numeric(gene.list.extended[gene,7]),0,0)) != 0], unique(GENIE.gene.list.extended.trimmed$SEQ_ASSAY_ID[GENIE.gene.list.extended.trimmed$Hugo_Symbol == gene]), "TCGA_LGG", "TCGA_GBM",c("STPv2", "STPv3", "STPv3.2", "STPv4", "STPL", "NPP")[gene.list.extended[gene, c(12,13,14,15,16,17)] != 0])
        idx <- c(c("1", "2", "3", "3.1", "STPv2", "STPv3", "STPv3.2", "STPv4", "STPL", "NPP")[(as.numeric(gene.list.extended[gene, c(3,5,8,8,12,13,14,15,16,17)])+as.numeric(gene.list.extended[gene, c(4,6,9,9)])+c(0,as.numeric(gene.list.extended[gene,7]),0,0)) != 0], unique(GENIE.gene.list.extended.trimmed$SEQ_ASSAY_ID[GENIE.gene.list.extended.trimmed$Hugo_Symbol == gene]), "TCGA_LGG", "TCGA_GBM")
        idx <- idx[!is.na(idx)]
        
        }
        if (gene == "CDKN2A"){
          print("checking focal coverage")
          print(gene)
          print(idx)
          #print(dfci.idx)
        }
        #idx <- unique(mutations_combined$PANEL_VERSION[mutations_combined$BEST_EFF_GENE == gene])
        #rev.idx <- !idx
        #versions <- names(table(master.sheet$PANEL_VERSION))
        
                ## gets the number of samples that were sequenced with the versions containing that gene
        #filtered <- subset(master.sheet, SAMPLE_ACCESSION_NBR %in% samples & PANEL_VERSION %in% versions[rev.idx])
        #filtered <- master.sheet_all
        #print(dim(filtered))
        #print(sum(samples %in% filtered$SAMPLE_ACCESSION_NBR))
        filtered <- subset(classifier, SAMPLE_ACCESSION_NBR %in% samples & PANEL_VERSION %in% idx)
        #filtered <- subset(overlap.samples, Sample.Identifier %in% samples & Sequence.Assay.ID %in% idx)
        #print(gene)
        if (gene == "TERT_Promoter"){
          print("checking EWSR1")
          print(idx)
          print(dim(subset(classifier, SAMPLE_ACCESSION_NBR %in% samples & PANEL_VERSION %in% idx)))
          #print(versions)
          print(dim(filtered))
          print(table(filtered$cohort))
          filtered <- subset(filtered, SAMPLE_ACCESSION_NBR %in% tcga.tertp.status$Case[!is.na(tcga.tertp.status$TERT.promoter.status)] | filtered$cohort != "TCGA")
          print(dim(filtered))
          #print(head(samples))
          #print(versions[rev.idx])
          print(length(samples))
          print("BL-13-A24307" %in% samples)
          print(head(filtered))
          print(table(filtered$cohort))
          #print(genes)
          print(nrow(filtered))
        } 
        if (gene == "CDKN2A"){
          print("checking CDKN2A")
          print(idx)
          #print(versions)
          print(dim(filtered))
          #filtered <- subset(filtered, SAMPLE_ACCESSION_NBR %in% tcga.tertp.status$Case[!is.na(tcga.tertp.status$TERT.promoter.status)] | filtered$cohort != "TCGA")
          #print(dim(filtered))
          #print(head(samples))
          #print(versions[rev.idx])
          print(length(samples))
          print("BL-13-A24307" %in% samples)
          print(head(filtered))
          print(table(filtered$cohort))
          #print(genes)
          print(nrow(filtered))
        } 
        
        nrow(filtered)
        
        
    }
    
    
    ## creates dataframe with denominator counts for each gene
    print("finding denominators")
    #print(unique(genes))
    denominators <- sapply(unique(genes), TrueDenominator)
    denominators.df <- data.frame(denominators, names(denominators))
    colnames(denominators.df)[2] <- "genes"
    print(head(denominators.df))
    #denominators.df <- denominators.df[denominators.df$denominators > length(samples)/10, ]
    print(length(samples)/10)
    
    ## merges two dataframes, then create column for plotting
    df.merged <- merge(df, denominators.df, by = "genes")
    df.merged$percent <- df.merged$Freq / df.merged$denominators
    df.merged <- df.merged[order(df.merged$percent, decreasing = T), ]
    df.merged$genes <- factor(df.merged$genes, levels = df.merged$genes)
    print(head(df.merged))
    print("plotting")
    if (length(unique(df.merged$genes))>500) {
      print("truncating data for plotting")
      plot <- ggplot(data = df.merged[df.merged$percent > 0.05, ], aes(x=genes, y = percent)) + geom_bar(stat = "identity") + ggtitle(label =  title) + rameen_theme
    }
    else{
    print(title)
    plot <- ggplot(data = df.merged, aes(x=genes, y = percent)) + geom_bar(stat = "identity") + ggtitle(label =  title) + rameen_theme
    }
    print(plot)
    return(df.merged)

}

EntropyHelper <- function(counts){
    ## takes in the number of samples in each split, returns the entropy
    
    ## Args:
        ## counts: vector of at least length 2 with the number of samples in each division
    
    ## error checking
    
    if (length(counts) < 2){
        #stop("Insufficient number of counts")
      print("Insufficient number of counts")
      0
    }
    
    if (any(is.na(counts))){
        stop("NAs detected")
    }
    
    if (any(counts < 0)){
        stop("Must be non-negative")
    }
    
    ## turn counts into proportions, exclude 0s, then calculate entropy
    proportions <- counts / sum(counts)
    
    proportions <- proportions[proportions != 0]
    
    sum(sapply(proportions, function(x){-x * log2(x)}))
    
}




CalculateEntropy <- function(tree, detailed = NA){
    ## takes in a matrix of samples and features, then determines which column segregregates the data to minimize the entropy
    
    ## Args:
        ## tree: a matrix or data frame. Columns are features of interest, rows are samples
        ## detailed: optional argument that will produce detailed breakdown of the entropy contribution of each gene to a 
        ## column of interest
        
        
  
    ## error checking
    if (ncol(tree) < 3){
        stop("Insufficient number of columns")
    }
    
    if (nrow(tree) < 3){
        print(dim(tree))
        #print(tree)
        stop("Insufficient number of rows")
    }
    
      tree <- apply(tree, 2, as.numeric)
      ## determines if any columns need to be skipped
      #print("determining skips")
      #sums <- colSums(tree)
      sums <- colSums(tree, na.rm = TRUE)
      skip.idx <- sums == 0| sums == colSums(!is.na(tree))
      #skip.idx <- sums < 0.05*colSums(!is.na(tree))| sums > 0.95*colSums(!is.na(tree))
      #print(sum(!skip.idx))
      tree <- tree[, !skip.idx]
      #print(skip.idx)
      #print(dim(tree))
    
      if (sum(!skip.idx)<2){
      cumulative.information.gain <- 0
      print("setting cum info gain to 0 -- v1")
      }
      
    else{
    ## determines which counts of "i" will be marked
    #print("detailed?")
    #print(missing(detailed))
    if (!missing(detailed)){
        #detailed.idx <- which(colnames(tree) %in% detailed)
        detailed.idx <- 1:length(colnames(tree))
        detailed.counts <- list()
    }else{
        detailed.idx <- c()
    }
    
    ## calculate entropy
    ## stores entropy value for each feature of interest
    cumulative.information.gain <- c()
    sums2 <- colSums(tree, na.rm = TRUE)
    skip.idx.2 <- sums2 < 0.05*colSums(!is.na(tree))| sums2 > 0.95*colSums(!is.na(tree))
    #print(skip.idx.2)
    print(colnames(tree)[skip.idx.2])
    #print(dim(tree))
    #print(ncol(tree))
    #for (i in (1:ncol(tree))[!skip.idx.2]){
    for (i in (1:ncol(tree))){
        #print(colnames(tree)[i])
        
        ## stores the cumulative information gain from all pairwise comparisons
        total.information.gain <- 0
        #branch.entropy <- 0
        
        ## stores individual pair-level information for features for plotting
        detailed.information.gain <- c()
        
        for(j in 1:ncol(tree)){
            if (i == j){
                ## skip
            }else{
                ## split data based on identified feature, figure out total on each side
                x <- table(as.numeric(tree[, i]), as.numeric(tree[, j]))
                ## gets number of samples assigned to each branch
                branch.totals <- rowSums(x)
                #branch.percents <- branch.totals / nrow(tree)
                branch.percents <- branch.totals / sum(branch.totals)
                ## calculates entropy within each split
                if (i == 1 & j > 2) {
                  #print(entropies)
                  #print(j)
                  #print(colnames(tree)[c(i,j)])
                  #print(x)
                  #print(branch.percents)
                }
                entropies <- vector("numeric", nrow(x))
                
                for (k in 1:nrow(x)){
                    entropies[k] <- EntropyHelper(x[k, ])
                }
                
                
                
                ## multiply entropy by percent to get weighted value
                entropy <- sum(branch.percents * entropies)
                
                ## Information gain is calculated as entropy of initial split minus entropy of each feature
                #prior.counts <- table(tree[, j])
                prior.counts <- colSums(x)
                prior.entropy <- EntropyHelper(prior.counts)
                information.gain <- prior.entropy - entropy
                
                #branch.entropy <- branch.entropy + prior.entropy
                total.information.gain <- total.information.gain + information.gain
                if (i == 1 & j > 2) {
                  #print(entropy)
                  #print(information.gain)
                  #print(total.information.gain)
                }
                ## pulls out entropy contribution from individual features of interest
                if (i %in% detailed.idx){
                    detailed.information.gain <- c(detailed.information.gain, information.gain)
                    names(detailed.information.gain)[length(detailed.information.gain)] <- colnames(tree)[j]
                }
            }
        }
        
        
            
        cumulative.information.gain[i] <- total.information.gain
        #cumulative.entropy[i] <- branch.entropy
        names(cumulative.information.gain)[i] <- colnames(tree)[i]
        
        
        if (i %in% detailed.idx){
            detailed.information.gain <- sort(detailed.information.gain, decreasing = T)
            detailed.counts[[length(detailed.counts) + 1]] <- detailed.information.gain
            names(detailed.counts)[length(detailed.counts)] <- colnames(tree)[i]

        }
    
    }
    print(which(is.na(cumulative.information.gain)))
    #cumulative.information.gain <- cumulative.information.gain[!skip.idx.2]
    }
      
    #cum.ig.entropy <- c(cumulative.information.gain, cumulative.entropy)
    if (!missing(detailed)){
        print(which(is.na(detailed.counts)))
        new.detailed.idx <- which(names(detailed.counts) %in% names(sort(cumulative.information.gain, decreasing = TRUE))[1:4])
        return(list(cumulative.information.gain, detailed.counts[new.detailed.idx], tree))
    }else{
        #print(max(cumulative.information.gain))
        return(cumulative.information.gain) 
        #return(cum.ig.entropy)
    }
    
}


BranchPlotter <- function(entropy.data, hit.number = 4){
    ## Takes an entropy object and plots the data at each branch point
    
    ## Args:
    ## entropy data: object created by CalculateEntropy
    ## hit.number: the number of features to show detailed information on
    ## incidence: boolean that determines whether incidence plot is generated
    if (!is.list(entropy.data)){
        summary.plot.data <- data.frame(entropy.data, names(entropy.data), stringsAsFactors = FALSE)
    }else{
        summary.plot.data <- data.frame(entropy.data[[1]], names(entropy.data[[1]]), stringsAsFactors = FALSE)
    }
    
    ## plot total information gain across all features
    colnames(summary.plot.data) <- c("information_gain", "feature")
    summary.plot.data$feature <- factor(summary.plot.data$feature, levels = summary.plot.data$feature[order(summary.plot.data$information_gain, decreasing = TRUE)])
    summary.plot.data <- summary.plot.data[summary.plot.data$information_gain > 0 & summary.plot.data$feature %in% levels(summary.plot.data$feature)[1:20], ]
    summary.plot <- ggplot(data = summary.plot.data, aes(x = feature, y = information_gain)) + geom_bar(stat = "identity") + 
        labs(title = "Features ranked by net information gain") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    detailed.plots <- names(entropy.data[[2]])
    print(detailed.plots)
    return.plots <- list(summary.plot)
    
    ## plots specific information gain contribution of each gene for those features specificied in initial entropy function call
    if (hit.number != 0){
        for (i in 1:length(detailed.plots)){
            subplot.idx <- which(names(entropy.data[[2]]) == detailed.plots[i])
            subplot.data <- data.frame(entropy.data[[2]][[subplot.idx]], names(entropy.data[[2]][[subplot.idx]]), stringsAsFactors = FALSE)
            colnames(subplot.data) <- c("information_gain", "feature")
            subplot.data$feature <- factor(subplot.data$feature, levels = subplot.data$feature[order(subplot.data$information_gain, decreasing = TRUE)])
            subplot.data <- subplot.data[subplot.data$information_gain > 0.01, ]
            plot.title <- paste("Information gain by feature for", detailed.plots[i])
            print(plot.title)
            return.plots[[i + 1]] <- ggplot(data = subplot.data, aes(x = feature, y = information_gain)) + geom_bar(stat = "identity") + 
                labs(title = plot.title) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(limits = c(0, 1))
        }
    }
    
    ## generates plot showing incidence of top features    
    incidence.plot.data <- colSums(entropy.data[[3]], na.rm = T)
    incidence.plot.data <- data.frame(names(incidence.plot.data), as.numeric(incidence.plot.data), stringsAsFactors = FALSE)
    colnames(incidence.plot.data) <- c("feature", "count")
    incidence.plot.data <- incidence.plot.data[incidence.plot.data$feature %in% detailed.plots, ]
    incidence.plot.data$feature <- factor(incidence.plot.data$feature, levels = levels(summary.plot.data$feature))
    incidence.plot <- ggplot(data = incidence.plot.data, aes (x = feature, y = count)) + geom_bar(stat = "identity") +
        labs(title = "Number of samples with alteration") + scale_y_continuous(limits = c(0, nrow(entropy.data[[3]])))
    return.plots[[length(return.plots) + 1]] <- incidence.plot
    
    return(return.plots)
}


OptimizedHelper <- function(tree, idx.1, idx.2){
    ## functionalizes inner for loop for feed in to apply
    x <- table(as.numeric(tree[, idx.1]), as.numeric(tree[, idx.2]))
    ## gets proportion of total samples in each branch
    branch.totals <- rowSums(x)
    #branch.percents <- branch.totals / nrow(tree)
    branch.percents <- branch.totals / sum(branch.totals)
    
    ## calculates entropy within each split
    entropies <- vector("numeric", nrow(x))
    for (i in 1:nrow(x)){
        entropies[i] <- EntropyHelper(x[i, ])
    }
    ## multiply entropy by percent to get weighted value
    entropy <- sum(branch.percents * entropies)
    ## Information gain is calculated as entropy of initial split minus entropy of each feature
    #prior.counts <- table(tree[, idx.2])
    prior.counts <- colSums(x)
    prior.entropy <- EntropyHelper(prior.counts)
    information.gain <- prior.entropy - entropy
    if (is.na(information.gain)){
      print(colnames(tree)[idx.1])
      print(colnames(tree)[idx.2])
      print(dim(tree))
      print(x)
      print(entropies)
      print(entropy)
      print(prior.counts)
      print(prior.entropy)
      print(information.gain)
      stop()
    }
    return(information.gain)
}

OptimizedCalculateEntropy <- function(tree, correction.method = "none", detailed = NA){
    ## trimmed down version for permutations testing
    ## requires tree to be input as numeric object, with full/empty columns removed
   num.cols <- ncol(tree)
 
    ## calculate entropy
    cumulative.entropy <- rep(0, num.cols)
    cumulative.entropy <- sapply(1:num.cols, function(i){sum(sapply((1:num.cols)[-i], function(x){OptimizedHelper(tree, i, x)}))})
    # for (i in 1:num.cols){
    #     inner.idx <- 1:num.cols
    #     inner.idx <- inner.idx[inner.idx != i]
    #     #print(paste("i:" , colnames(tree)[i]))
    #     cumulative.entropy[i] <- sum(sapply(inner.idx, function(x){OptimizedHelper(tree, i, x)}))
    #     if (is.na(cumulative.entropy[i])){
    #       print(paste("i:" , colnames(tree)[i]))
    #       notreal <- (colnames(tree)[is.na(sapply(inner.idx, function(x){OptimizedHelper(tree, i, x)}))])
    #       print(notreal)
    #       print(table(as.numeric(tree[, i]), as.numeric(tree[, notreal[1]])))
    #       stop()
    #     }
    # }
    
    return(cumulative.entropy)    
}

PermuteDecisionTree <- function(tree, permutation.number){
    ## Takes a decision tree, permutes the data, then recalculates entropy to determine whether cutoffs are significant
    
    ## Args:
    ## tree: a decision tree
    ## permutation number: number of iterations to permute through
     sums <- colSums(tree, na.rm = TRUE)
     skip.idx <- sums == 0| sums == colSums(!is.na(tree))
     tree <- tree[, !skip.idx]
    #tree <- tree[, colSums(tree) != 0 & colSums(tree) != nrow(tree)]
    sample.number <- nrow(tree)
    permuted.values <- rep(0, permutation.number)
    permuted.names <- rep("", permutation.number)
    #permuted.trees <- permatswap(tree, mtype = "prab", times = permutation.number)
    permuted.trees <- simulate(nullmodel(tree, "curveball"), nsim = permutation.number)
    permuted.trees <- apply(permuted.trees, 3, function(x){as.data.frame(x)})
    print("finished Permuting")
    for (i in 1:permutation.number){
        print(paste("permutation #", i))
        #permuted.tree <- apply(tree, 2, function(x){(sample(x, size = sample.number, FALSE))})
        #permuted.tree <- permuted.trees$perm[[i]]
        permuted.tree <- permuted.trees[[i]]
        all.values <- as.numeric(OptimizedCalculateEntropy(permuted.tree))
        #permuted.values[i] <- max(as.numeric(OptimizedCalculateEntropy(permuted.tree)))
        permuted.values[i] <- max(all.values) #, na.rm = T)
        max.idx <- which(all.values == max(all.values)) #, na.rm = T))
        if (length(colnames(tree[max.idx]))==0){
          print("checking permute!")
          print(i)
          print(str(tree))
          print(colnames(tree))
          print(all.values)
          print(permuted.values[i])
          print(max.idx)
          print(colnames(tree[max.idx]))
          stop()
        }
        if (length(permuted.names[i]) != length(colnames(tree[max.idx]))){
          # print("lengths unequal!")
          # print(i)
          # print(all.values)
          # print(max.idx)
          # print(colnames(tree[max.idx]))
          max.idx <- max.idx[1]
        }
        permuted.names[i] <- colnames(tree[max.idx])
    }
    names(permuted.values) <- permuted.names
    permuted.values
    
}


DecisionTree <- function(tree, depth, parent, permutations){
    
    ## Takes input to decision tree analysis and recursively draws tree to a pre-specified depth
    
    ## Args:
    ## tree: data.frame with variables as columns
    ## depth: number of levels for tree
    ## parent: recursive argument to keep track of parent node. Initialize with tumor type 
    ## permutations: number of permutations to perform at each split. 0 to skip permutations testing 
    
    
    if (depth < 1){
        ## do nothing
    }else{
        
        ## Calculate highest informative feature     
        print(paste("Calculating values for level ", depth, sep = ""))
        print(dim(tree))
        #print(colnames(tree))
        output <- CalculateEntropy(tree)
        if (length(output)==1){
          #stop building tree
          print("stopping")
          return(NULL)
        }
        else{
          feature <- names(sort(output, decreasing = TRUE)[1])
        if (depth == 3 & feature == "10q_loss"){
          feature <- names(sort(output, decreasing = TRUE)[2])
          feature.value <- sort(output, decreasing = TRUE)[2]
        }
        else{
        feature <- names(sort(output, decreasing = TRUE)[1])
        feature.value <- sort(output, decreasing = TRUE)[1]
        }
        #print(sort(output, decreasing = TRUE))
        
        ## compare information gain to permuted information gain 
        if (permutations != 0){
            permuted.value <- max(branch1.permuted <- PermuteDecisionTree(tree, permutations))
        }else{
            permuted.value <- "NA"
        }
        print(paste("Information gain is ", feature.value, " permuted max is ", permuted.value, " for 5 permutations when splitting based on ", feature, sep = ""))    
        
        ## data structure to record tree relationships as well as top differential features at each branch
        #print("recording tree relationships")
        feature.column <- which(colnames(tree) == feature)
        #splits <- paste(table(tree[,feature.column])[1], table(tree[,feature.column])[2], sep = "-")
        current <- c(feature, parent, depth, "tbd")
        names(current) <- c("feature", "parent", "depth", "side")
        natural.history <- paste(parent, feature, sep = "--")
        print(natural.history)
        
        # print("finding enriched factors")
        #enriched.factors <- MaxFisherValue(tree, feature)
        enriched.factors <- MaxFisherValue(tree.or[rownames(tree), ], feature)
        enriched.list <- list(enriched.factors)
        names(enriched.list) <- natural.history
        #(splits
        
        ## recurisvely call on split data
        print("recursively calling")
        print(natural.history)
        print(dim(tree))
        print(feature)
        print(table(tree[,feature.column]))
        print(sum(is.na(tree[, feature.column])))
        print("left")
        print(dim(tree[tree[, feature.column] == 1 & !is.na(tree[ ,feature.column]), -feature.column]))
        print("right")
        # print(dim(tree[tree[, feature.column] == 0 & !is.na(tree[ ,feature.column]), -feature.column]))
        # print(colSums(tree[tree[, feature.column] == 0, -feature.column], na.rm = TRUE))
         print(dim(tree[tree[, feature.column] != 1 & !is.na(tree[ ,feature.column]), -feature.column]))
        if (table(tree[,feature.column])[1] > 2){
          print(colSums(tree[tree[, feature.column] != 1, -feature.column], na.rm = TRUE))
          right <- DecisionTree(tree[tree[, feature.column] != 1 & !is.na(tree[ ,feature.column]), -feature.column], depth - 1, natural.history, permutations)
        }
         else {
           print("less than 2 row in right branch")
           right <- NULL
         }
         print(dim(tree[tree[, feature.column] == 1 & !is.na(tree[ ,feature.column]), -feature.column]))
         if (table(tree[,feature.column])[2] > 2){
           print(colSums(tree[tree[, feature.column] == 1, -feature.column], na.rm = TRUE))
           left <- DecisionTree(tree[tree[, feature.column] == 1 & !is.na(tree[ ,feature.column]), -feature.column], depth - 1, natural.history, permutations)
         }
         else {
           print("less than 2 row in left branch")
           left <- NULL
         }
         # right <- DecisionTree(tree[tree[, feature.column] == 0 & !is.na(tree[ ,feature.column]), -feature.column], depth - 1, natural.history, permutations)
        
        
        ## return data
        if (is.null(right) & is.null(left)){
            ## don't return subsequent calls
            return (list(current, enriched.list))
        }else{
            print("left & right:")
            print(length(right))
            print(str(right))
            print(length(left))
            print(str(left))
            ## check for correct indexing
            if (!is.null(left)){
              if (length(left[[1]]) > 4){
                left[[1]][1, 4] <- "left"
              }else{
                left[[1]][4] <- "left"
              }
            }
            
            if (!is.null(right)){
            if (length(right[[1]]) > 4){
                right[[1]][1, 4] <- "right"
            }else{
                right[[1]][4] <- "right"
            }
            }
            
            combined.ordering <- rbind(current, left[[1]], right[[1]])
            combined.enriched <- c(left[[2]], right[[2]], enriched.list)
            return(list(combined.ordering, combined.enriched) )
        }
        }
    }
}



MaxFisherValue <- function(data, splitter){
    ## Determines the most significantly different features of a decision tree when split along a given column
    
    ## Args:
    
    ## data: a data frame
    ## splitter: column name for data to be split upon
    #print(dim(data))
    idx <- 1:ncol(data)
    # print(splitter)
    # print("idx")
    # print(length(idx))
    split.idx <- which(colnames(data) == splitter)
    skip.idx <- colSums(data, na.rm = TRUE) == 0 | colSums(data, na.rm = TRUE) == colSums(!is.na(data))
    #print(skip.idx)
    idx <- idx[!skip.idx]
    #print(length(idx))
    
    values <- data.frame(rep(0, length(idx)), rep(0, length(idx)), rep("name", length(idx)), stringsAsFactors = FALSE)
    colnames(values) <- c("OR", "p.value", "gene")
    print(splitter)
    print(dim(data))
    skips <- 0
    for (i in idx){
        x <- i-skips
        #print(i)
        #print(colnames(data)[i])
        # if (length(unique(data[,split.idx])) + length(unique(data[,i])) > 6){
        #   #print("using simulated p-val")
        #   result <- fisher.test(data[, split.idx], data[, i], simulate.p.value = TRUE)
        #   #print(result)
        #   values[i, 1] <- "N/A"
        #   values[i, 2] <- result$p.value
        #   values[i, 3] <- colnames(data)[i]
        # }
        #else{
        if (ncol(data)<40 & colnames(data)[i] == "PDGFRA"){
          #print(cbind(data[, split.idx], data[, i]))
          #print(levels(data[, split.idx]))
          #print(levels(data[, i]))
        }
        #print(split.idx)
        dt <- table(data[, split.idx], data[, i])
        #print(dt)
        if ((sum(rowSums(dt) ==0, colSums(dt) ==0) != 0) | nrow(dt)<2 | ncol(dt)<2){
          print(dim(data))
          print(colnames(data)[split.idx])
          print(colnames(data)[i])
          print(table(data[, split.idx], data[, i]))
          skips <- skips + 1
          print(paste("skipping OR calculations for:", colnames(data)[i]))
        }
        else{
        # print(colnames(data)[split.idx])
        # print(colnames(data)[i])  
        # print(table(data[, split.idx], data[, i]))
        result <- fisher.test(data[, split.idx], data[, i])
        values[x, 1] <- result$estimate
        values[x, 2] <- result$p.value
        values[x, 3] <- colnames(data)[i]
        #}
        }
    }
    if (sum(is.na(values$p.val))>0){
      #print(values)
      print(sum(is.na(values$p.val)))
      print(skips)
      #print(values[!(is.infinite(values$OR) | values$gene == "name" ), ])
      print(paste("NAs IN VALUES for node:", colnames(data)[split.idx]))
    }
    ## return top 2 by significance for >1 and <1 OR
    values <- values[values$gene != splitter, ]
    values$q.value <- p.adjust(values$p.value, "fdr")
    values <- values[order(values$q.value), ]
    #print(values)
    #values <- values[!(is.infinite(values$OR) | values$gene == "name" ), ]
    values <- values[!(values$gene == "name" ), ]
    #print(values)
    # positives <- which(values$OR > 1)[1:4]
    # negatives <- which(values$OR < 1)[1:4]
    positives <- which(values$OR > 1 & values$q.value < 0.1)
    negatives <- which(values$OR < 1 & values$q.value < 0.1)
    values[c(positives, negatives), ]
    #print("finished fisher")
}



