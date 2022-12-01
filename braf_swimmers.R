## Pinky Langat 2021
## This code generates a swimmers plot based on clinical treatment and outcomes data

#install.packages("plotly")
#install.packages("swimplot")
#install.packages("ggforce")

library(swimplot)
library(ggplot2) 
library(ggforce)
#library(dplyr, warn.conflicts=FALSE)   # Useful for manipulating the dataframes
#library(reshape2) # Reformmating dataframes
#ibrary(grid)
#library(plotly) # Allows us to make the swimmer plot interactive
#library(knitr)


## Load data files
braf_swimmers.tx  <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/outcomes/braf_swimmers_treatments.csv", stringsAsFactors = F)

braf_swimmers.obs  <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/outcomes/braf_swimmers_observation.csv", stringsAsFactors = F)

braf_swimmers.events <- read.csv("/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/outcomes/braf_swimmers_events.csv", stringsAsFactors = F)

## Format general treatment categories
braf_swimmers.tx$Treatment_Specific <- braf_swimmers.tx$Treatment
braf_swimmers.tx$Treatment <- "Other"
braf_swimmers.tx$Treatment[is.na(braf_swimmers.tx$Treatment_Specific)] <- NA
braf_swimmers.tx$Treatment[grepl("Observation", braf_swimmers.tx$Treatment_Specific)] <- "Observation"
braf_swimmers.tx$Treatment[grepl("BRAF", braf_swimmers.tx$Treatment_Specific)] <- "BRAFi/MEKi"
braf_swimmers.tx$Treatment[grepl("MEK", braf_swimmers.tx$Treatment_Specific)] <- "BRAFi/MEKi"
braf_swimmers.tx$Treatment[grepl("XRT", braf_swimmers.tx$Treatment_Specific)] <- "XRT+TMZ or XRT"
braf_swimmers.tx$Treatment[braf_swimmers.tx$Treatment_Specific == "TMZ"] <- "TMZ"

## Format ID
braf_swimmers.obs$ID <- braf_swimmers.tx$ID[match(braf_swimmers.obs$SAMPLE_ACCESSION_NBR, braf_swimmers.tx$SAMPLE_ACCESSION_NBR)]
braf_swimmers.events$ID <- braf_swimmers.tx$ID[match(braf_swimmers.events$SAMPLE_ACCESSION_NBR, braf_swimmers.tx$SAMPLE_ACCESSION_NBR)]

## Format preBRAF times
braf_swimmers.tx$End_Time_preBRAF <- -braf_swimmers.tx$Start_Time
braf_swimmers.obs$Start_Obs_preBRAF <- -braf_swimmers.obs$End_Obs_BRAF0
braf_swimmers.obs$End_Obs_preBRAF <- -braf_swimmers.obs$Start_Obs_BRAF0
braf_swimmers.events$Time_preBRAF <- -braf_swimmers.events$Time_BRAF0

## Filter for pre & post BRAF only
braf_swimmers_postBRAF.tx <- braf_swimmers.tx[braf_swimmers.tx$Start_Time > -0.01,]
braf_swimmers_postBRAF.obs <- braf_swimmers.obs[braf_swimmers.obs$Start_Obs_BRAF0 > -0.01,]
braf_swimmers_postBRAF.events <- braf_swimmers.events[braf_swimmers.events$Time_BRAF0 > -0.01,]

braf_swimmers_preBRAF.tx <- braf_swimmers.tx[braf_swimmers.tx$Start_Time < 0,]
braf_swimmers_preBRAF.obs <- braf_swimmers.obs[braf_swimmers.obs$Start_Obs_BRAF0 < 0,]
braf_swimmers_preBRAF.events <- braf_swimmers.events[braf_swimmers.events$Time_BRAF0 < 0,]

## Set ID order
ID_order <- rev(c("PXA-1","PXA-2","PXA-3","PXA-4","PXA-5","PXA-6","GBM-1","GBM-2","GBM-3","GBM-4*","PA-1**","PA-2","Astro-1"))

### POST BRAF
## 1a. Super basic look at samples you have and total treatment length
swimmer_plot(df=braf_swimmers_postBRAF.tx,id='ID',end='End_Time',fill='lightblue',width=0.85)

## 1b. Treatment plot with specified colors, obs line, and starting at zero for post-BRAF
tx_plot <- swimmer_plot(df=braf_swimmers_postBRAF.tx, id='ID', end='End_Time', name_fill='Treatment', col="white", size=0.25, alpha=0.95,width=0.5, id_order=ID_order) + scale_fill_manual(name="Treatment", values=c("BRAFi/MEKi"="#e41a1c", "XRT+TMZ or XRT"="#377eb8", "TMZ"="#4daf4a", "Other"="grey"), breaks=c("BRAFi/MEKi", "XRT+TMZ or XRT", "TMZ", "Other")) + swimmer_lines(df_lines=braf_swimmers_postBRAF.obs, id='ID', start='Start_Obs_BRAF0', end='End_Obs_BRAF0', col='black', size=0.5, linetype="dotted") + scale_y_continuous(name="Time (Months)",breaks=seq(0,96,by=6)) + theme(axis.text=element_text(size =12), axis.title=element_text(size=14)) + scale_x_discrete(position="top") #+ scale_y_log10(breaks=c(10, 20, 30, 40, 50, 75, 100, 200, 300)) #+ theme(aspect.ratio = 0.6) #+theme(panel.border=element_blank(), axis.line=element_line(colour="black"))
tx_plot
# grey #656565

## 2. Add events to plot
ev_plot <- tx_plot + swimmer_points(df_points=braf_swimmers_postBRAF.events,id='ID',time='Time_BRAF0',name_shape = 'Event',size=4,fill='white',col='black', adj.y=0) + scale_shape_manual(name="Event",values=c(Progression=21,Death=18,Remission=23),breaks=c('Progression','Remission','Death'))
ev_plot

## 3. Add continuation arrows
swim_plot <- ev_plot + swimmer_arrows(df_arrows=braf_swimmers_postBRAF.tx, id='ID', arrow_start = 'End_Time', cont= 'Continued_treatment', arrow_positions = c(0.1, 2), length=0.1, cex=0.5, type='open')
swim_plot

### PRE BRAF
## 1a. Super basic look at samples you have and total treatment length
swimmer_plot(df=braf_swimmers_preBRAF.tx,id='SAMPLE_ACCESSION_NBR',end='End_Time_preBRAF',fill='lightblue',width=0.85)

## 1b. Treatment plot with specified colors, obs line, and starting at zero for post-BRAF
braf_swimmers_preBRAF.tx$mask = 0
braf_swimmers_preBRAF.tx$mask[braf_swimmers_preBRAF.tx$End_Time_preBRAF > 100] = 1
braf_swimmers_preBRAF.obs$mask = 0
braf_swimmers_preBRAF.obs$mask[braf_swimmers_preBRAF.obs$Start_Obs_preBRAF > 100] = 1

tx_plot2 <- swimmer_plot(df=braf_swimmers_preBRAF.tx, id='ID', end='End_Time_preBRAF', name_fill='Treatment', col="black", size=0.25, alpha=0.95,width=0.5, id_order=ID_order) + scale_fill_manual(name="Treatment", values=c("BRAFi/MEKi"="#e41a1c", "XRT+TMZ or XRT"="#377eb8", "TMZ"="#4daf4a", "Other"="grey"), breaks=c("BRAFi/MEKi", "XRT+TMZ or XRT", "TMZ", "Other"))  + theme(axis.text=element_text(size =12), axis.title=element_text(size=14)) + scale_y_continuous(name="Time (Months)") + facet_zoom(ylim = c(0, 100))#+ facet_grid(. ~mask, scales="free", space="free")#,breaks=seq(0,96,by=6))  #+ scale_y_log10(breaks=c(10, 20, 30, 40, 50, 75, 100, 200, 300)) #+ theme(aspect.ratio = 0.6) #+theme(panel.border=element_blank(), axis.line=element_line(colour="black")) 
tx_plot2

tx_plot2 <- swimmer_plot(df=braf_swimmers_preBRAF.tx, id='ID', end='End_Time_preBRAF', name_fill='Treatment', col="white", size=0.25, alpha=0.95,width=0.5, id_order=ID_order) + scale_fill_manual(name="Treatment", values=c("BRAFi/MEKi"="#e41a1c", "XRT+TMZ or XRT"="#377eb8", "TMZ"="#4daf4a", "Other"="grey"), breaks=c("BRAFi/MEKi", "XRT+TMZ or XRT", "TMZ", "Other"))  + theme(axis.text=element_text(size =12), axis.title=element_text(size=14)) + scale_y_continuous(name="Time (Months)", limits = c(0,240), breaks = seq(0, 240, by = 24)) + swimmer_lines(df_lines=braf_swimmers_preBRAF.obs, id='ID', start='Start_Obs_preBRAF', end='End_Obs_preBRAF', col='black', size=0.5, linetype="dotted") + scale_x_discrete(position="top") #+scale_y_continuous(breaks=c(0,12,24,36,48,60,72,84,100,200,300,400))#+ facet_grid(. ~mask, scales="free") #,breaks=seq(0,96,by=6))  #+ scale_y_log10(breaks=c(10, 20, 30, 40, 50, 75, 100, 200, 300)) #+ theme(aspect.ratio = 0.6) #+theme(panel.border=element_blank(), axis.line=element_line(colour="black")) 
tx_plot2
# grey #656565
#

## 2. Add events to plot
ev_plot2 <- tx_plot2 + swimmer_points(df_points=braf_swimmers_preBRAF.events,id='ID',time='Time_preBRAF',name_shape = 'Event',size=4,fill='white',col='black', adj.y=0) + scale_shape_manual(name="Event",values=c(Progression=21,Death=18,Remission=23),breaks=c('Progression','Remission','Death'))
ev_plot2

## 3. Add continuation arrows
swim_plot2 <- ev_plot2 + swimmer_arrows(df_arrows=braf_swimmers_preBRAF.tx, id='ID', arrow_start = 'End_Time_preBRAF', cont= 'Continued_treatment', arrow_positions = c(0.1, 2), length=0.1, cex=0.5, type='open') + swimmer_lines(df_lines=braf_swimmers_preBRAF.obs, id='ID', start='Start_Obs_preBRAF', end='End_Obs_preBRAF', col='black', size=0.5, linetype="dotted") + swimmer_lines(df_lines=braf_swimmers_preBRAF.obs, id='ID', start='Start_Obs_preBRAF', end='End_Obs_preBRAF', col='black', size=0.5, linetype="dotted") + facet_grid(braf_swimmers_preBRAF.obs ~ mask, scales="free", space="free")
swim_plot2




####
### Building it up/testing features of swimmers plot ###
## a. Treatment with default colors
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.75,width=.75)
tx_plot

## b. Treatment plot with set colors
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.75,width=.75) + scale_fill_brewer(type = "qual", palette = "Paired")
tx_plot

## c. Treatment plot with set colors, stratified by path subtype
tx_plot <- swimmer_plot(df=braf_swimmers_postBRAF.tx,id='ID',end='End_Trt_BRAF0',name_fill='Treatment', col="black",alpha=0.75,width=.75, stratify='Subtype') + scale_fill_brewer(type = "qual", palette = "Paired")
tx_plot

## d. Treatment plot with set colors and obs line
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.75,width=.75) + scale_fill_brewer(type = "qual", palette = "Paired") + swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs', end='End_Obs', name_col='Observation', size=.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey"))
tx_plot

## e. Treatment plot with specific colors and obs line
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.9,width=.75) + swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs', end='End_Obs', name_col='Observation', size=.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#9e0142", "BRAF"="#d53e4f","MEK"="#f46d43",  "XRT+TMZ"="#fee08b", "TMZ"="#ffffbf", "XRT"="#e6f598","VEGF"="#abdda4", "PD-1"="#66c2a5", "PD-1, VEGF"="#66c2a5", "PD-1, XRT+TMZ"="#66c2a5", "ATM, XRT"="#3288bd", "Oncolytic virus"="#5e4fa2",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "TMZ", "XRT", "VEGF", "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "ATM, XRT", "Oncolytic virus"))
tx_plot

## f. Treatment plot with specific colors #2 and obs line
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.9,width=.75) + swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs', end='End_Obs', name_col='Observation', size=.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#d53e4f", "BRAF"="#f46d43","MEK"="#fdae61",  "XRT+TMZ"="#ffffbf", "TMZ"="#e6f598", "XRT"="#abdda4","VEGF"="#3288bd", "PD-1"="#5e4fa2", "PD-1, VEGF"="#5e4fa2", "PD-1, XRT+TMZ"="#5e4fa2", "ATM, XRT"="grey", "Oncolytic virus"="grey",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "TMZ", "XRT", "VEGF", "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "ATM, XRT", "Oncolytic virus"))
tx_plot

## g. Treatment plot with specific colors #3 and obs line
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.9,width=.75) + swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs', end='End_Obs', name_col='Observation', size=.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#e31a1c", "BRAF"="#ff7f00","MEK"="#fdbf6f",  "XRT+TMZ"="#ffff99", "XRT"="#33a02c", "TMZ"="#b2df8a", "PD-1"="#a6cee3", "PD-1, VEGF"="#a6cee3", "PD-1, XRT+TMZ"="#a6cee3", "VEGF"="#cab2d6", "ATM, XRT"="grey", "Oncolytic virus"="grey",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "XRT", "TMZ", "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "VEGF", "ATM, XRT", "Oncolytic virus"))
tx_plot

## h. Treatment plot with specific colors #4 and obs line
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.9,width=.75) + swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs', end='End_Obs', name_col='Observation', size=.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#c51b7d", "BRAF"="#de77ae","MEK"="#f1b6da",  "XRT+TMZ"="#e6f598", "XRT"="#abdda4", "TMZ"="#80cdc1",  "PD-1"="#abd9e9", "PD-1, VEGF"="#abd9e9", "PD-1, XRT+TMZ"="#abd9e9", "VEGF"="#b2abd2", "ATM, XRT"="grey", "Oncolytic virus"="grey",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "XRT", "TMZ",  "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "VEGF", "ATM, XRT", "Oncolytic virus"))
tx_plot

## i. Treatment plot with specific colors #5 and obs line
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.9,width=.75) + swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs', end='End_Obs', name_col='Observation', size=.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#c51b7d", "BRAF"="#fc4e2a","MEK"="#fd8d3c",  "XRT+TMZ"="#e6f598", "XRT"="#abdda4", "TMZ"="#80cdc1",  "PD-1"="#abd9e9", "PD-1, VEGF"="#abd9e9", "PD-1, XRT+TMZ"="#abd9e9", "VEGF"="#b2abd2", "ATM, XRT"="grey", "Oncolytic virus"="grey",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "XRT", "TMZ",  "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "VEGF", "ATM, XRT", "Oncolytic virus"))
tx_plot

## j. Treatment plot with specific colors #6 and obs line
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.9,width=.75) + swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs', end='End_Obs', name_col='Observation', size=.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#ce1256", "BRAF"="#e7298a","MEK"="#df65b0",  "XRT+TMZ"="#e6f598", "XRT"="#abdda4", "TMZ"="#80cdc1",  "PD-1"="#abd9e9", "PD-1, VEGF"="#abd9e9", "PD-1, XRT+TMZ"="#abd9e9", "VEGF"="#74add1", "ATM, XRT"="grey", "Oncolytic virus"="grey",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "XRT", "TMZ",  "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "VEGF", "ATM, XRT", "Oncolytic virus"))
tx_plot

## k. Treatment plot with specific colors #6 and stratified starting at BRAF tx
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="lightgrey",alpha=0.9,width=.75, stratify='Subtype') + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#ce1256", "BRAF"="#e7298a","MEK"="#df65b0",  "XRT+TMZ"="#e6f598", "XRT"="#abdda4", "TMZ"="#80cdc1",  "PD-1"="#abd9e9", "PD-1, VEGF"="#abd9e9", "PD-1, XRT+TMZ"="#abd9e9", "VEGF"="#74add1", "ATM, XRT"="grey", "Oncolytic virus"="grey",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "XRT", "TMZ",  "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "VEGF", "ATM, XRT", "Oncolytic virus"))
tx_plot

## l. Treatment plot with specific colors #6 and obs line and starting at BRAF tx
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt_BRAF0',name_fill='Treatment', col="lightgrey",alpha=0.9,width=.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#ce1256", "BRAF"="#e7298a","MEK"="#df65b0",  "XRT+TMZ"="#e6f598", "XRT"="#abdda4", "TMZ"="#80cdc1",  "PD-1"="#abd9e9", "PD-1, VEGF"="#abd9e9", "PD-1, XRT+TMZ"="#abd9e9", "VEGF"="#74add1", "ATM, XRT"="grey", "Oncolytic virus"="grey",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "XRT", "TMZ",  "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "VEGF", "ATM, XRT", "Oncolytic virus")) + swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs_BRAF0', end='End_Obs_BRAF0', name_col='Observation', size=.75)
tx_plot

## m. Treatment plot with specific colors #6 and obs line and starting at BRAF tx
tx_plot <- swimmer_plot(df=braf_swimmers.tx,id='ID', start='start', end='end', name_fill='Treatment', col="grey",alpha=0.9,width=0.75) + scale_color_manual(name="Observation",values=c("Observation"="darkgrey")) + scale_fill_manual(name="Treatment", values=c("BRAF+MEK"="#ce1256", "BRAF"="#e7298a","MEK"="#df65b0",  "XRT+TMZ"="#e6f598", "XRT"="#abdda4", "TMZ"="#80cdc1",  "PD-1"="#abd9e9", "PD-1, VEGF"="#abd9e9", "PD-1, XRT+TMZ"="#abd9e9", "VEGF"="#74add1", "ATM, XRT"="grey", "Oncolytic virus"="grey",na.value=NA),breaks=c("BRAF+MEK","BRAF","MEK",  "XRT+TMZ", "XRT", "TMZ",  "PD-1", "PD-1, VEGF", "PD-1, XRT+TMZ", "VEGF", "ATM, XRT", "Oncolytic virus")) #+ swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs_BRAF0', end='End_Obs_BRAF0', name_col='Observation', size=.75)
tx_plot

## n. Treatment plot with specific colors and obs line and starting at BRAF tx
tx_plot <- swimmer_plot(df=braf_swimmers.tx, id='ID', end='end_time', name_fill='Treatment', col="white", alpha=0.95,width=0.5) + scale_fill_manual(name="Treatment", values=c("BRAFi/MEKi"="#e41a1c", "XRT+TMZ or XRT"="#377eb8", "TMZ"="#4daf4a", "Other"="grey", na.value='lightgrey'),breaks=c("BRAFi/MEKi", "XRT+TMZ or XRT", "TMZ", "Other"))  #+ swimmer_lines(df_lines=braf_swimmers.obs, id='ID', start='Start_Obs', end='End_Obs', col='#656565', size=.5, linetype="dotted") #+ scale_y_log10(breaks=c(10, 20, 30, 40, 50, 75, 100, 200, 300)) #+ theme(aspect.ratio = 0.6)
tx_plot

### Adding events
## a. Basic events
ev_plot <- tx_plot + swimmer_points(df_points=braf_swimmers.events,id='ID',time='Time',name_shape = 'Event',size=2.5,fill='black',col='black')
ev_plot

## b. Basic events with manual shapes
ev_plot <- tx_plot + swimmer_points(df_points=braf_swimmers.events,id='ID',time='Time',name_shape = 'Event',size=2.5,fill='black',col='black') + scale_shape_manual(name="Event",values=c(Progression=2,Death=23,'In Remission'=1),breaks=c('Progression','Death','In Remission'))
ev_plot

## c. Basic events with manual shapes slightly below treatments
ev_plot <- tx_plot + swimmer_points(df_points=braf_swimmers.events,id='ID',time='Time',name_shape = 'Event',size=2,fill='grey35',col='grey35', adj.y=-0.15) + scale_shape_manual(name="Event",values=c(Progression=2,Death=23,'In Remission'=1),breaks=c('Progression','Death','In Remission'))
ev_plot

## d. Basic events with manual shapes slightly below treatments start at BRAF
ev_plot <- tx_plot + swimmer_points(df_points=braf_swimmers.events,id='ID',time='Time',name_shape = 'Event',size=2,fill='black',col='black', adj.y=0) + scale_shape_manual(name="Event",values=c(Progression=1,Death=23,'In Remission'=5),breaks=c('Progression','In Remission','Death'))
ev_plot

## Adding arrows
ev_plot+
  swimmer_arrows(df_arrows=braf_swimmers.tx,id='ID',arrow_start='End_Trt',
                 cont = 'Continued_treatment',name_col='Treatment',type =
                   "open",cex=1)

swim_plot <- ev_plot +
  swimmer_arrows(df_arrows=braf_swimmers.tx,id='ID',arrow_start='End_Trt',
                 cont = 'Continued_treatment',name_col='Treatment',show.legend = FALSE,type =
                   "open",cex=1) + scale_color_discrete(drop=FALSE)
swim_plot <- swim_plot +
  #scale_fill_manual(name="Treatment",values=c("Arm A" = "#e41a1c", "Arm B"="#377eb8","Off Treatment"='#4daf4a'))+
  #scale_color_manual(name="Treatment",values=c("Arm A"="#e41a1c", "Arm B" ="#377eb8","Off Treatment"='#4daf4a')) +
  scale_shape_manual(name="Event",values=c(Progression=24,Death=17,'In Remission'=19),breaks=c('Progression','Death','In Remission'))
swim_plot

## Stratifying 
swim_plot_stratify <- swimmer_plot(df=braf_swimmers.tx,id='ID',end='End_Trt',name_fill='Treatment', col="white",alpha=0.75,width=.8,base_size = 10,stratify= c('Subtype')) + coord_flip()
swim_plot_stratify

## Troubleshooting the gaps in events
Gap_data <- data.frame(patient_ID=c('ID:3','ID:1','ID:1','ID:1','ID:2',
                                    'ID:2','ID:2','ID:3','ID:3','ID:1'),
                       start=c(10,1,2,7,2,10,14,5,0,0),
                       end=c(20,2,4,10,7,14,22,7,3,-2),
                       treatment=c("A","B","C","A","A","C","A","B","C","C"))

swimmer_plot(df=Gap_data,id='patient_ID',name_fill="treatment",col='white', width=0.25) + scale_y_continuous(breaks=c(-5:25))

### Space to check things quicks
table(all.braf.mastersheet$Subtype)

tcga.braf.mastersheet <- all.braf.mastersheet[grepl("TCGA", all.braf.mastersheet$SAMPLE_ACCESSION_NBR),]
write.csv(tcga.braf.mastersheet, "/Users/pclangat/Dropbox (Partners HealthCare)/Lab projects current/Basic Science/BRAF glioma/3-Data/rna-seq/tcga_braf_glioma_mastersheet.csv")
