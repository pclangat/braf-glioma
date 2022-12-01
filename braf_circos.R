#install.packages("circlize")

library(circlize)
library(RColorBrewer)


#### Test example code
## Set colors
col = c('#E41A1C', '#A73C52', '#6B5F88', '#3780B3', '#3F918C', '#47A266',
        '#53A651', '#6D8470', '#87638F', '#A5548D', '#C96555', '#ED761C',
        '#FF9508', '#FFC11A', '#FFEE2C', '#EBDA30', '#CC9F2C', '#AD6428',
        '#BB614F', '#D77083', '#F37FB8', '#DA88B3', '#B990A6', '#999999')

## Initialize (blank)
circos.initializeWithIdeogram(plotType = NULL)

## Chr labels
circos.trackPlotRegion(ylim = c(0, 1), bg.border = NA, track.height = 0.05,
                       panel.fun = function(x, y) {
                         xlim = get.cell.meta.data("xlim")
                         chr = get.cell.meta.data("sector.index")
                         circos.text(mean(xlim), 0.5, labels = chr, facing = "clockwise", niceFacing = TRUE)
                       })

## Random data generated
bed = generateRandomBed(nr = 200, fun = function(k) runif(k))

## Plot genomic ?coverage by chromosome
circos.genomicTrackPlotRegion(bed, bg.border = NA, panel.fun = function(region, value, ...) {
  i = get.cell.meta.data("sector.numeric.index")
  circos.genomicLines(region, value, area = TRUE, border = NA, baseline = 0, col = col[i])
}, track.height = 0.1)

## Plot color blocks for each chromosome
circos.trackPlotRegion(ylim = c(0, 1), bg.border = col, bg.col = col, panel.fun = function(x, y) {
  chr = get.cell.meta.data("sector.index")
  circos.axis(h = "bottom", labels = NULL, sector.index = chr, direction = "inside")
}, track.height = 0.05)

## Define fusion regions and link on circos
region1 = data.frame(chr = c("chr1", "chr1"),
                     start = c(12345678, 22222222),
                     #end = c(12345678, 22222222))
                     end = c(12345678, 42222222))
region2 = data.frame(chr = c("chr1", "chr1"),
                     start = c(87654321, 99999999),
                     end = c(87654321, 99999999))  

circos.genomicLink(region1, region2, h = 0.2)

#### End of test

## Define fusion regions and link on circos for CHR7 for ADULT GLIOMA
c("KIAA1549", "FAM131B", "AGK", "UBE2H","GNAI1")
# peds c("BCAS1", "CCDC6", "GIT2", "PTPRZ11")
region1 = data.frame(chr = c("chr7", "chr7", "chr7", "chr7", "chr7"),
                     start = c(138831376, 143353399, 141551409, 129830731, 80134830),
                     end = c(138981625, 143382303, 141662152, 129952959, 80226180))
region_braf = data.frame(chr = c("chr7", "chr7", "chr7", "chr7", "chr7"),
                     start = c(140713327,140713327,140713327,140713327,140713327),
                     end = c(140924928,140924928,140924928,140924928,140924928)) 

#circos.initializeWithIdeogram(plotType = NULL)
circos.clear()

## Settings for inset of figure (I think)
par(mar = c(1, 1, 1, 1), new = TRUE)

## Generate chromosome arch
#circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2), clock.wise = FALSE,
#circos.par("canvas.xlim" = c(-1, 1), "canvas.ylim" = c(-1, 1), clock.wise = FALSE,
           #cell.padding = c(0, 0, 0, 0), gap.degree = 180)
circos.par(clock.wise = FALSE, cell.padding = c(0, 0, 0, 0), gap.degree = 180)
circos.initializeWithIdeogram(chromosome.index = "chr7", plotType = c("ideogram", "axis"))
#circos.info()
text(0, 0.6, "chr7")

#circos.genomicLink(region1, region_braf, h = 0.5)
#circos.genomicLink(region1, region_braf, col = c(brewer.pal(4, "Set1")), lwd = c(9,3,2,1), h=c(0.5,0.2,0.2,0.2))
#circos.genomicLink(region1, region_braf, col = c(brewer.pal(4, "Set3")))
#circos.genomicLink(region1, region_braf, col = c(brewer.pal(5, "Set3")), lwd=c(8,3,3,1,1))
#circos.genomicText(region1, value=NULL, y = 3, labels=c("KIAA1549", "FAM131B", "AGK", "UBE2H","GNAI1"),col = c(brewer.pal(5, "Set3")), cex=0.6)
circos.genomicLink(region1, region_braf, col = c(brewer.pal(6, "Paired")[c(6,4,2,3,1)]), lwd=c(8,3.5,3.5,2,2))

circos.clear()

## Define fusion regions and link on circos for CHR7-20 for PED GLIOMA
#adult c("KIAA1549", "FAM131B", "AGK", "UBE2H","GNAI1")
c("KIAA1549", "PTPRZ11", "CCDC6","GIT2", "BCAS1")
region1.peds = data.frame(chr = c("chr7", "chr7", "chr10", "chr12",  "chr20"),
                     start = c(138831376, 121873160, 59788746, 109929791, 53943540),
                     end = c(138981625, 122062035, 59906555, 109996367, 54070593))
region_braf.peds = data.frame(chr = c("chr7", "chr7", "chr7", "chr7", "chr7"),
                         start = c(140713327,140713327,140713327,140713327,140713327),
                         end = c(140924928,140924928,140924928,140924928,140924928)) 

#circos.par(clock.wise = FALSE, cell.padding = c(0, 0, 0, 0))
#circos.initializeWithIdeogram(chromosome.index = "chr7", plotType = c("ideogram", "axis"))
#circos.initializeWithIdeogram(plotType = c("ideogram", "axis"))
#circos.initializeWithIdeogram()
circos.par(gap.degree = 2)
circos.initializeWithIdeogram(chromosome.index = paste0("chr", c(20:1,'X','Y',22,21)))
#text(0, 0.6, "chr7")
#circos.genomicText(region1.peds, value=NULL, y = 1, labels=c("KIAA1549", "PTPRZ11", "CCDC6","GIT2", "BCAS1"),col = c(brewer.pal(5, "Set3")), cex=0.6)
#circos.genomicLink(region1.peds, region_braf.peds, col = c(brewer.pal(6, "Set3")[c(1,3,4,5,6)]), lwd=c(10,2,2,2,2))
circos.genomicLink(region1.peds, region_braf.peds, col = c(brewer.pal(6, "Paired")[c(6,4,3,2,1)]), lwd=c(10,2,2,2,2))
#circos.genomicLink(region1.peds, region_braf.peds, col = c('#e31a1c', '#1f78b4', '#a6cee3',))

circos.clear()
