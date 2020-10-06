#### Clearing workspace, loading libraries, setting seed ----

## Clear environment and set working directory
rm(list = ls())
setwd("/media/projects1/Fabian/Oral microbiome/StrainRecognitionFC")

## Load libraries
library("Phenoflow")
library("flowViz")
library("flowFDA")
library("flowAI")
#library("vegan")
library("ggplot2")
library("RColorBrewer")
#library("ggrepel")
#library("ape")
#library("gridExtra")
#library("grid")
#library("scales")
library("cowplot")
library("reshape2")
library("dplyr")
library("tidyverse")

set.seed(777)


#### Loading data ----

Datapath <- "/media/projects1/Fabian/Oral microbiome/StrainRecognitionFC/data"
fcsfiles <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = TRUE)
flowData <- flowCore::read.flowSet(files = fcsfiles, transformation = FALSE, emptyValue = F)


#### Transformation of data ----
attributes(flowData)

flowData_transformed <- transform(flowData,
                                  `FSC-A` = asinh(`FSC-A`),
                                  `SSC-A` = asinh(`SSC-A`),
                                  `BL1-A` = asinh(`BL1-A`),
                                  `BL3-A` = asinh(`BL3-A`),
                                  `FSC-H` = asinh(`FSC-H`),
                                  `SSC-H` = asinh(`SSC-H`),
                                  `BL1-H` = asinh(`BL1-H`),
                                  `BL3-H` = asinh(`BL3-H`),
                                  `FSC-W` = asinh(`FSC-W`),
                                  `SSC-W` = asinh(`SSC-W`),
                                  `BL1-W` = asinh(`BL1-W`),
                                  `BL3-W` = asinh(`BL3-W`))
param = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-W", "SSC-W", "BL1-W", "BL3-W")


#### Extrating metadata from .fcs files ----
#sampleNames(flowData_transformed) <- substring(sampleNames(flowData_transformed), 0, nchar(sampleNames(flowData_transformed))-4)       # nchar takes a character vector and returns the number of characters in the vector

# Extracting the actual metadata
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed), "_"), rbind)))
colnames(metadata) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain", "Replicate", "Dilution", "Stain", "Well")

metadata$Dilution <- as.numeric(metadata$Dilution)

name <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = FALSE)
Sample_Info <- cbind(name, metadata)

#### Quality control data ----

## Contamination
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1:90)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 12-05-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(91:122)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 Mixes 12-05-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(123:164)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 Co-cultures 13-05-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(165:206)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 Co-cultures 14-05-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(207:296)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 14-05-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(297:386)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC MM 25-05-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(387:413)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC MM Mixes 25-05-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(414:503)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC MM 27-05-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(504:537)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 16-06-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(538:573)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC MM 16-06-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(574:613)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 Co-cultures 17-06-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(614:655)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC MM Co-cultures 17-06-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(656:675)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 Co-cultures 18-06-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(676:696)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC MM Co-cultures 18-06-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(697:729)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC MM 18-06-2020 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

## Singlets
p_singlets <- xyplot(`BL1-H`~`BL1-A`, data = flowData_transformed[c(1:729)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0,16))),
       axix = axis.default, nbin = 125, main = "QC singlets (BL1-A ~ BL1-H)", xlab = "BL1-A", ylab = "BL1-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)


#### Select data ----
flowData_transformed_sel <- flowData_transformed[c(1:122, 297:413, 656:696)]
metadata_sel <- metadata[c(1:122, 297:413, 656:696), ]

#### Gating ----

### BHI2
## Gating out extra background based on SSC-BL1
## Since the difference between control and filtered sample is not very encouraging, reducing BG through SSC-BL1 will be skipped
# Constructing gate
# sqrcut1 <- matrix(c(5, 12, 14, 10, 5,
#                     2, 5, 13.5, 13.5, 8.5), ncol = 2, nrow = 5)
# colnames(sqrcut1) <- c("BL1-A", "SSC-A")
# polyGate1 <- polygonGate(.gate = sqrcut1, filterId = "Reduced background BL1-SSC BHI2")
# 
# # Gating quality check
# xyplot(`SSC-A`~`BL1-A`, data = flowData_transformed_sel[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], filter = polyGate1,
#        scales = list(y = list(limits = c(0, 16)),
#                      x = list(limits = c(0, 16))),
#        axis = axis.default, nbin = 125, main = "Quality check gating BL1-SSC BHI2", xlab = "BL1-A", ylab = "SSC-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

## Gating cells based on BL1-BL3
sqrcut2 <- matrix(c(6.7, 11, 14, 10, 6.7,
                    2, 4.5, 13, 10.5, 5.3), ncol = 2, nrow = 5)
colnames(sqrcut2) <- c("BL1-A", "BL3-A")
polyGate2 <- polygonGate(.gate = sqrcut2, filterId = "Cells BHI2")

# Gating quality check pure cultures
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], filter = polyGate2,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check mixes
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119)], filter = polyGate2,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2 mixes", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check co-cultures
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(656, 659, 662, 665, 668, 671, 674, 675)], filter = polyGate2,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2 co-cultures", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

### MM
## Gating out extra background based on SSC-BL1
## Since reducing BG for BHI2 through SSC-BL1 was not performed, this will also not be done for MM
# # Constructing gate
# sqrcut3 <- matrix(c(5.5, 13, 15, 11.2, 5.5,
#                     2, 5, 14.5, 14.5, 9.6), ncol = 2, nrow = 5)
# colnames(sqrcut3) <- c("BL1-A", "SSC-A")
# polyGate3 <- polygonGate(.gate = sqrcut3, filterId = "Reduced background BL1-SSC MM")
# 
# # Gating quality check
# xyplot(`SSC-A`~`BL1-A`, data = flowData_transformed[c(297, 303, 309, 315, 317, 319, 325, 327, 333, 339, 345, 351, 357, 363, 369, 375, 381)], filter = polyGate3,
#        scales = list(y = list(limits = c(0, 16)),
#                      x = list(limits = c(4, 16))),
#        axis = axis.default, nbin = 125, main = "Quality check gating BL1-SSC MM", xlab = "BL1-A", ylab = "SSC-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

## Gating cells based on BL1-BL3
sqrcut4 <- matrix(c(7, 13.5, 14.8, 10, 7,
                    2, 5, 13.5, 10.5, 5.3), ncol = 2, nrow = 5)
colnames(sqrcut4) <- c("BL1-A", "BL3-A")
polyGate4 <- polygonGate(.gate = sqrcut4, filterId = "Cells MM")

# Gating quality check
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(297, 303, 309, 315, 317, 319, 325, 327, 333, 339, 345, 351, 357, 363, 369, 375, 381, 413)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check mixes
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(387, 389, 391, 393, 395, 397, 399, 401, 403, 405, 407, 409, 411)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM mixes", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check co-cultures
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(678, 680, 683, 686, 689, 692, 695, 696)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM co-cultures", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Check extra population for So
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(357:362)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM So", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)
# Sample A has an extra population --> leave it out for further analysis

# Check BG for Pi
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(333:338)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM Pi", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)
# Sample A looks a bit off, additionally, the duplicate measurement for sample A is off --> leave it out for further analysis
# Sample B also has a big difference between replicatre measurements, but better compared to sample A --> leave B out

### General gate for both media
## Gating cells based on BL1-BL3
sqrcut5 <- matrix(c(6.7, 13.3, 14.5, 10, 6.7,
                    2, 5, 13.5, 10.5, 5.3), ncol = 2, nrow = 5)
colnames(sqrcut5) <- c("BL1-A", "BL3-A")
polyGate5 <- polygonGate(.gate = sqrcut5, filterId = "Cells")

# Gating quality check pure cultures BHI2
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check mixes BHI2
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2 mixes", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check co-cultures BHI2
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(656, 659, 662, 665, 668, 671, 674, 675)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2 co-cultures", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check MM
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(297, 303, 309, 315, 317, 319, 325, 327, 333, 339, 345, 351, 357, 363, 369, 375, 381, 413)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check mixes MM
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(387, 389, 391, 393, 395, 397, 399, 401, 403, 405, 407, 409, 411)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM mixes", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check co-cultures MM
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(678, 680, 683, 686, 689, 692, 695, 696)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM co-cultures", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

## Gating quality check all samples
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_sel, filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

### Subsetting dataset
flowData_transformed_gated <- Subset(flowData_transformed, polyGate5)


#### Cell concentrations ----

### Cell counts
cells <- flowCore::filter(flowData_transformed, polyGate5)
TotalCount <- summary(cells)
TotalCount <- toTable(TotalCount)

### Extracting volumes
# Volumes are in µL
vol <- as.numeric(flowCore::fsApply(flowData_transformed, FUN = function(x) x@description$`$VOL`))/1000

### Concentrations
# Calculating concentrations --> Concentrations will be in cells/µL
cell_concentrations <- data.frame(Sample_name = flowCore::sampleNames(flowData_transformed),
                                  Strain = metadata$Strain,
                                  Concentration = (TotalCount$true*metadata$Dilution)/vol)


#### Phenotypic diversity analysis ----

### Normalization of data
summary <- fsApply(x = flowData_transformed_gated, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max = max(summary[, "BL1-A"])
max2 = max(summary[, "BL1-H"])
max3 = max(summary[, "BL1-W"])
mytrans <- function(x) x/max
mytrans2 <- function(x) x/max2
mytrans3 <- function(x) x/max3
flowData_transformed_2 <- transform(flowData_transformed_gated,
                                    `FSC-A` = mytrans(`FSC-A`), 
                                    `SSC-A` = mytrans(`SSC-A`), 
                                    `BL1-A` = mytrans(`BL1-A`), 
                                    `BL3-A` = mytrans(`BL3-A`),
                                    `FSC-H` = mytrans2(`FSC-H`), 
                                    `SSC-H` = mytrans2(`SSC-H`), 
                                    `BL1-H` = mytrans2(`BL1-H`), 
                                    `BL3-H` = mytrans2(`BL3-H`),
                                    `FSC-W` = mytrans3(`FSC-W`), 
                                    `SSC-W` = mytrans3(`SSC-W`), 
                                    `BL1-W` = mytrans3(`BL1-W`), 
                                    `BL3-W` = mytrans3(`BL3-W`))

### Calculating fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed_2, param, nbin = 128, 
                    bw = 0.01, normalize = function(x) x)

### Calculate ecological parameters
Diversity.fbasis <- Diversity(fbasis, d = 3, plot = FALSE, R = 999)
#Evenness.fbasis <- Evenness(fbasis, d = 3, plot = FALSE)
#Structural.organization.fbasis <- So(fbasis, d = 3, plot = FALSE)
#Coef.var.fbasis <- CV(fbasis, d = 3, plot = FALSE)

######## Script for plot is not updated yet
## Plot ecological parameters
palphadiv <- ggplot(data = Diversity.fbasis, aes(x = as.character(metadata$Concentration), y = D2))+
        geom_point(size = 4, alpha = 0.7)+
        geom_line()+
        facet_grid(metadata$Compound ~ ., scales = "free")+       # Creates different subplots
        theme_bw()+
        labs(y = "Phenotypic diversity (D2)", x = "Concentration", title = "Phenotypic diversity analysis")+
        geom_errorbar(aes(ymin = D2-sd.D2, ymax = D2+sd.D2), width = 0.05)
print(palphadiv)

### Export ecological data to .csv file in the chosen directory
#write.csv2(file="results.metrics.csv",
#           cbind(Diversity.fbasis, Evenness.fbasis,
#                 Structural.organization.fbasis,
#                 Coef.var.fbasis))

### Beta-diversity assessment of fingerprint
beta.div <- beta_div_fcm(fbasis, ord.type="PCoA")

######## Script for plot is not updated yet
# Plot ordination
plot_beta_fcm(beta.div, color = metadata$State, shape = as.factor(metadata$Time), labels = list("State", "Timepoint")) + 
        theme_bw() +
        geom_point(size = 8, alpha = 0.5)

#### Training Random Forest classifier ----

# Define parameters on which RF will build its model
paramRF = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-W", "SSC-W", "BL1-W", "BL3-W")

# Extract total number of events per fcs file in order to calculate accuracy of model to predict strain
vol2 <- data.frame(Sample_name = flowCore::sampleNames(flowData_transformed),
                   Volume = vol)
TotalCount <- left_join(TotalCount, vol2, by = c("sample" = "Sample_name"))
TotalCount <- left_join(TotalCount, cell_concentrations, by = c("sample" = "Sample_name"))
write.csv2(file = "TotalCount.csv", TotalCount)

### Model for So, Fn and Pg in BHI2
# Sample selection So, Fn and Pg
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(23:28, 319:324)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Fn", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(61:66, 357:362)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check So", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(31:36, 327:332)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Pg", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> Fn (replicate C), So (replicate A), Pg (replicate B) grown in BHI2
fcs_names_SoFnPg <- c("20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs")

Sample_Info_SoFnPg <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPg)

Model_RF_SoFnPg <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFnPg], Sample_Info_SoFnPg, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
# Test whether normalization of data influences performance of the model
#Model_RF_SoFnPg_N <- Phenoflow::RandomF_FCS(flowData_transformed_2[fcs_names_SoFnPg], Sample_Info_SoFnPg, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
# -> Performance is identical

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_SoFnPg <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_B1.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co1_B_1000_SG_C1.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co1_B_1000_SG_D1.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co1_C_1000_SG_E1.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co1_C_1000_SG_F1.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co1_A_1000_SG_A1.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co1_B_1000_SG_B1.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co1_C_1000_SG_C1.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co2_A_1000_SG_A2.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co2_A_1000_SG_B2.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co2_B_1000_SG_C2.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co2_B_1000_SG_D2.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co2_C_1000_SG_G1.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co2_C_1000_SG_H1.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co2_A_1000_SG_A2.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co2_B_1000_SG_B2.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co2_C_1000_SG_D1.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co3_A_1000_SG_A3.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co3_A_1000_SG_B3.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co3_B_1000_SG_E2.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co3_B_1000_SG_F2.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co3_C_1000_SG_G2.fcs",
                           "20200617_Fabian_14strainID_BHI2_24h_Co3_C_1000_SG_H2.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co3_A_1000_SG_A3.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co3_B_1000_SG_C2.fcs",
                           "20200618_Fabian_14strainID_BHI2_48h_Co3_C_1000_SG_D2.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix1_NA_1000_SG_A8.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix1_NA_1000_SG_A9.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix2_NA_1000_SG_A10.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix2_NA_1000_SG_A11.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix3_NA_1000_SG_B8.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix3_NA_1000_SG_B9.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix8_NA_1000_SG_D10.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix8_NA_1000_SG_D11.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix9_NA_1000_SG_E8.fcs",
                           "20200512_Fabian_14strainID_BHI2_NA_Mix9_NA_1000_SG_E9.fcs")

flowData_topre_BHI2_SoFnPg <- flowData_transformed_gated[fcs_topre_BHI2_SoFnPg]
test_pred_BHI2_SoFnPg <- RandomF_predict(x = Model_RF_SoFnPg[[1]], new_data =  flowData_topre_BHI2_SoFnPg, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFnPg <- left_join(test_pred_BHI2_SoFnPg, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoFnPg <- left_join(test_pred_BHI2_SoFnPg, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoFnPg <- test_pred_BHI2_SoFnPg %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoFnPg.csv", test_pred_BHI2_SoFnPg)

### Model for So and Fn in BHI2
# Select the fcs files based on which the model will be trained --> Fn (replicate C), So (replicate A) grown in BHI2
fcs_names_SoFn <- c("20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                    "20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs")

Sample_Info_SoFn <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFn)
Model_RF_SoFn <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFn], Sample_Info_SoFn, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_SoFn <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_B1.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co1_B_1000_SG_C1.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co1_B_1000_SG_D1.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co1_C_1000_SG_E1.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co1_C_1000_SG_F1.fcs",
                         "20200618_Fabian_14strainID_BHI2_48h_Co1_A_1000_SG_A1.fcs",
                         "20200618_Fabian_14strainID_BHI2_48h_Co1_B_1000_SG_B1.fcs",
                         "20200618_Fabian_14strainID_BHI2_48h_Co1_C_1000_SG_C1.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co2_A_1000_SG_A2.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co2_A_1000_SG_B2.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co2_B_1000_SG_C2.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co2_B_1000_SG_D2.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co2_C_1000_SG_G1.fcs",
                         "20200617_Fabian_14strainID_BHI2_24h_Co2_C_1000_SG_H1.fcs",
                         "20200618_Fabian_14strainID_BHI2_48h_Co2_A_1000_SG_A2.fcs",
                         "20200618_Fabian_14strainID_BHI2_48h_Co2_B_1000_SG_B2.fcs",
                         "20200618_Fabian_14strainID_BHI2_48h_Co2_C_1000_SG_D1.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix1_NA_1000_SG_A8.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix1_NA_1000_SG_A9.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix2_NA_1000_SG_A10.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix2_NA_1000_SG_A11.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix3_NA_1000_SG_B8.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix3_NA_1000_SG_B9.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix8_NA_1000_SG_D10.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix8_NA_1000_SG_D11.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix9_NA_1000_SG_E8.fcs",
                         "20200512_Fabian_14strainID_BHI2_NA_Mix9_NA_1000_SG_E9.fcs")

flowData_topre_BHI2_SoFn <- flowData_transformed_gated[fcs_topre_BHI2_SoFn]
test_pred_BHI2_SoFn <- RandomF_predict(x = Model_RF_SoFn[[1]], new_data =  flowData_topre_BHI2_SoFn, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFn <- left_join(test_pred_BHI2_SoFn, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoFn <- left_join(test_pred_BHI2_SoFn, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoFn <- test_pred_BHI2_SoFn %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoFn.csv", test_pred_BHI2_SoFn)


### Model for So, Fn and Pi in BHI2 for comparison in performance compared to FACSVerse
fcs_names_SoFnPi <- c("20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Pi_A_1000_SG_A1.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs")

Sample_Info_sb_SoFnPi <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPi)
Model_RF_SoFnPi <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFnPi], Sample_Info_sb_SoFnPi, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)


### Model for So, Fn, Pg and Vp
# Sample selection Vp
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(85:90, 381:386)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Vp", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> Fn (replicate C), So (replicate A), Pg (replicate B), Vp (replicate B) grown in BHI2
fcs_names_SoFnPgVp <- c("20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                        "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
                        "20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                        "20200512_Fabian_14strainID_BHI2_24h_Vp_B_1000_SG_A10.fcs")

Sample_Info_SoFnPgVp <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPgVp)
Model_RF_SoFnPgVp <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFnPgVp], Sample_Info_SoFnPgVp, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_SoFnPgVp <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_B1.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co1_B_1000_SG_C1.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co1_B_1000_SG_D1.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co1_C_1000_SG_E1.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co1_C_1000_SG_F1.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co1_A_1000_SG_A1.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co1_B_1000_SG_B1.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co1_C_1000_SG_C1.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co2_A_1000_SG_A2.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co2_A_1000_SG_B2.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co2_B_1000_SG_C2.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co2_B_1000_SG_D2.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co2_C_1000_SG_G1.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co2_C_1000_SG_H1.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co2_A_1000_SG_A2.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co2_B_1000_SG_B2.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co2_C_1000_SG_D1.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co3_A_1000_SG_A3.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co3_A_1000_SG_B3.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co3_B_1000_SG_E2.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co3_B_1000_SG_F2.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co3_C_1000_SG_G2.fcs",
                             "20200617_Fabian_14strainID_BHI2_24h_Co3_C_1000_SG_H2.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co3_A_1000_SG_A3.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co3_B_1000_SG_C2.fcs",
                             "20200618_Fabian_14strainID_BHI2_48h_Co3_C_1000_SG_D2.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix1_NA_1000_SG_A8.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix1_NA_1000_SG_A9.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix2_NA_1000_SG_A10.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix2_NA_1000_SG_A11.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix3_NA_1000_SG_B8.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix3_NA_1000_SG_B9.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix8_NA_1000_SG_D10.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix8_NA_1000_SG_D11.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix9_NA_1000_SG_E8.fcs",
                             "20200512_Fabian_14strainID_BHI2_NA_Mix9_NA_1000_SG_E9.fcs")

flowData_topre_BHI2_SoFnPgVp <- flowData_transformed_gated[fcs_topre_BHI2_SoFnPgVp]
test_pred_BHI2_SoFnPgVp <- RandomF_predict(x = Model_RF_SoFnPgVp[[1]], new_data =  flowData_topre_BHI2_SoFnPgVp, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFnPgVp <- left_join(test_pred_BHI2_SoFnPgVp, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoFnPgVp <- left_join(test_pred_BHI2_SoFnPgVp, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoFnPgVp <- test_pred_BHI2_SoFnPgVp %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoFnPgVp.csv", test_pred_BHI2_SoFnPgVp)

### Model for So, Ssal, Ssan, Smi, Sg, Vp
# Sample selection Ssal
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(67:72, 363:368)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Ssal", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection Ssan
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(73:78, 369:374)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Ssan", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection Smi
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(49:54, 345:350)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Smi", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection Sg
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(43:48, 339:344)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Sg", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Ssal (replicate A), Ssan (replicate B), Smi (replicate A), Sg (replicate C), Vp (replicate B) grown in BHI2
fcs_names_SoSsalSsanSmiSgVp <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                                 "20200512_Fabian_14strainID_BHI2_24h_Ssal_A_1000_SG_F2.fcs",
                                 "20200512_Fabian_14strainID_BHI2_24h_Ssan_B_1000_SG_C3.fcs",
                                 "20200512_Fabian_14strainID_BHI2_24h_Smi_A_1000_SG_F7.fcs",
                                 "20200512_Fabian_14strainID_BHI2_24h_Sg_C_1000_SG_G5.fcs",
                                 "20200512_Fabian_14strainID_BHI2_24h_Vp_B_1000_SG_A10.fcs")

Sample_Info_SoSsalSsanSmiSgVp <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoSsalSsanSmiSgVp)
Model_RF_SoSsalSsanSmiSgVp <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoSsalSsanSmiSgVp], Sample_Info_SoSsalSsanSmiSgVp, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_SoSsalSsanSmiSgVp <- c("20200617_Fabian_14strainID_BHI2_24h_Co4_A_1000_SG_C3.fcs",
                                      "20200617_Fabian_14strainID_BHI2_24h_Co4_A_1000_SG_D3.fcs",
                                      "20200617_Fabian_14strainID_BHI2_24h_Co4_B_1000_SG_E3.fcs",
                                      "20200617_Fabian_14strainID_BHI2_24h_Co4_B_1000_SG_F3.fcs",
                                      "20200617_Fabian_14strainID_BHI2_24h_Co4_C_1000_SG_G3.fcs",
                                      "20200617_Fabian_14strainID_BHI2_24h_Co4_C_1000_SG_H3.fcs",
                                      "20200618_Fabian_14strainID_BHI2_48h_Co4_A_1000_SG_B3.fcs",
                                      "20200618_Fabian_14strainID_BHI2_48h_Co4_B_1000_SG_C3.fcs",
                                      "20200618_Fabian_14strainID_BHI2_48h_Co4_C_1000_SG_D3.fcs",
                                      "20200512_Fabian_14strainID_BHI2_NA_Mix4_NA_1000_SG_B10.fcs",
                                      "20200512_Fabian_14strainID_BHI2_NA_Mix4_NA_1000_SG_B11.fcs")

flowData_topre_BHI2_SoSsalSsanSmiSgVp <- flowData_transformed_gated[fcs_topre_BHI2_SoSsalSsanSmiSgVp]
test_pred_BHI2_SoSsalSsanSmiSgVp <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgVp[[1]], new_data =  flowData_topre_BHI2_SoSsalSsanSmiSgVp, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgVp <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVp, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoSsalSsanSmiSgVp <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVp, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoSsalSsanSmiSgVp <- test_pred_BHI2_SoSsalSsanSmiSgVp %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoSsalSsanSmiSgVp.csv", test_pred_BHI2_SoSsalSsanSmiSgVp)


### Model for So, Ssal, Ssan, Smi, Sg, Vp, Av, An
# Sample selection Av
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(13:18, 309:314)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Av", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection An
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(7:12, 303:308)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check An", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Ssal (replicate A), Ssan (replicate B), Smi (replicate A), Sg (replicate C), Vp (replicate B), Av (replicate A), An (replicate A, B would be better but mix contains cells from replicate A) grown in BHI2
fcs_names_SoSsalSsanSmiSgVpAvAn <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Ssal_A_1000_SG_F2.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Ssan_B_1000_SG_C3.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Smi_A_1000_SG_F7.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Sg_C_1000_SG_G5.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Vp_B_1000_SG_A10.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Av_A_1000_SG_E1.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_An_A_1000_SG_G7.fcs")

Sample_Info_SoSsalSsanSmiSgVpAvAn <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoSsalSsanSmiSgVpAvAn)
Model_RF_SoSsalSsanSmiSgVpAvAn <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoSsalSsanSmiSgVpAvAn], Sample_Info_SoSsalSsanSmiSgVpAvAn, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_SoSsalSsanSmiSgVpAvAn <- c("20200617_Fabian_14strainID_BHI2_24h_Co4_A_1000_SG_C3.fcs",
                                          "20200617_Fabian_14strainID_BHI2_24h_Co4_A_1000_SG_D3.fcs",
                                          "20200617_Fabian_14strainID_BHI2_24h_Co4_B_1000_SG_E3.fcs",
                                          "20200617_Fabian_14strainID_BHI2_24h_Co4_B_1000_SG_F3.fcs",
                                          "20200617_Fabian_14strainID_BHI2_24h_Co4_C_1000_SG_G3.fcs",
                                          "20200617_Fabian_14strainID_BHI2_24h_Co4_C_1000_SG_H3.fcs",
                                          "20200618_Fabian_14strainID_BHI2_48h_Co4_A_1000_SG_B3.fcs",
                                          "20200618_Fabian_14strainID_BHI2_48h_Co4_B_1000_SG_C3.fcs",
                                          "20200618_Fabian_14strainID_BHI2_48h_Co4_C_1000_SG_D3.fcs",
                                          "20200512_Fabian_14strainID_BHI2_NA_Mix4_NA_1000_SG_B10.fcs",
                                          "20200512_Fabian_14strainID_BHI2_NA_Mix4_NA_1000_SG_B11.fcs")

flowData_topre_BHI2_SoSsalSsanSmiSgVpAvAn <- flowData_transformed_gated[fcs_topre_BHI2_SoSsalSsanSmiSgVpAvAn]
test_pred_BHI2_SoSsalSsanSmiSgVpAvAn <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgVpAvAn[[1]], new_data = flowData_topre_BHI2_SoSsalSsanSmiSgVpAvAn, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgVpAvAn <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVpAvAn, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoSsalSsanSmiSgVpAvAn <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVpAvAn, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoSsalSsanSmiSgVpAvAn <- test_pred_BHI2_SoSsalSsanSmiSgVpAvAn %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoSsalSsanSmiSgVpAvAn.csv", test_pred_BHI2_SoSsalSsanSmiSgVpAvAn)














# Run Random Forest classifier to predict the Strain based on the single-cell FCM data
# fcs files for co-culture 1 and 2 are selected
fcs_topre <- c("20200617_Fabian_14strainID_MM_24h_Co2_A_100_SG_A7.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co2_B_100_SG_C7.fcs",
               "20200617_Fabian_14strainID_BHI2_24h_Co2_A_1000_SG_A2.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_A_100_SG_A8.fcs",
               "20200617_Fabian_14strainID_BHI2_24h_Co2_B_1000_SG_C2.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_B_100_SG_C7.fcs",
               "20200617_Fabian_14strainID_BHI2_24h_Co2_C_1000_SG_G1.fcs",
               "20200618_Fabian_14strainID_BHI2_48h_Co2_A_1000_SG_A2.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_B_100_SG_C7.fcs",
               "20200618_Fabian_14strainID_BHI2_48h_Co2_B_1000_SG_B2.fcs",
               "20200618_Fabian_14strainID_BHI2_48h_Co2_C_1000_SG_D1.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co1_B_100_SG_C6.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co1_C_100_SG_E6.fcs",
               "20200618_Fabian_14strainID_BHI2_48h_Co1_A_1000_SG_A1.fcs",
               "20200618_Fabian_14strainID_BHI2_48h_Co1_B_1000_SG_B1.fcs",
               "20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
               "20200617_Fabian_14strainID_BHI2_24h_Co1_B_1000_SG_C1.fcs",
               "20200617_Fabian_14strainID_BHI2_24h_Co1_C_1000_SG_E1.fcs",
               "20200618_Fabian_14strainID_BHI2_48h_Co1_C_1000_SG_C1.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co1_B_100_SG_A7.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co1_C_100_SG_B7.fcs")

fcs_topre_MM <- c("20200617_Fabian_14strainID_MM_24h_Co2_A_100_SG_A7.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co2_B_100_SG_C7.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_A_100_SG_A8.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_B_100_SG_C7.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_B_100_SG_C7.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co1_B_100_SG_C6.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co1_C_100_SG_E6.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co1_B_100_SG_A7.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co1_C_100_SG_B7.fcs")

flowData_topre_MM <- flowData_transformed_2[fcs_topre_MM]

test_pred_MM <- RandomF_predict(x = Model_RF[[1]], new_data =  flowData_topre_MM, cleanFCS = FALSE)

# fcs files for artificial mixes are selected
fcs_topre_BHI2_AM <- c("20200512_Fabian_14strainID_BHI2_NA_Mix1_NA_1000_SG_A8.fcs",
                       "20200512_Fabian_14strainID_BHI2_NA_Mix1_NA_1000_SG_A9.fcs",
                       "20200512_Fabian_14strainID_BHI2_NA_Mix2_NA_1000_SG_A10.fcs",
                       "20200512_Fabian_14strainID_BHI2_NA_Mix2_NA_1000_SG_A11.fcs",
                       "20200512_Fabian_14strainID_BHI2_NA_Mix8_NA_1000_SG_D10.fcs",
                       "20200512_Fabian_14strainID_BHI2_NA_Mix8_NA_1000_SG_D11.fcs",
                       "20200512_Fabian_14strainID_BHI2_NA_Mix9_NA_1000_SG_E8.fcs",
                       "20200512_Fabian_14strainID_BHI2_NA_Mix9_NA_1000_SG_E9.fcs")

flowData_topre_BHI2_AM <- flowData_transformed_gated[fcs_topre_BHI2_AM]

test_pred_BHI2_AM <- RandomF_predict(x = Model_RF_SoFnPg[[1]], new_data =  flowData_topre_BHI2_AM, cleanFCS = FALSE)
test_pred_BHI2_AM

## Model based on two strains (Fn, So)
# Select the fcs files based on which the model will be trained --> Fn (replicate C), So (replicate A), Pg (replicate B)
fcs_names_sofn <- c("20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
               "20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs")

# Sample info has to contain a column called 'name' which matches the sammplenames of the fcs files
Sample_Info_sb_sofn <- Sample_Info %>% dplyr::filter(name %in% fcs_names_sofn)

Model_RF_sofn <- Phenoflow::RandomF_FCS(flowData_transformed_2[fcs_names_sofn], Sample_Info_sb_sofn, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

# Predict artificial mixes with 2 strain model
test_pred_BHI2_AM_sofn <- RandomF_predict(x = Model_RF_sofn[[1]], new_data =  flowData_topre_BHI2_AM, cleanFCS = FALSE)
test_pred_BHI2_AM_sofn
