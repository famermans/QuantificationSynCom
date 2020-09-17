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
                                  `BL3-H` = asinh(`BL3-H`))
param = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "FSC-H", "SSC-H", "BL1-H", "BL3-H")


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
TotalCount <- summary(cells); TotalCount <- toTable(TotalCount)

### Extracting volumes
# Volumes are in µL
vol <- as.numeric(flowCore::fsApply(flowData_transformed, FUN = function(x) x@description$`$VOL`))/1000

### Concentrations
# Calculating concentrations --> Concentrations will be in cells/µL
cell_concentrations <- data.frame(Samples = flowCore::sampleNames(flowData_transformed),
                                  Strain = metadata$Strain,
                                  Concentration = (TotalCount$true*metadata$Dilution)/vol)


#### Phenotypic diversity analysis ----

### Normalization of data
summary <- fsApply(x = flowData_transformed_gated, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max = max(summary[, "BL1-A"])
max2 = max(summary[, "BL1-H"])
mytrans <- function(x) x/max
mytrans2 <- function(x) x/max2
flowData_transformed_2 <- transform(flowData_transformed_gated,
                                    `FSC-A` = mytrans(`FSC-A`), 
                                    `SSC-A` = mytrans(`SSC-A`), 
                                    `BL1-A` = mytrans(`BL1-A`), 
                                    `BL3-A` = mytrans(`BL3-A`),
                                    `FSC-H` = mytrans2(`FSC-H`), 
                                    `SSC-H` = mytrans2(`SSC-H`), 
                                    `BL1-H` = mytrans2(`BL1-H`), 
                                    `BL3-H` = mytrans2(`BL3-H`))

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

# Sample selection
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

# Define parameters on which RF will build its model
paramRF = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "FSC-H", "SSC-H", "BL1-H", "BL3-H")

## Model for So, Fn and Pg in BHI2
# Select the fcs files based on which the model will be trained --> Fn (replicate C), So (replicate A), Pg (replicate B)
fcs_names <- c("20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
               "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
               "20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs")

# Sample info has to contain a column called 'name' which matches the sammplenames of the fcs files
Sample_Info_sb <- Sample_Info %>% dplyr::filter(name %in% fcs_names)

#Model_RF <- Phenoflow::RandomF_FCS(flowData_transformed_2[fcs_names], Sample_Info_sb, target_label = "Strain", downsample = 1000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

# Using Jasmine's function
Model_RF <- RandomF_FCS_tmp(flowData_transformed_2[fcs_names], Sample_Info_sb, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
Model_RF_NotNormalized <- RandomF_FCS_tmp(flowData_transformed[fcs_names], Sample_Info_sb, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)


## Model for So, Fn and Pi in BHI2
fcs_names_pi <- c("20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                  "20200512_Fabian_14strainID_BHI2_24h_Pi_A_1000_SG_A1.fcs",
                  "20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs")

Sample_Info_sb_pi <- Sample_Info %>% dplyr::filter(name %in% fcs_names_pi)

Model_RF_NotNormalized_Pi <- RandomF_FCS_tmp(flowData_transformed[fcs_names_pi], Sample_Info_sb_pi, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
Model_RF_Pi <- RandomF_FCS_tmp(flowData_transformed_2[fcs_names_pi], Sample_Info_sb_pi, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

# Run Random Forest classifier to predict the Strain based on the single-cell FCM data
# Choose the fcs files in which the model will be used/tested
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

flowData_topre <- flowData_transformed_2[fcs_topre]


test_pred <- RandomF_predict(x = Model_RF[[1]], new_data =  flowData_topre, cleanFCS = FALSE)
test_pred
