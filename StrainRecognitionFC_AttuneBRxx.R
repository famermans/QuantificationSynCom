#### Clearing workspace, loading libraries, setting seed ----

## Clear environment and set working directory
rm(list = ls())
setwd("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM")

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
library("ggcyto")

set.seed(777)


#### Loading data ----

Datapath <- "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/data"
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

metadata$Well <- substr(metadata$Well, start = 1, stop = nchar(metadata$Well)-4)

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
p_singlets <- xyplot(`BL1-H`~`BL1-W`, data = flowData_transformed[c(1:729)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0,16))),
       axix = axis.default, nbin = 125, main = "QC singlets (BL1-A ~ BL1-W)", xlab = "BL1-W", ylab = "BL1-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)


#### Select data ----
flowData_transformed_sel <- flowData_transformed[c(1:122, 297:413, 656:696)]
metadata_sel <- metadata[c(1:122, 297:413, 656:696), ]

flowData_transformed_BHI <- flowData_transformed[c(1:122, 656:675)]
metadata_BHI <- metadata[c(1:122, 656:675), ]

rownames(metadata_BHI) <- NULL

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

### Singlet analysis pure cultures
p_singlets_pure_FSC <- xyplot(`FSC-W`~`FSC-H`, data = flowData_transformed_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                              scales = list(y = list(limits = c(0, 10)),
                                            x = list(limits = c(5, 16))),
                              axis = axis.default, nbin = 125, main = "Singlet analysis pure cultures FSC", xlab = "FSC-H", ylab = "FSC-W",
                              par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_singlets_pure_SSC <- xyplot(`SSC-W`~`SSC-H`, data = flowData_transformed_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                              scales = list(y = list(limits = c(0, 10)),
                                            x = list(limits = c(5, 16))),
                              axis = axis.default, nbin = 125, main = "Singlet analysis pure cultures SSC", xlab = "SSC-H", ylab = "SSC-W",
                              par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

print(p_singlets_pure_FSC)
print(p_singlets_pure_SSC)


p_singlets_pure_FSC_A <- xyplot(`FSC-A`~`FSC-H`, data = flowData_transformed_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                                scales = list(y = list(limits = c(0, 16)),
                                              x = list(limits = c(5, 16))),
                                axis = axis.default, nbin = 125, main = "Singlet analysis pure cultures FSC", xlab = "FSC-H", ylab = "FSC-A",
                                par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_singlets_pure_SSC_A <- xyplot(`SSC-A`~`SSC-H`, data = flowData_transformed_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                                scales = list(y = list(limits = c(0, 16)),
                                              x = list(limits = c(5, 16))),
                                axis = axis.default, nbin = 125, main = "Singlet analysis pure cultures SSC", xlab = "SSC-H", ylab = "SSC-A",
                                par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

print(p_singlets_pure_FSC_A)
print(p_singlets_pure_SSC_A)

xyplot(`SSC-A`~`BL1-A`, data = flowData_transformed_gated[c(19, 21)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(5, 16))),
       axis = axis.default, nbin = 125, main = "Check blancs", xlab = "BL1-A", ylab = "SSC-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_gated[c(19, 21)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(5, 16))),
       axis = axis.default, nbin = 125, main = "Check blancs", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(19, 21)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(5, 16))),
       axis = axis.default, nbin = 125, main = "Check blancs", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Use none gated data
p_singlets_pure_FSC_NG <- xyplot(`FSC-W`~`FSC-H`, data = flowData_transformed[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                                 scales = list(y = list(limits = c(0, 10)),
                                               x = list(limits = c(5, 16))),
                                 axis = axis.default, nbin = 125, main = "Singlet analysis pure cultures FSC", xlab = "FSC-H", ylab = "FSC-W",
                                 par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_singlets_pure_SSC_NG <- xyplot(`SSC-W`~`SSC-H`, data = flowData_transformed[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                                 scales = list(y = list(limits = c(0, 10)),
                                               x = list(limits = c(5, 16))),
                                 axis = axis.default, nbin = 125, main = "Singlet analysis pure cultures SSC", xlab = "SSC-H", ylab = "SSC-W",
                                 par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

print(p_singlets_pure_FSC_NG)
print(p_singlets_pure_SSC_NG)


p_singlets_pure_FSC_A_NG <- xyplot(`FSC-A`~`FSC-H`, data = flowData_transformed[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                                   scales = list(y = list(limits = c(0, 16)),
                                                 x = list(limits = c(5, 16))),
                                   axis = axis.default, nbin = 125, main = "Singlet analysis pure cultures FSC", xlab = "FSC-H", ylab = "FSC-A",
                                   par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_singlets_pure_SSC_A_NG <- xyplot(`SSC-A`~`SSC-H`, data = flowData_transformed[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                                   scales = list(y = list(limits = c(0, 16)),
                                                 x = list(limits = c(5, 16))),
                                   axis = axis.default, nbin = 125, main = "Singlet analysis pure cultures SSC", xlab = "SSC-H", ylab = "SSC-A",
                                   par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

print(p_singlets_pure_FSC_A_NG)
print(p_singlets_pure_SSC_A_NG)

## Density plots

p_density_SSC_W <- autoplot(flowData_transformed_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], "SSC-W")+
        labs(title = "Density plot SSC-W",
             x = "SSC-W",
             y = "Density")
print(p_density_SSC_W)

# p_density_SSC_A <- autoplot(flowData_transformed_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], "SSC-A")+
#         labs(title = "Density plot SSC-A",
#              x = "SSC-A",
#              y = "Density")
# print(p_density_SSC_A)
# 
# p_density_SSC_H <- autoplot(flowData_transformed_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], "SSC-H")+
#         labs(title = "Density plot SSC-H",
#              x = "SSC-H",
#              y = "Density")
# print(p_density_SSC_H)


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
# Sample selection So
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(61:66, 357:362)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check So", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection Fn
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(23:28, 319:324)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Fn", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection Pg
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(31:36, 327:332)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Pg", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Fn (replicate C), Pg (replicate B) grown in BHI2
fcs_names_SoFnPg <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs")

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
# Select the fcs files based on which the model will be trained --> So (replicate A), Fn (replicate C) grown in BHI2
fcs_names_SoFn <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                    "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs")

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
# Sample selection Pi
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(37:42, 333:338)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Pi", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Fn (replicate C), Pi (replicate A) grown in BHI2
fcs_names_SoFnPi <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Pi_A_1000_SG_A1.fcs")

Sample_Info_SoFnPi <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPi)
Model_RF_SoFnPi <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFnPi], Sample_Info_SoFnPi, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)


### Model for So, Fn, Pg and Vp
# Sample selection Vp
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(85:90, 381:386)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Vp", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Fn (replicate C), Pg (replicate B), Vp (replicate B) grown in BHI2
fcs_names_SoFnPgVp <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                        "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                        "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
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
# REMARK: not enough cells for An to train model --> left out of the calculations

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


### Model for So, Ssal, Ssan, Smi, Sg, Smu, Ssob
# Sample selection Smu
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(55:60, 351:356)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Smu", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection Ssob
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(79:84, 375:380)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Ssob", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Ssal (replicate A), Ssan (replicate B), Smi (replicate A), Sg (replicate C), Smu (replicate A), Ssob (replicate A) grown in BHI2
fcs_names_SoSsalSsanSmiSgSmuSsob <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Ssal_A_1000_SG_F2.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Ssan_B_1000_SG_C3.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Smi_A_1000_SG_F7.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Sg_C_1000_SG_G5.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Smu_A_1000_SG_B7.fcs",
                                     "20200512_Fabian_14strainID_BHI2_24h_Ssob_A_1000_SG_C7.fcs")

Sample_Info_SoSsalSsanSmiSgSmuSsob <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoSsalSsanSmiSgSmuSsob)
Model_RF_SoSsalSsanSmiSgSmuSsob <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoSsalSsanSmiSgSmuSsob], Sample_Info_SoSsalSsanSmiSgSmuSsob, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_SoSsalSsanSmiSgSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix5_NA_1000_SG_C8.fcs",
                                           "20200512_Fabian_14strainID_BHI2_NA_Mix5_NA_1000_SG_C9.fcs")

flowData_topre_BHI2_SoSsalSsanSmiSgSmuSsob <- flowData_transformed_gated[fcs_topre_BHI2_SoSsalSsanSmiSgSmuSsob]
test_pred_BHI2_SoSsalSsanSmiSgSmuSsob <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgSmuSsob[[1]], new_data = flowData_topre_BHI2_SoSsalSsanSmiSgSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgSmuSsob <- left_join(test_pred_BHI2_SoSsalSsanSmiSgSmuSsob, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoSsalSsanSmiSgSmuSsob <- left_join(test_pred_BHI2_SoSsalSsanSmiSgSmuSsob, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoSsalSsanSmiSgSmuSsob <- test_pred_BHI2_SoSsalSsanSmiSgSmuSsob %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoSsalSsanSmiSgSmuSsob.csv", test_pred_BHI2_SoSsalSsanSmiSgSmuSsob)


### Model for Aa, Fn, Pg, Pi, Smu, Ssob
# Sample selection Aa
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1:6, 297:302)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Aa", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> Aa (replicate A), Fn (replicate C), Pg (replicate B), Pi (replicate A), Smu (replicate A), Ssob (replicate A) grown in BHI2
fcs_names_AaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_24h_Aa_A_1000_SG_D7.fcs",
                               "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                               "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
                               "20200512_Fabian_14strainID_BHI2_24h_Pi_A_1000_SG_A1.fcs",
                               "20200512_Fabian_14strainID_BHI2_24h_Smu_A_1000_SG_B7.fcs",
                               "20200512_Fabian_14strainID_BHI2_24h_Ssob_A_1000_SG_C7.fcs")

Sample_Info_AaFnPgPiSmuSsob <- Sample_Info %>% dplyr::filter(name %in% fcs_names_AaFnPgPiSmuSsob)
Model_RF_AaFnPgPiSmuSsob <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_AaFnPgPiSmuSsob], Sample_Info_AaFnPgPiSmuSsob, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_AaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix6_NA_1000_SG_C10.fcs",
                                    "20200512_Fabian_14strainID_BHI2_NA_Mix6_NA_1000_SG_C11.fcs",
                                    "20200617_Fabian_14strainID_BHI2_24h_Co5_A_1000_SG_A4.fcs",
                                    "20200617_Fabian_14strainID_BHI2_24h_Co5_A_1000_SG_B4.fcs",
                                    "20200617_Fabian_14strainID_BHI2_24h_Co5_B_1000_SG_C4.fcs",
                                    "20200617_Fabian_14strainID_BHI2_24h_Co5_B_1000_SG_D4.fcs",
                                    "20200617_Fabian_14strainID_BHI2_24h_Co5_C_1000_SG_E4.fcs",
                                    "20200617_Fabian_14strainID_BHI2_24h_Co5_C_1000_SG_F4.fcs",
                                    "20200618_Fabian_14strainID_BHI2_48h_Co5_A_1000_SG_A4.fcs",
                                    "20200618_Fabian_14strainID_BHI2_48h_Co5_B_1000_SG_B4.fcs",
                                    "20200618_Fabian_14strainID_BHI2_48h_Co5_C_1000_SG_C4.fcs")

flowData_topre_BHI2_AaFnPgPiSmuSsob <- flowData_transformed_gated[fcs_topre_BHI2_AaFnPgPiSmuSsob]
test_pred_BHI2_AaFnPgPiSmuSsob <- RandomF_predict(x = Model_RF_AaFnPgPiSmuSsob[[1]], new_data = flowData_topre_BHI2_AaFnPgPiSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_AaFnPgPiSmuSsob <- left_join(test_pred_BHI2_AaFnPgPiSmuSsob, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_AaFnPgPiSmuSsob <- left_join(test_pred_BHI2_AaFnPgPiSmuSsob, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_AaFnPgPiSmuSsob <- test_pred_BHI2_AaFnPgPiSmuSsob %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsAaFnPgPiSmuSsob.csv", test_pred_BHI2_AaFnPgPiSmuSsob)


### Model for So, Ssal, Ssan, Smi, Sg, Vp, Aa, Fn, Pg, Pi, Smu, Ssob
# Select the fcs files based on which the model will be trained --> So (replicate A), Ssal (replicate A), Ssan (replicate B), Smi (replicate A), Sg (replicate C), Vp (replicate B), Aa (replicate A), Fn (replicate C), Pg (replicate B), Pi (replicate A), Smu (replicate A), Ssob (replicate A) grown in BHI2
fcs_names_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssal_A_1000_SG_F2.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssan_B_1000_SG_C3.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Smi_A_1000_SG_F7.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Sg_C_1000_SG_G5.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Vp_B_1000_SG_A10.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Aa_A_1000_SG_D7.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Pi_A_1000_SG_A1.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Smu_A_1000_SG_B7.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssob_A_1000_SG_C7.fcs")

Sample_Info_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob)
Model_RF_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob], Sample_Info_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix7_NA_1000_SG_D8.fcs",
                                                     "20200512_Fabian_14strainID_BHI2_NA_Mix7_NA_1000_SG_D9.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_A_1000_SG_A5.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_A_1000_SG_B5.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_B_1000_SG_C5.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_B_1000_SG_D5.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_C_1000_SG_G4.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_C_1000_SG_H4.fcs",
                                                     "20200618_Fabian_14strainID_BHI2_48h_Co6_A_1000_SG_A5.fcs",
                                                     "20200618_Fabian_14strainID_BHI2_48h_Co6_B_1000_SG_B5.fcs",
                                                     "20200618_Fabian_14strainID_BHI2_48h_Co6_C_1000_SG_D4.fcs")

flowData_topre_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- flowData_transformed_gated[fcs_topre_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob]
test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob[[1]], new_data = flowData_topre_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoSsalSsanSmiSgVpAaFnPgPiSmuSsob.csv", test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob)


### Model for So, Ssal, Ssan, Smi, Sg, Vp, Av, An, Aa, Fn, Pg, Pi, Smu, Ssob
# Select the fcs files based on which the model will be trained --> So (replicate A), Ssal (replicate A), Ssan (replicate B), Smi (replicate A), Sg (replicate C), Vp (replicate B), Av (replicate A), An (replicate A, B would be better but mixes were prepared with replicate A), Aa (replicate A), Fn (replicate C), Pg (replicate B), Pi (replicate A), Smu (replicate A), Ssob (replicate A) grown in BHI2
fcs_names_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssal_A_1000_SG_F2.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssan_B_1000_SG_C3.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Smi_A_1000_SG_F7.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Sg_C_1000_SG_G5.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Vp_B_1000_SG_A10.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Av_A_1000_SG_E1.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_An_A_1000_SG_G7.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Aa_A_1000_SG_D7.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Pi_A_1000_SG_A1.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Smu_A_1000_SG_B7.fcs",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssob_A_1000_SG_C7.fcs")

Sample_Info_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob)
Model_RF_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob], Sample_Info_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
# REMARK: not enough cells for An to train model --> left out of the calculations

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix7_NA_1000_SG_D8.fcs",
                                                     "20200512_Fabian_14strainID_BHI2_NA_Mix7_NA_1000_SG_D9.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_A_1000_SG_A5.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_A_1000_SG_B5.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_B_1000_SG_C5.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_B_1000_SG_D5.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_C_1000_SG_G4.fcs",
                                                     "20200617_Fabian_14strainID_BHI2_24h_Co6_C_1000_SG_H4.fcs",
                                                     "20200618_Fabian_14strainID_BHI2_48h_Co6_A_1000_SG_A5.fcs",
                                                     "20200618_Fabian_14strainID_BHI2_48h_Co6_B_1000_SG_B5.fcs",
                                                     "20200618_Fabian_14strainID_BHI2_48h_Co6_C_1000_SG_D4.fcs")

flowData_topre_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- flowData_transformed_gated[fcs_topre_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob]
test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob[[1]], new_data = flowData_topre_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob.csv", test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob)


#### Some tests with training with replicates ----
### Model for So, Fn and Pg in BHI2
# Select the fcs files based on which the model will be trained --> So, Fn, Pg grown in BHI2
fcs_names_SoFnPg2 <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D2.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_So_B_1000_SG_D3.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_So_B_1000_SG_D4.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_So_C_1000_SG_D5.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_So_C_1000_SG_D6.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Fn_A_1000_SG_B1.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Fn_A_1000_SG_B2.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Fn_B_1000_SG_B3.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Fn_B_1000_SG_B4.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B6.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Pg_A_1000_SG_E7.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Pg_A_1000_SG_E8.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E9.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Pg_C_1000_SG_E11.fcs",
                       "20200512_Fabian_14strainID_BHI2_24h_Pg_C_1000_SG_E12.fcs")

Sample_Info_SoFnPg2 <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPg2)

# If samples have a lower count than the 'downsample' argument, these samples will be left out and you get over/under representation of certain strains
Model_RF_SoFnPg2 <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFnPg2], Sample_Info_SoFnPg2, target_label = "Strain", downsample = 7957, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)


## Alternative approach for dealing with this over/under representation: using FSC_pool
# Select the fcs files based on which the model will be trained --> So, Fn, Pg grown in BHI2
flowData_transformed_test <- flowData_transformed_gated

flowData_pooled <- FCS_pool(x = flowData_transformed_test,
                            stub = c("20200512_Fabian_14strainID_BHI2_24h_So*",
                                     "20200512_Fabian_14strainID_BHI2_24h_Fn*",
                                     "20200512_Fabian_14strainID_BHI2_24h_Pg*"))

metadata_pooled <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_pooled), "_"), rbind)))
colnames(metadata_pooled) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain")

name_pooled <- flowCore::sampleNames(flowData_pooled)
Sample_Info_pooled <- cbind(name_pooled, metadata_pooled)

# Plot not working for some reason
xyplot(`BL3-A`~`BL1-A`, data = flowData_pooled[c(1:3)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Test pooled samples", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

cells2 <- flowCore::filter(flowData_pooled, polyGate5)
TotalCount2 <- summary(cells2)
TotalCount2 <- toTable(TotalCount2)


fcs_names_SoFnPg3 <- c("20200512_Fabian_14strainID_BHI2_24h_So*",
                       "20200512_Fabian_14strainID_BHI2_24h_Fn*",
                       "20200512_Fabian_14strainID_BHI2_24h_Pg*")

Sample_Info_SoFnPg3 <- Sample_Info_pooled %>% dplyr::filter(name_pooled %in% fcs_names_SoFnPg3)

Model_RF_SoFnPg3 <- Phenoflow::RandomF_FCS(flowData_pooled[fcs_names_SoFnPg3], Sample_Info_SoFnPg3, sample_col = "name_pooled", target_label = "Strain", downsample = 50000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)










########################################################
# fcs files for co-culture 1 and 2 in MM are selected
fcs_topre_MM <- c("20200617_Fabian_14strainID_MM_24h_Co2_A_100_SG_A7.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co2_B_100_SG_C7.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_A_100_SG_A8.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_B_100_SG_C7.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co2_B_100_SG_C7.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co1_B_100_SG_C6.fcs",
               "20200617_Fabian_14strainID_MM_24h_Co1_C_100_SG_E6.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co1_B_100_SG_A7.fcs",
               "20200618_Fabian_14strainID_MM_48h_Co1_C_100_SG_B7.fcs")

