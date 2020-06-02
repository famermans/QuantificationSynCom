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
                                  `BL3-A` = asinh(`BL3-A`))
param = c("FSC-A", "SSC-A", "BL1-A", "BL3-A")

# remove(flowData)

#### Extrating metadata from .fcs files ----
sampleNames(flowData_transformed) <- substring(sampleNames(flowData_transformed), 0, nchar(sampleNames(flowData_transformed))-4)       # nchar takes a character vector and returns the number of characters in the vector

# Extracting the actual metadata
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed), "_"), rbind)))
colnames(metadata) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain", "Replicate", "Dilution", "Stain", "Well")

#### Quality control (contamination) ----
## BHI2
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1:296)],
       scales = list(y = list(limits = c(0, 20)),
                     x = list(limits = c(0, 20))),
       axis = axis.default, nbin = 125, main = "QC BHI2 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

## MM
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(297:503)],
       scales = list(y = list(limits = c(0, 20)),
                     x = list(limits = c(0, 20))),
       axis = axis.default, nbin = 125, main = "QC MM (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

#### Gating ----

### BHI2
## Gating out extra background based on SSC-BL1
# Constructing gate
sqrcut1 <- matrix(c(5, 10, 14, 10, 5,
                    2, 5, 13.5, 13.5, 8.5), ncol = 2, nrow = 5)
colnames(sqrcut1) <- c("BL1-A", "SSC-A")
polyGate1 <- polygonGate(.gate = sqrcut1, filterId = "Reduced background BL1-SSC BHI2")

# Gating quality check
xyplot(`SSC-A`~`BL1-A`, data = flowData_transformed[c(19:22, 29, 30, 37, 55, 85)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating BL1-SSC BHI2", xlab = "BL1-A", ylab = "SSC-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

## Gating cells based on BL1-BL3
sqrcut2 <- matrix(c(5, 10, 14, 10, 5,
                    2, 3, 4, 5, 6), ncol = 2, nrow = 5)
colnames(sqrcut2) <- c("BL1-A", "BL3-A")
polyGate2 <- polygonGate(.gate = sqrcut2, filterId = "Cells BHI2")

# Gating quality check
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(21:24)], filter = polyGate2,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

### MM
## Gating out extra background based on SSC-BL1
# Constructing gate
sqrcut3 <- matrix(c(5, 10, 14, 10, 5,
                    2, 5, 13.5, 13.5, 8.5), ncol = 2, nrow = 5)
colnames(sqrcut3) <- c("BL1-A", "SSC-A")
polyGate3 <- polygonGate(.gate = sqrcut3, filterId = "Reduced background BL1-SSC MM")

# Gating quality check
xyplot(`SSC-A`~`BL1-A`, data = flowData_transformed[c(21:24)], filter = polyGate3,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating BL1-SSC MM", xlab = "BL1-A", ylab = "SSC-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

## Gating cells based on BL1-BL3
sqrcut4 <- matrix(c(5, 10, 14, 10, 5,
                    2, 3, 4, 5, 6), ncol = 2, nrow = 5)
colnames(sqrcut4) <- c("BL1-A", "BL3-A")
polyGate4 <- polygonGate(.gate = sqrcut4, filterId = "Cells MM")

# Gating quality check
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(21:24)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)
