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
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed,
       scales = list(y = list(limits = c(0, 20)),
                     x = list(limits = c(0, 20))),
       axis = axis.default, nbin = 125, main = "QC (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)