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
                                  `BL1-H` = asinh(`BL1-H`),
                                  `BL3-H` = asinh(`BL3-H`))
param = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "BL1-H", "BL3-H")


#### Extrating metadata from .fcs files ----
sampleNames(flowData_transformed) <- substring(sampleNames(flowData_transformed), 0, nchar(sampleNames(flowData_transformed))-4)       # nchar takes a character vector and returns the number of characters in the vector

# Extracting the actual metadata
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed), "_"), rbind)))
colnames(metadata) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain", "Replicate", "Dilution", "Stain", "Well")

metadata$Dilution <- as.numeric(metadata$Dilution)

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
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], filter = polyGate2,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check mixes BHI2
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119)], filter = polyGate2,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2 mixes", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check co-cultures BHI2
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(656, 659, 662, 665, 668, 671, 674, 675)], filter = polyGate2,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2 co-cultures", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check MM
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(297, 303, 309, 315, 317, 319, 325, 327, 333, 339, 345, 351, 357, 363, 369, 375, 381, 413)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check mixes MM
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(387, 389, 391, 393, 395, 397, 399, 401, 403, 405, 407, 409, 411)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM mixes", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check co-cultures MM
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(678, 680, 683, 686, 689, 692, 695, 696)], filter = polyGate4,
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
mytrans <- function(x) x/max
flowData_transformed_2 <- transform(flowData_transformed_gated,
                                    `FSC-A` = mytrans(`FSC-A`), 
                                    `SSC-A` = mytrans(`SSC-A`), 
                                    `BL1-A` = mytrans(`BL1-A`), 
                                    `BL3-A` = mytrans(`BL3-A`))

### Calculating fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed_2, param, nbin = 128, 
                    bw = 0.01, normalize = function(x) x)

### Calculate ecological parameters
Diversity.fbasis <- Diversity(fbasis, d = 3, plot = FALSE, R = 999)
#Evenness.fbasis <- Evenness(fbasis, d = 3, plot = FALSE)
#Structural.organization.fbasis <- So(fbasis, d = 3, plot = FALSE)
#Coef.var.fbasis <- CV(fbasis, d = 3, plot = FALSE)

########################################## Check script from here

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

# Plot ordination
plot_beta_fcm(beta.div, color = metadata$State, shape = as.factor(metadata$Time), labels = list("State", "Timepoint")) + 
        theme_bw() +
        geom_point(size = 8, alpha = 0.5)

#Prediction of relative abundances in the mixed cultures

#### Random Forest
# Build random forest for Fn, Pi, and So

# Select the fcs files based on which the model will be trained
fcs_names <- c("Fn_10_3_20191015_161647.fcs", "Pi_10_3_20191015_161828.fcs", "So_10_3_20191015_160057.fcs")

# Sample info has to contain a column called 'name' which matches the sammplenames of the fcs files
name <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = FALSE)
Sample_Info <- cbind(name, metadata)
Sample_Info_sb <- Sample_Info %>% dplyr::filter(name %in% fcs_names)

### Random forest model
## Parameters of the model:
## downsample: amount of cells that will be used in each sample in the calculations of the model
## param: which parameters should be taken into account when constructing the model
## plot_fig: whether a figure of the results should be constructed
Model_RF <- Phenoflow::RandomF_FCS(flowData_transformed[fcs_names], Sample_Info_sb, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

# Run Random Forest classifier to predict the Strain based on the single-cell FCM data
#Choose the fcs files in which the model will be used/tested
fcs_topre <- c("Fn_10_3_20191017_145607.fcs",
               "So_Fn_10_3_20191015_162124.fcs",
               "So_Fn_Pi_10_3_20191015_162306.fcs",
               "So_Fn_Pi_Vp_10_3_20191015_162448.fcs",
               "So_10_3_20191017_143923.fcs",
               "Pi_10_3_20191017_145747.fcs",
               "Vp_10_3_20191017_145142.fcs",
               "1_10_3_20191016_141301.fcs",
               "1_10_3_20191017_143044.fcs",
               "2_10_3_20191016_141435.fcs",
               "2_10_3_20191017_143159.fcs")
flowData_topre <- flowData_transformed[fcs_topre]

#### Ignore this comment: The model's performance highly decreased when trasformed data were used - to be discussed with Ruben --> which transformation? asin transformation, normalization, ...?
# Arguments of RandomF_predict function?
test_pred3 <- RandomF_predict(x = Model_RF[[1]], new_data =  flowData_topre, cleanFCS = FALSE)
test_pred3
