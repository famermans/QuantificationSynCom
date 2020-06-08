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
                                  `BL3-A` = asinh(`BL3-A`),
                                  `BL1-H` = asinh(`BL1-H`))
param = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "BL1-H")


#### Extrating metadata from .fcs files ----
sampleNames(flowData_transformed) <- substring(sampleNames(flowData_transformed), 0, nchar(sampleNames(flowData_transformed))-4)       # nchar takes a character vector and returns the number of characters in the vector

# Extracting the actual metadata
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed), "_"), rbind)))
colnames(metadata) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain", "Replicate", "Dilution", "Stain", "Well")


#### Quality control data ----
### BHI2
## Contamination 
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1:296)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC BHI2 (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

## Singlets
xyplot(`BL1-H`~`BL1-A`, data = flowData_transformed[c(1:296)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0,16))),
       axix = axis.default, nbin = 125, main = "QC singlets BHI2 (BL1-A ~ BL1-H)", xlab = "BL1-A", ylab = "BL1-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

### MM
## Contamination
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(297:503)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "QC MM (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

## Singlets
xyplot(`BL1-H`~`BL1-A`, data = flowData_transformed[c(297:503)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0,16))),
       axix = axis.default, nbin = 125, main = "QC singlets MM (BL1-A ~ BL1-H)", xlab = "BL1-A", ylab = "BL1-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)


#### Gating ----

### BHI2
## Gating out extra background based on SSC-BL1
# Constructing gate
sqrcut1 <- matrix(c(5, 12, 14, 10, 5,
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
sqrcut2 <- matrix(c(6.7, 11, 14, 10, 6.7,
                    2, 4.5, 13, 10.5, 5.3), ncol = 2, nrow = 5)
colnames(sqrcut2) <- c("BL1-A", "BL3-A")
polyGate2 <- polygonGate(.gate = sqrcut2, filterId = "Cells BHI2")

# Gating quality check
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(19:22, 29, 30, 37, 55, 85)], filter = polyGate2,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells BHI2", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

### MM
## Gating out extra background based on SSC-BL1
# Constructing gate -> remark: gating was not optimized for this medium
sqrcut3 <- matrix(c(5, 12, 14, 10, 5,
                    2, 5, 13.5, 13.5, 8.5), ncol = 2, nrow = 5)
colnames(sqrcut3) <- c("BL1-A", "SSC-A")
polyGate3 <- polygonGate(.gate = sqrcut3, filterId = "Reduced background BL1-SSC MM")

# Gating quality check
xyplot(`SSC-A`~`BL1-A`, data = flowData_transformed[c(309, 315:318, 325:326, 333, 345)], filter = polyGate3,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating BL1-SSC MM", xlab = "BL1-A", ylab = "SSC-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

## Gating cells based on BL1-BL3
sqrcut4 <- matrix(c(6.7, 13.3, 14.5, 10, 6.7,
                    2, 5, 13.5, 10.5, 5.3), ncol = 2, nrow = 5)
colnames(sqrcut4) <- c("BL1-A", "BL3-A")
polyGate4 <- polygonGate(.gate = sqrcut4, filterId = "Cells MM")

# Gating quality check
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(309, 315:318, 325:326, 333, 345)], filter = polyGate4,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells MM", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

### General gate for both media
## Gating cells based on BL1-BL3
sqrcut5 <- matrix(c(6.7, 13.3, 14.5, 10, 6.7,
                    2, 5, 13.5, 10.5, 5.3), ncol = 2, nrow = 5)
colnames(sqrcut5) <- c("BL1-A", "BL3-A")
polyGate5 <- polygonGate(.gate = sqrcut5, filterId = "Cells")

# Gating quality check

xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(19:22, 29, 30, 37, 55, 85, 309, 315:318, 325:326, 333, 345)], filter = polyGate5,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check all samples
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1:503)], filter = polyGate5,
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
# Calculating concentrations --> since we diluted 1000 times, concentrations will be in cells/µL
cell_concentrations <- data.frame(Samples = flowCore::sampleNames(flowData_transformed),
                                  Strain = metadata$Strain,
                                  Concentration = (TotalCount$true*1000)/vol)


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
