# 1. Clearing workspace, loading libraries, setting seed ----

# Clear environment and set working directory
rm(list = ls())
setwd("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM")

# Load libraries
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

seed <- 777
set.seed(seed)

fourtycolors <- c("#F0E68C", "#FFE4B5", "#BDB76B", "#00CED3", "#20B2AA", "#00FFFF", "#0000FF", "#191970", "#000080", "#FFD700",
                  "#FFFF00", "#CCCC33", "#D2691E", "#FF69B5", "#8A2BE2", "#FF66FF", "#FF1493", "#FF00FF", "#4B0082", "#BA55D3",
                  "#800080", "#800000", "#CC0000", "#F08080", "#FF0000", "#FF3333", "#FF5530", "#FF6C00", "#FFA500", "#4682B4",
                  "#008080", "#008000", "#32CD32", "#336633", "#808000", "#2E8B57", "#00FF7F", "#7CFC00", "#778899", "#ADFF2F")

# Load functions
source(file = "/Projects1/Fabian/paper_theme_fab.R")
source(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/Functions/count_capital_letters.R")


# 2. Load and transform data ----

# Load data
Datapath <- "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/data"
fcsfiles <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = TRUE)
flowData <- flowCore::read.flowSet(files = fcsfiles, transformation = FALSE, emptyValue = F)


# Transformation of data
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


# Extract metadata from filenames
metadata <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_transformed), "_"), rbind)))
colnames(metadata) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain", "Replicate", "Dilution", "Stain", "Well")

metadata$Dilution <- as.numeric(metadata$Dilution)
metadata$Well <- substr(metadata$Well, start = 1, stop = nchar(metadata$Well)-4)

name <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = FALSE)
Sample_Info <- cbind(name, metadata)

# Select data
flowData_transformed_sel <- flowData_transformed[c(1:122, 297:413, 656:696)]
metadata_sel <- metadata[c(1:122, 297:413, 656:696), ]

flowData_transformed_BHI <- flowData_transformed[c(1:122, 574:613, 656:675)]
metadata_BHI <- metadata[c(1:122, 574:613, 656:675), ]
Sample_Info_BHI <- Sample_Info[c(1:122, 574:613, 656:675), ]

rownames(metadata_BHI) <- NULL
rownames(Sample_Info_BHI) <- NULL

# 3. Quality control ----

# Contamination
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

# Singlets
p_singlets <- xyplot(`BL1-H`~`BL1-W`, data = flowData_transformed[c(1:729)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0,16))),
       axix = axis.default, nbin = 125, main = "QC singlets (BL1-A ~ BL1-W)", xlab = "BL1-W", ylab = "BL1-H",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# Stability flow
p_stable_flow <- xyplot(`BL1-H`~`Time`, data = flowData_transformed[c(1:90)],
                        scales = list(y = list(limits = c(4, 14)),
                                      x = list(limits = c(0, 52000))),
                        axix = axis.default, nbin = 125, main = "Stability fluorescence over time (BHI2 12-05-2020)", xlab = "Time [ms]", ylab = "BL1-H",
                        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

# 4. Gating ----

# Gating cells based on BL1-BL3
sqrcut1 <- matrix(c(6.7, 14, 14, 10, 6.7,
                    2, 2, 13, 10.5, 5.3), ncol = 2, nrow = 5)
colnames(sqrcut1) <- c("BL1-A", "BL3-A")
polyGate1 <- polygonGate(.gate = sqrcut1, filterId = "Cells")

# Gating quality check pure cultures
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_BHI[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(4, 15))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells (axenic cultures)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check mixes
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_BHI[c(91, 93, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, 115, 117, 119)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(4, 15))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells (mock mixes)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check co-cultures
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_BHI[c(123, 126, 129, 132, 135, 138, 141, 142)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(4, 15))),
       axis = axis.default, nbin = 125, main = "Quality check gating cells (co-cultures)", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating quality check all samples
p_gating <- xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_BHI, filter = polyGate1,
                   scales = list(y = list(limits = c(0, 15)),
                                 x = list(limits = c(4, 15))),
                   axis = axis.default, nbin = 125, main = "Quality check gating cells BHI", xlab = "BL1-A", ylab = "BL3-A",
                   par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Subset data
flowData_transformed_BHI_gated <- Subset(flowData_transformed_BHI, polyGate1)

# Quality control after gating

# Density plots
p_density_SSC_W <- autoplot(flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], "SSC-W")+
  labs(title = "Density plot SSC-W",
       x = "SSC-W",
       y = "Density")
print(p_density_SSC_W)

p_density_SSC_A <- autoplot(flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], "SSC-A")+
        labs(title = "Density plot SSC-A",
             x = "SSC-A",
             y = "Density")
print(p_density_SSC_A)

p_density_SSC_H <- autoplot(flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], "SSC-H")+
        labs(title = "Density plot SSC-H",
             x = "SSC-H",
             y = "Density")
print(p_density_SSC_H)


# 5. Singlet analysis ----

# Axenic cultures
p_singlets_pure_FSC <- xyplot(`FSC-W`~`FSC-H`, data = flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                              scales = list(y = list(limits = c(0, 10)),
                                            x = list(limits = c(5, 16))),
                              axis = axis.default, nbin = 125, main = "Singlet analysis FSC (axenic cultures)", xlab = "FSC-H", ylab = "FSC-W",
                              par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_singlets_pure_SSC <- xyplot(`SSC-A`~`SSC-H`, data = flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                              scales = list(y = list(limits = c(5, 15)),
                                            x = list(limits = c(5, 15))),
                              axis = axis.default, nbin = 125, main = "Singlet analysis SSC (axenic cultures)", xlab = "SSC-H", ylab = "SSC-A",
                              par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_singlets_pure_SSC_W <- xyplot(`SSC-W`~`SSC-H`, data = flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                              scales = list(y = list(limits = c(0, 10)),
                                            x = list(limits = c(5, 16))),
                              axis = axis.default, nbin = 125, main = "Singlet analysis SSC (axenic cultures)", xlab = "SSC-H", ylab = "SSC-W",
                              par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_singlets_pure_BL1 <- xyplot(`BL1-A`~`BL1-H`, data = flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                              scales = list(y = list(limits = c(5, 15)),
                                            x = list(limits = c(5, 15))),
                              axis = axis.default, nbin = 125, main = "Singlet analysis BL1 (axenic cultures)", xlab = "BL1-H", ylab = "BL1-A",
                              par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_singlets_pure_BL1_W <- xyplot(`BL1-W`~`BL1-H`, data = flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)],
                              scales = list(y = list(limits = c(0, 8)),
                                            x = list(limits = c(5, 15))),
                              axis = axis.default, nbin = 125, main = "Singlet analysis BL1 (axenic cultures)", xlab = "BL1-H", ylab = "BL1-W",
                              par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Gating singlets
sqrcut_singlets <- matrix(c(7, 13., 13.7, 7.7,
                            7, 14, 14, 7), ncol = 2, nrow = 4)
colnames(sqrcut_singlets) <- c("SSC-H", "SSC-A")
polyGateSinglets <- polygonGate(.gate = sqrcut_singlets, filterId = "Singlets")

# Gating quality check
xyplot(`SSC-A`~`SSC-H`, data = flowData_transformed_BHI_gated[c(1, 7, 13, 19, 21, 23, 29, 31, 37, 43, 49, 55, 61, 67, 73, 79, 85)], filter = polyGateSinglets,
       scales = list(y = list(limits = c(5, 15)),
                     x = list(limits = c(5, 15))),
       axis = axis.default, nbin = 125, main = "Gating of singlets (axenic cultures)", xlab = "SSC-H", ylab = "SSC-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

p_QC_singlets <- xyplot(`SSC-A`~`SSC-H`, data = flowData_transformed_BHI_gated, filter = polyGateSinglets,
                        scales = list(y = list(limits = c(5, 15)),
                                      x = list(limits = c(5, 15))),
                        axis = axis.default, nbin = 125, main = "Gating of singlets (axenic cultures)", xlab = "SSC-H", ylab = "SSC-A",
                        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

flowData_transformed_BHI_singlets <- Subset(flowData_transformed_BHI_gated, polyGateSinglets)

singlets <- flowCore::filter(flowData_transformed_BHI_gated, polyGateSinglets)
SingletCount <- summary(singlets)
SingletCount <- toTable(SingletCount)

singlets_full <- cbind(SingletCount, metadata_BHI)
vol <- as.numeric(flowCore::fsApply(flowData_transformed_BHI_gated, FUN = function(x) x@description$`$VOL`))/1000
singlets_full <- cbind(singlets_full, vol)
names(singlets_full)[names(singlets_full) == "count"] <- "events"
names(singlets_full)[names(singlets_full) == "true"] <- "singlets"
names(singlets_full)[names(singlets_full) == "false"] <- "multiplets"

singlets_full <- singlets_full[, c(3:6, 13:16, 18, 19)]
singlets_full$concentration <- ((singlets_full$singlets*singlets_full$Dilution)/singlets_full$vol)*1000

singlets_mean <- aggregate(concentration ~ Strain + Replicate + Timepoint, data = singlets_full, FUN = mean)
names(singlets_mean)[names(singlets_mean) == "concentration"] <- "concentration_singlets"
singlets_percent_mean <- aggregate(percent ~ Strain + Replicate + Timepoint, data = singlets_full, FUN = mean)
singlets_mean$percent <- singlets_percent_mean$percent

# Run 6. before running this part!
singlets_mean$concentration_total <- cell_concentrations_mean$Concentration

singlets_subset <- subset(singlets_mean, Strain == "So" | Strain == "Fn" | Strain == "Pg" | Strain == "Vp" | Strain == "Mix1" | Strain == "Mix2" | Strain == "Mix3" | Strain == "Mix8" | Strain == "Mix9" | Strain == "Co1" | Strain == "Co2" | Strain == "Co3")
singlets_subset2 <- singlets_subset
singlets_subset2$Merge <- paste(singlets_subset2$Strain, singlets_subset2$Replicate, sep = "_")
singlets_theoretical_mocks <- merge(volumes_mocks, singlets_subset2, by = "Merge")
singlets_theoretical_mocks <- singlets_theoretical_mocks[, c(2:6, 10:12)]
colnames(singlets_theoretical_mocks)[colnames(singlets_theoretical_mocks) == "Strain.x"] <- "Strain"
colnames(singlets_theoretical_mocks)[colnames(singlets_theoretical_mocks) == "Replicate.x"] <- "Replicate"
colnames(singlets_theoretical_mocks)[colnames(singlets_theoretical_mocks) == "concentration_singlets"] <- "concentration_axenic_singlets"
colnames(singlets_theoretical_mocks)[colnames(singlets_theoretical_mocks) == "concentration_total"] <- "concentration_axenic_total"

singlets_theoretical_mocks$counts_axenic_singlets <- singlets_theoretical_mocks$Volume*singlets_theoretical_mocks$concentration_axenic_singlets
singlets_theoretical_mocks$concentration_mock_singlets <- singlets_theoretical_mocks$counts_axenic_singlets/singlets_theoretical_mocks$Total_volume
singlets_theoretical_mocks$counts_axenic_total <- singlets_theoretical_mocks$Volume*singlets_theoretical_mocks$concentration_axenic_total
singlets_theoretical_mocks$concentration_mock_total <- singlets_theoretical_mocks$counts_axenic_total/singlets_theoretical_mocks$Total_volume
singlets_theoretical_mocks <- subset(singlets_theoretical_mocks, Mix_ID == "Mix1" | Mix_ID == "Mix2" | Mix_ID == "Mix3" | Mix_ID == "Mix8" | Mix_ID == "Mix9")
row.names(singlets_theoretical_mocks) <- NULL

singlets_theoretical_mocks <- singlets_theoretical_mocks[, c(1, 2, 10, 12)]
saveRDS(object = singlets_theoretical_mocks, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/singlets_theoretical_mocks.rds")

singlets_theoretical_mocks_sumsinglets <- aggregate(concentration_mock_singlets ~ Mix_ID, data = singlets_theoretical_mocks, FUN = sum)
singlets_theoretical_mocks_sumtotal <- aggregate(concentration_mock_total ~ Mix_ID, data = singlets_theoretical_mocks, FUN = sum)
singlets_mocks <- merge(singlets_theoretical_mocks_sumsinglets, singlets_theoretical_mocks_sumtotal, by = "Mix_ID")
singlets_mocks$percentage_theoretical <- (singlets_mocks$concentration_mock_singlets/singlets_mocks$concentration_mock_total)*100
colnames(singlets_mocks)[colnames(singlets_mocks) == "concentration_mock_singlets"] <- "concentration_singlets_theoretical"
colnames(singlets_mocks)[colnames(singlets_mocks) == "concentration_mock_total"] <- "concentration_total_theoretical"

singlets_mean_mocks <- subset(singlets_mean, Strain == "Mix1" | Strain == "Mix2" | Strain == "Mix3" | Strain == "Mix8" | Strain == "Mix9")
colnames(singlets_mean_mocks)[colnames(singlets_mean_mocks) == "Strain"] <- "Mix_ID"
colnames(singlets_mean_mocks)[colnames(singlets_mean_mocks) == "concentration_singlets"] <- "concentration_singlets_measured"
colnames(singlets_mean_mocks)[colnames(singlets_mean_mocks) == "concentration_total"] <- "concentration_total_measured"
colnames(singlets_mean_mocks)[colnames(singlets_mean_mocks) == "percent"] <- "percentage_measured"
singlets_mean_mocks <- singlets_mean_mocks[, c(1, 4:6)]
row.names(singlets_mean_mocks) <- NULL

singlets_mocks <- merge(singlets_mocks, singlets_mean_mocks, by = "Mix_ID")
singlets_melted <- reshape2::melt(singlets_mocks, id.vars = c("Mix_ID"))
singlets_melted$type <- c(rep(c("Theoretical"), 15), rep(c("Measured"), 15))
singlets_melted$type2 <- c(rep(c("Concentration"), 10), rep(c("Percentage"), 5), rep(c("Concentration"), 5), rep(c("Percentage"), 5), rep(c("Concentration"), 5))
singlets_melted$Mix_ID <- gsub("Mix1", "Mock 1", singlets_melted$Mix_ID)
singlets_melted$Mix_ID <- gsub("Mix2", "Mock 2", singlets_melted$Mix_ID)
singlets_melted$Mix_ID <- gsub("Mix3", "Mock 3", singlets_melted$Mix_ID)
singlets_melted$Mix_ID <- gsub("Mix8", "Mock 8", singlets_melted$Mix_ID)
singlets_melted$Mix_ID <- gsub("Mix9", "Mock 9", singlets_melted$Mix_ID)
singlets_melted_concentration <- subset(singlets_melted, type2 == "Concentration")
singlets_melted_concentration$population <- c(rep(c(rep(c("Singlets"), 5), rep(c("Total"), 5)), 2))
singlets_melted_percentage <- subset(singlets_melted, type2 == "Percentage")

p_singlets_mocks_concentration <- ggplot(data = singlets_melted_concentration, aes(x = Mix_ID, y = value, color = population, shape = type)) +
  geom_point(size = 7, alpha = 0.6) +
  labs(y = "Concentration (cells/mL)", x = NULL, color = NULL, shape = NULL) +
  paper_theme_fab +
  scale_color_manual(values = c("Singlets" = "darkred", "Total" = "blue3"),
                     labels = c("Singlets" = "Singlets", "Total"= "All")) +
  scale_shape_manual(values = c("Measured" = 16, "Theoretical" = 17),
                     labels = c("Measured" = "Measured", "Theoretical" = "Calculated"))
print(p_singlets_mocks_concentration)

p_singlets_mocks_percentage <- ggplot(data = singlets_melted_percentage, aes(x = Mix_ID, y = value, color = type)) +
  geom_point(size = 7, alpha = 0.6) +
  labs(y = "Relative abundance (%)", x = NULL, color = NULL) +
  paper_theme_fab +
  scale_y_continuous(limits = c(50, 80)) +
  scale_color_manual(values = c("Measured" = "darkred", "Theoretical" = "blue3"),
                     labels = c("Measured" = "Measured", "Theoretical" = "Calculated"))
print(p_singlets_mocks_percentage)

singlets_cocult <- subset(singlets_mean, Strain == "Co1" | Strain == "Co2" | Strain == "Co3")
singlets_cocult_mean <- aggregate(concentration_singlets ~ Strain + Replicate + Timepoint, data = singlets_cocult, FUN = mean)
singlets_cocult_mean2 <- aggregate(concentration_total ~ Strain + Replicate + Timepoint, data = singlets_cocult, FUN = mean)
singlets_cocult_mean$concentration_total <- singlets_cocult_mean2$concentration_total
rm(singlets_cocult_mean2)
singlets_cocult_mean$percentage <- (singlets_cocult_mean$concentration_singlets/singlets_cocult_mean$concentration_total)*100
singlets_cocult_mean$ID <- paste(singlets_cocult_mean$Strain, singlets_cocult_mean$Replicate, sep = "_")
singlets_cocult_mean$ID <- gsub("Co1_", "Co-culture 1 ", singlets_cocult_mean$ID)
singlets_cocult_mean$ID <- gsub("Co2_", "Co-culture 2 ", singlets_cocult_mean$ID)
singlets_cocult_mean$ID <- gsub("Co3_", "Co-culture 3 ", singlets_cocult_mean$ID)
singlets_cocult_mean$ID <- as.factor(singlets_cocult_mean$ID)

p_singlets_cocult_percentage <- ggplot(data = singlets_cocult_mean, aes(x = Timepoint, y = percentage)) +
  geom_point(size = 7) +
  labs(y = "Relative abundance (%)", x = NULL) +
  facet_wrap(~ ID, ncol = 3) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(30, 90), breaks = seq(30, 90, by = 10))
print(p_singlets_cocult_percentage)

# 6. Cell concentrations ----
## 6.1. Calculate concentrations ----

# Cell counts
cells <- flowCore::filter(flowData_transformed_BHI, polyGate1)
TotalCount <- summary(cells)
TotalCount <- toTable(TotalCount)

# Extracting volumes
# Volumes are in µL
vol <- as.numeric(flowCore::fsApply(flowData_transformed_BHI_gated, FUN = function(x) x@description$`$VOL`))/1000

# Concentrations (cells/mL)
cell_concentrations <- data.frame(Sample_name = flowCore::sampleNames(flowData_transformed_BHI_gated),
                                  Strain = metadata_BHI$Strain,
                                  Replicate = metadata_BHI$Replicate,
                                  Timepoint = metadata_BHI$Timepoint,
                                  Concentration = ((TotalCount$true*metadata_BHI$Dilution)/vol)*1000)


## 6.2. Calculations theoretical concentrations mocks ----
# Calculate mean cell concentration per replicate
cell_concentrations_mean <- aggregate(Concentration ~ Strain + Replicate + Timepoint, data = cell_concentrations, FUN = mean)
cell_concentrations_mean$Merge <- paste(cell_concentrations_mean$Strain, cell_concentrations_mean$Replicate, sep = "_")

volumes_mocks <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/20200512_Mocks.csv", header = T, sep = ";", stringsAsFactors = T)
volumes_mocks$Merge <- paste(volumes_mocks$Strain, volumes_mocks$Replicate, sep = "_")

theoretical_mocks <- merge(volumes_mocks, cell_concentrations_mean, by = "Merge")
theoretical_mocks <- theoretical_mocks[, c(2:6, 10)]
colnames(theoretical_mocks)[colnames(theoretical_mocks) == "Strain.x"] <- "Strain"
colnames(theoretical_mocks)[colnames(theoretical_mocks) == "Replicate.x"] <- "Replicate"
colnames(theoretical_mocks)[colnames(theoretical_mocks) == "Concentration"] <- "Concentration_Axenic"

theoretical_mocks$Counts <- theoretical_mocks$Volume*theoretical_mocks$Concentration_Axenic
theoretical_mocks$Concentration_Mock <- theoretical_mocks$Counts/theoretical_mocks$Total_volume

theoretical_mocks <- theoretical_mocks[, c(1, 2, 8)]
saveRDS(object = theoretical_mocks, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/theoretical_mocks.rds")


# 7. Phenotypic diversity analysis ----
# Normalization of data
summary <- fsApply(x = flowData_transformed_BHI_gated, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max = max(summary[, "BL1-H"])
mytrans <- function(x) x/max

flowData_transformed_BHI_norm <- transform(flowData_transformed_BHI_gated,
                                    `FSC-A` = mytrans(`FSC-A`), 
                                    `SSC-A` = mytrans(`SSC-A`), 
                                    `BL1-A` = mytrans(`BL1-A`), 
                                    `BL3-A` = mytrans(`BL3-A`),
                                    `FSC-H` = mytrans(`FSC-H`), 
                                    `SSC-H` = mytrans(`SSC-H`), 
                                    `BL1-H` = mytrans(`BL1-H`), 
                                    `BL3-H` = mytrans(`BL3-H`))

# Calculating fingerprint with bw = 0.01
#fbasis <- flowBasis(flowData_transformed_BHI_norm, param, nbin = 128, 
#                    bw = 0.01, normalize = function(x) x)
#saveRDS(object = fbasis, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/fbasis.rds")
fbasis <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/fbasis.rds")

# Calculate ecological parameters
#Diversity.fbasis <- Diversity(fbasis, d = 3, plot = FALSE, R = 999)
#Evenness.fbasis <- Evenness(fbasis, d = 3, plot = FALSE)
#Structural.organization.fbasis <- So(fbasis, d = 3, plot = FALSE)
#Coef.var.fbasis <- CV(fbasis, d = 3, plot = FALSE)
#saveRDS(object = Diversity.fbasis, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/Diversity_fbasis.rds")
Diversity.fbasis <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/Diversity_fbasis.rds")

# Plot ecological parameters
p_alphadiv <- ggplot(data = Diversity.fbasis, aes(x = as.character(metadata_BHI$Strain), y = D2))+
        geom_point(size = 4, alpha = 0.7)+
        geom_line()+
        facet_grid(metadata_BHI$Replicate ~ ., scales = "free")+       # Creates different subplots
        theme_bw()+
        labs(y = "Phenotypic diversity (D2)", x = "Sample", title = "Phenotypic diversity analysis")+
        geom_errorbar(aes(ymin = D2-sd.D2, ymax = D2+sd.D2), width = 0.05)
print(p_alphadiv)

# Beta-diversity assessment of fingerprint
beta.div <- beta_div_fcm(fbasis, ord.type = "PCoA")

# Plot ordination
p_betadiv <- plot_beta_fcm(beta.div, color = metadata_BHI$Strain, shape = as.factor(metadata_BHI$Replicate), labels = list("Strain", "Replicate")) +
  geom_point(size = 8, alpha = 0.5) +
  scale_color_manual(values = fourtycolors) +
  scale_shape_manual(values = c(0, 1, 2, 5)) +
  theme_bw()
print(p_betadiv)

var.pcoa <- vegan::eigenvals(beta.div)/sum(vegan::eigenvals(beta.div))
PcoA <- as.data.frame(beta.div$points)
names(PcoA)[1:2] <- c("Axis1", "Axis2")

p_betadiv_2 <- ggplot(PcoA, aes(x = Axis1, y = Axis2, color = metadata_BHI$Strain, shape = metadata_BHI$Replicate))+
  geom_point(size = 8, alpha = 0.7) +
  labs(title = "Ordination of phenotypic fingerprints",
       x = paste0("Axis1 (", round(100 * var.pcoa[1], 1), "%)"),
       y = paste0("Axis2 (", round(100 * var.pcoa[2], 1), "%)"),
       color = "Strain",
       shape = "Replicate") +
  scale_color_manual(values = fourtycolors) +
  scale_shape_manual(values = c(15:18)) +
  theme_bw()
print(p_betadiv_2)


# 8. Random Forest ----
## 8.1. Training with one replicate ----

# Define parameters on which RF will build its model
paramRF = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "FSC-H", "SSC-H", "BL1-H", "BL3-H")
paramRF2 = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-W", "SSC-W", "BL1-W", "BL3-W")

# Extract total number of events per fcs file in order to calculate accuracy of model to predict strain
vol2 <- data.frame(Sample_name = flowCore::sampleNames(flowData_transformed_BHI),
                   Volume = vol)
TotalCount2 <- left_join(TotalCount, vol2, by = c("sample" = "Sample_name"))
TotalCount2 <- left_join(TotalCount2, cell_concentrations, by = c("sample" = "Sample_name"))
write.csv2(file = "TotalCount.csv", TotalCount2)

### Model for So and Fn in BHI2
# Select the fcs files based on which the model will be trained --> So (replicate A), Fn (replicate C) grown in BHI2
fcs_names_SoFn <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                    "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs")

Sample_Info_SoFn <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFn)
#Model_RF_SoFn <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoFn], Sample_Info_SoFn, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoFn, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFn.rds")
Model_RF_SoFn <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFn.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoFn <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
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

flowData_topre_BHI_SoFn <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoFn]
test_pred_BHI_SoFn <- RandomF_predict(x = Model_RF_SoFn[[1]], new_data =  flowData_topre_BHI_SoFn, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFn <- left_join(test_pred_BHI_SoFn, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoFn <- left_join(test_pred_BHI2_SoFn, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoFn <- test_pred_BHI2_SoFn %>%
  mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoFn.csv", test_pred_BHI2_SoFn)


### Model for So, Fn and Pg in BHI2
# Sample selection So
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_BHI[c(61:66)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(4, 15))),
       axis = axis.default, nbin = 125, main = "Quality check So", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection Fn
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_BHI[c(23:28)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(4, 15))),
       axis = axis.default, nbin = 125, main = "Quality check Fn", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Sample selection Pg
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_BHI[c(31:36)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(4, 15))),
       axis = axis.default, nbin = 125, main = "Quality check Pg", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Fn (replicate C), Pg (replicate B) grown in BHI2
fcs_names_SoFnPg <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs")

Sample_Info_SoFnPg <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPg)

# Train random forest classifier
#Model_RF_SoFnPg <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoFnPg], Sample_Info_SoFnPg, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
# Test whether normalization of data influences performance of the model -> Performance is identical
#Model_RF_SoFnPg_N <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_norm[fcs_names_SoFnPg], Sample_Info_SoFnPg, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoFnPg, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPg.rds")
Model_RF_SoFnPg <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPg.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoFnPg <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
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

flowData_topre_BHI_SoFnPg <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoFnPg]
test_pred_BHI_SoFnPg <- RandomF_predict(x = Model_RF_SoFnPg[[1]], new_data =  flowData_topre_BHI_SoFnPg, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFnPg <- left_join(test_pred_BHI_SoFnPg, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoFnPg <- left_join(test_pred_BHI2_SoFnPg, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoFnPg <- test_pred_BHI2_SoFnPg %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoFnPg.csv", test_pred_BHI2_SoFnPg)


### Model for So, Fn and Pi in BHI2 for comparison in performance compared to FACSVerse
# Sample selection Pi
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(37:42)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "Quality check Pi", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Fn (replicate C), Pi (replicate A) grown in BHI2
fcs_names_SoFnPi <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                      "20200512_Fabian_14strainID_BHI2_24h_Pi_A_1000_SG_A1.fcs")

Sample_Info_SoFnPi <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPi)
#Model_RF_SoFnPi <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoFnPi], Sample_Info_SoFnPi, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoFnPi, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPi.rds")
Model_RF_SoFnPi <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPi.rds")

### Model for So, Fn, Pg and Vp
# Sample selection Vp
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(85:90)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 15)),
                     x = list(limits = c(4, 15))),
       axis = axis.default, nbin = 125, main = "Quality check Vp", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

# Select the fcs files based on which the model will be trained --> So (replicate A), Fn (replicate C), Pg (replicate B), Vp (replicate B) grown in BHI2
fcs_names_SoFnPgVp <- c("20200512_Fabian_14strainID_BHI2_24h_So_A_1000_SG_D1.fcs",
                        "20200512_Fabian_14strainID_BHI2_24h_Fn_C_1000_SG_B5.fcs",
                        "20200512_Fabian_14strainID_BHI2_24h_Pg_B_1000_SG_E10.fcs",
                        "20200512_Fabian_14strainID_BHI2_24h_Vp_B_1000_SG_A10.fcs")

Sample_Info_SoFnPgVp <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPgVp)
#Model_RF_SoFnPgVp <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoFnPgVp], Sample_Info_SoFnPgVp, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoFnPgVp, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPgVp.rds")
Model_RF_SoFnPgVp <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPgVp.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoFnPgVp <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
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

flowData_topre_BHI_SoFnPgVp <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoFnPgVp]
test_pred_BHI_SoFnPgVp <- RandomF_predict(x = Model_RF_SoFnPgVp[[1]], new_data =  flowData_topre_BHI_SoFnPgVp, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFnPgVp <- left_join(test_pred_BHI_SoFnPgVp, vol2, by = c("Sample" = "Sample_name"))
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
#Model_RF_SoSsalSsanSmiSgVp <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoSsalSsanSmiSgVp], Sample_Info_SoSsalSsanSmiSgVp, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoSsalSsanSmiSgVp, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgVp.rds")
Model_RF_SoSsalSsanSmiSgVp <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgVp.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoSsalSsanSmiSgVp <- c("20200617_Fabian_14strainID_BHI2_24h_Co4_A_1000_SG_C3.fcs",
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

flowData_topre_BHI_SoSsalSsanSmiSgVp <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoSsalSsanSmiSgVp]
test_pred_BHI_SoSsalSsanSmiSgVp <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgVp[[1]], new_data =  flowData_topre_BHI_SoSsalSsanSmiSgVp, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgVp <- left_join(test_pred_BHI_SoSsalSsanSmiSgVp, vol2, by = c("Sample" = "Sample_name"))
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
#Model_RF_SoSsalSsanSmiSgVpAvAn <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoSsalSsanSmiSgVpAvAn], Sample_Info_SoSsalSsanSmiSgVpAvAn, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
# REMARK: not enough cells for An to train model --> left out of the calculations
#saveRDS(object = Model_RF_SoSsalSsanSmiSgVpAvAn, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgVpAvAn.rds")
Model_RF_SoSsalSsanSmiSgVpAvAn <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgVpAvAn.rds")


## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoSsalSsanSmiSgVpAvAn <- c("20200617_Fabian_14strainID_BHI2_24h_Co4_A_1000_SG_C3.fcs",
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

flowData_topre_BHI_SoSsalSsanSmiSgVpAvAn <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoSsalSsanSmiSgVpAvAn]
test_pred_BHI_SoSsalSsanSmiSgVpAvAn <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgVpAvAn[[1]], new_data = flowData_topre_BHI_SoSsalSsanSmiSgVpAvAn, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgVpAvAn <- left_join(test_pred_BHI_SoSsalSsanSmiSgVpAvAn, vol2, by = c("Sample" = "Sample_name"))
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
#Model_RF_SoSsalSsanSmiSgSmuSsob <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoSsalSsanSmiSgSmuSsob], Sample_Info_SoSsalSsanSmiSgSmuSsob, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoSsalSsanSmiSgSmuSsob, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgSmuSsob.rds")
Model_RF_SoSsalSsanSmiSgSmuSsob <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgSmuSsob.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoSsalSsanSmiSgSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix5_NA_1000_SG_C8.fcs",
                                           "20200512_Fabian_14strainID_BHI2_NA_Mix5_NA_1000_SG_C9.fcs")

flowData_topre_BHI_SoSsalSsanSmiSgSmuSsob <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoSsalSsanSmiSgSmuSsob]
test_pred_BHI_SoSsalSsanSmiSgSmuSsob <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgSmuSsob[[1]], new_data = flowData_topre_BHI_SoSsalSsanSmiSgSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgSmuSsob <- left_join(test_pred_BHI_SoSsalSsanSmiSgSmuSsob, vol2, by = c("Sample" = "Sample_name"))
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
#Model_RF_AaFnPgPiSmuSsob <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_AaFnPgPiSmuSsob], Sample_Info_AaFnPgPiSmuSsob, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_AaFnPgPiSmuSsob, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_AaFnPgPiSmuSsob.rds")
Model_RF_AaFnPgPiSmuSsob <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_AaFnPgPiSmuSsob.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_AaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix6_NA_1000_SG_C10.fcs",
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

flowData_topre_BHI_AaFnPgPiSmuSsob <- flowData_transformed_BHI_gated[fcs_topre_BHI_AaFnPgPiSmuSsob]
test_pred_BHI_AaFnPgPiSmuSsob <- RandomF_predict(x = Model_RF_AaFnPgPiSmuSsob[[1]], new_data = flowData_topre_BHI_AaFnPgPiSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_AaFnPgPiSmuSsob <- left_join(test_pred_BHI_AaFnPgPiSmuSsob, vol2, by = c("Sample" = "Sample_name"))
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
#Model_RF_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob], Sample_Info_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob.rds")
Model_RF_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix7_NA_1000_SG_D8.fcs",
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

flowData_topre_BHI_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob]
test_pred_BHI_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob[[1]], new_data = flowData_topre_BHI_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob <- left_join(test_pred_BHI_SoSsalSsanSmiSgVpAaFnPgPiSmuSsob, vol2, by = c("Sample" = "Sample_name"))
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
#Model_RF_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- Phenoflow::RandomF_FCS(flowData_transformed_BHI_gated[fcs_names_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob], Sample_Info_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, target_label = "Strain", downsample = 10000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
# REMARK: not enough cells for An to train model --> left out of the calculations
#saveRDS(object = Model_RF_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob.rds")
Model_RF_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob.rds")


## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix7_NA_1000_SG_D8.fcs",
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

flowData_topre_BHI_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob]
test_pred_BHI_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- RandomF_predict(x = Model_RF_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob[[1]], new_data = flowData_topre_BHI_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- left_join(test_pred_BHI_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- left_join(test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob <- test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob %>%
        mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob.csv", test_pred_BHI2_SoSsalSsanSmiSgVpAvAnAaFnPgPiSmuSsob)


## 8.2. Training with replicates ----
# downsample is set to 12328 events, as this is the lowest number of events for a pooled strain (An) --> make models comparable

### Model for So and Fn in BHI2
# Select the fcs files based on which the model will be trained --> So, Fn grown in BHI2
flowData_pooled_SoFn <- FCS_pool(x = flowData_transformed_BHI_gated,
                                     stub = c("20200512_Fabian_14strainID_BHI2_24h_So",
                                              "20200512_Fabian_14strainID_BHI2_24h_Fn"))

metadata_pooled_SoFn <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_pooled_SoFn), "_"), rbind)))
colnames(metadata_pooled_SoFn) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain")

name_pooled_SoFn <- flowCore::sampleNames(flowData_pooled_SoFn)
Sample_Info_pooled_SoFn <- cbind(name_pooled_SoFn, metadata_pooled_SoFn)

# Plot not working for some reason
# xyplot(`BL3-A`~`BL1-A`, data = flowData_pooled_SoFn,
#        scales = list(y = list(limits = c(0, 15)),
#                      x = list(limits = c(5, 15))),
#        axis = axis.default, nbin = 125, main = "Test pooled samples", xlab = "BL1-A", ylab = "BL3-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

cells_SoFn_pooled <- flowCore::filter(flowData_pooled_SoFn, polyGate1)
TotalCount_SoFn_pooled <- summary(cells_SoFn_pooled)
TotalCount_SoFn_pooled <- toTable(TotalCount_SoFn_pooled)

fcs_names_SoFn_pooled <- c("20200512_Fabian_14strainID_BHI2_24h_So",
                           "20200512_Fabian_14strainID_BHI2_24h_Fn")

#Model_RF_SoFn_pooled <- Phenoflow::RandomF_FCS(flowData_pooled_SoFn[fcs_names_SoFn_pooled], Sample_Info_pooled_SoFn, sample_col = "name_pooled_SoFn", target_label = "Strain", downsample = 12328, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoFn_pooled, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFn_pooled.rds")
Model_RF_SoFn_pooled <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFn_pooled.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoFn <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
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

flowData_topre_BHI_SoFn <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoFn]
test_pred_BHI_SoFn_pooled <- RandomF_predict(x = Model_RF_SoFn_pooled[[1]], new_data = flowData_topre_BHI_SoFn, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFn_pooled <- left_join(test_pred_BHI_SoFn_pooled, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoFn_pooled <- left_join(test_pred_BHI2_SoFn_pooled, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoFn_pooled <- test_pred_BHI2_SoFn_pooled %>%
  mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoFn_pooled.csv", test_pred_BHI2_SoFn_pooled)


### Model for So, Fn and Pg in BHI2
# Select the fcs files based on which the model will be trained --> So, Fn, Pg grown in BHI2
flowData_pooled_SoFnPg <- FCS_pool(x = flowData_transformed_BHI_gated,
                                     stub = c("20200512_Fabian_14strainID_BHI2_24h_So",
                                              "20200512_Fabian_14strainID_BHI2_24h_Fn",
                                              "20200512_Fabian_14strainID_BHI2_24h_Pg"))

metadata_pooled_SoFnPg <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_pooled_SoFnPg), "_"), rbind)))
colnames(metadata_pooled_SoFnPg) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain")

name_pooled_SoFnPg <- flowCore::sampleNames(flowData_pooled_SoFnPg)
Sample_Info_pooled_SoFnPg <- cbind(name_pooled_SoFnPg, metadata_pooled_SoFnPg)

# Plot not working for some reason
# xyplot(`BL3-A`~`BL1-A`, data = flowData_pooled_SoFnPg,
#        scales = list(y = list(limits = c(0, 15)),
#                      x = list(limits = c(5, 15))),
#        axis = axis.default, nbin = 125, main = "Test pooled samples", xlab = "BL1-A", ylab = "BL3-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

cells_SoFnPg_pooled <- flowCore::filter(flowData_pooled_SoFnPg, polyGate1)
TotalCount_SoFnPg_pooled <- summary(cells_SoFnPg_pooled)
TotalCount_SoFnPg_pooled <- toTable(TotalCount_SoFnPg_pooled)

fcs_names_SoFnPg_pooled <- c("20200512_Fabian_14strainID_BHI2_24h_So",
                             "20200512_Fabian_14strainID_BHI2_24h_Fn",
                             "20200512_Fabian_14strainID_BHI2_24h_Pg")

#Model_RF_SoFnPg_pooled <- Phenoflow::RandomF_FCS(flowData_pooled_SoFnPg[fcs_names_SoFnPg_pooled], Sample_Info_pooled_SoFnPg, sample_col = "name_pooled_SoFnPg", target_label = "Strain", downsample = 12328, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoFnPg_pooled, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPg_pooled.rds")
Model_RF_SoFnPg_pooled <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPg_pooled.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoFnPg <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
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

flowData_topre_BHI_SoFnPg <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoFnPg]
test_pred_BHI_SoFnPg_pooled <- RandomF_predict(x = Model_RF_SoFnPg_pooled[[1]], new_data = flowData_topre_BHI_SoFnPg, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFnPg_pooled <- left_join(test_pred_BHI_SoFnPg_pooled, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoFnPg_pooled <- left_join(test_pred_BHI2_SoFnPg_pooled, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoFnPg_pooled <- test_pred_BHI2_SoFnPg_pooled %>%
  mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoFnPg_pooled.csv", test_pred_BHI2_SoFnPg_pooled)


### Model for So, Fn, Pg and Vp in BHI2
# Select the fcs files based on which the model will be trained --> So, Fn, Pg, Vp grown in BHI2
flowData_pooled_SoFnPgVp <- FCS_pool(x = flowData_transformed_BHI_gated,
                            stub = c("20200512_Fabian_14strainID_BHI2_24h_So",
                                     "20200512_Fabian_14strainID_BHI2_24h_Fn",
                                     "20200512_Fabian_14strainID_BHI2_24h_Pg",
                                     "20200512_Fabian_14strainID_BHI2_24h_Vp"))

metadata_pooled_SoFnPgVp <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_pooled_SoFnPgVp), "_"), rbind)))
colnames(metadata_pooled_SoFnPgVp) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain")

name_pooled_SoFnPgVp <- flowCore::sampleNames(flowData_pooled_SoFnPgVp)
Sample_Info_pooled_SoFnPgVp <- cbind(name_pooled_SoFnPgVp, metadata_pooled_SoFnPgVp)

# Plot not working for some reason
# xyplot(`BL3-A`~`BL1-A`, data = flowData_pooled_SoFnPgVp,
#        scales = list(y = list(limits = c(0, 15)),
#                      x = list(limits = c(5, 15))),
#        axis = axis.default, nbin = 125, main = "Test pooled samples", xlab = "BL1-A", ylab = "BL3-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

cells_SoFnPgVp_pooled <- flowCore::filter(flowData_pooled_SoFnPgVp, polyGate1)
TotalCount_SoFnPgVp_pooled <- summary(cells_SoFnPgVp_pooled)
TotalCount_SoFnPgVp_pooled <- toTable(TotalCount_SoFnPgVp_pooled)

fcs_names_SoFnPgVp_pooled <- c("20200512_Fabian_14strainID_BHI2_24h_So",
                               "20200512_Fabian_14strainID_BHI2_24h_Fn",
                               "20200512_Fabian_14strainID_BHI2_24h_Pg",
                               "20200512_Fabian_14strainID_BHI2_24h_Vp")

#Model_RF_SoFnPgVp_pooled <- Phenoflow::RandomF_FCS(flowData_pooled_SoFnPgVp[fcs_names_SoFnPgVp_pooled], Sample_Info_pooled_SoFnPgVp, sample_col = "name_pooled_SoFnPgVp", target_label = "Strain", downsample = 12328, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SoFnPgVp_pooled, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPgVp_pooled.rds")
Model_RF_SoFnPgVp_pooled <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SoFnPgVp_pooled.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SoFnPgVp <- c("20200617_Fabian_14strainID_BHI2_24h_Co1_A_1000_SG_A1.fcs",
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

flowData_topre_BHI_SoFnPgVp <- flowData_transformed_BHI_gated[fcs_topre_BHI_SoFnPgVp]
test_pred_BHI_SoFnPgVp_pooled <- RandomF_predict(x = Model_RF_SoFnPgVp_pooled[[1]], new_data = flowData_topre_BHI_SoFnPgVp, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SoFnPgVp_pooled <- left_join(test_pred_BHI_SoFnPgVp_pooled, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SoFnPgVp_pooled <- left_join(test_pred_BHI2_SoFnPgVp_pooled, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SoFnPgVp_pooled <- test_pred_BHI2_SoFnPgVp_pooled %>%
  mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSoFnPgVp_pooled.csv", test_pred_BHI2_SoFnPgVp_pooled)


## Models for artificial mixtures

# Mix4
### Model for An, Av, Sg, Smi, So, Ssal, Ssan and Vp in BHI2
# Select the fcs files based on which the model will be trained
flowData_pooled_AnAvSgSmiSoSsalSsanVp <- FCS_pool(x = flowData_transformed_BHI_gated,
                                                  stub = c("20200512_Fabian_14strainID_BHI2_24h_An",
                                                           "20200512_Fabian_14strainID_BHI2_24h_Av",
                                                           "20200512_Fabian_14strainID_BHI2_24h_Sg",
                                                           "20200512_Fabian_14strainID_BHI2_24h_Smi",
                                                           "20200512_Fabian_14strainID_BHI2_24h_So",
                                                           "20200512_Fabian_14strainID_BHI2_24h_Ssal",
                                                           "20200512_Fabian_14strainID_BHI2_24h_Ssan",
                                                           "20200512_Fabian_14strainID_BHI2_24h_Vp"))

metadata_pooled_AnAvSgSmiSoSsalSsanVp <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_pooled_AnAvSgSmiSoSsalSsanVp), "_"), rbind)))
colnames(metadata_pooled_AnAvSgSmiSoSsalSsanVp) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain")

name_pooled_AnAvSgSmiSoSsalSsanVp <- flowCore::sampleNames(flowData_pooled_AnAvSgSmiSoSsalSsanVp)
Sample_Info_pooled_AnAvSgSmiSoSsalSsanVp <- cbind(name_pooled_AnAvSgSmiSoSsalSsanVp, metadata_pooled_AnAvSgSmiSoSsalSsanVp)

# Plot not working for some reason
# xyplot(`BL3-A`~`BL1-A`, data = flowData_pooled_AnAvSgSmiSoSsalSsanVp,
#        scales = list(y = list(limits = c(0, 15)),
#                      x = list(limits = c(5, 15))),
#        axis = axis.default, nbin = 125, main = "Test pooled samples", xlab = "BL1-A", ylab = "BL3-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

cells_AnAvSgSmiSoSsalSsanVp_pooled <- flowCore::filter(flowData_pooled_AnAvSgSmiSoSsalSsanVp, polyGate1)
TotalCount_AnAvSgSmiSoSsalSsanVp_pooled <- summary(cells_AnAvSgSmiSoSsalSsanVp_pooled)
TotalCount_AnAvSgSmiSoSsalSsanVp_pooled <- toTable(TotalCount_AnAvSgSmiSoSsalSsanVp_pooled)

fcs_names_AnAvSgSmiSoSsalSsanVp_pooled <- c("20200512_Fabian_14strainID_BHI2_24h_An",
                                            "20200512_Fabian_14strainID_BHI2_24h_Av",
                                            "20200512_Fabian_14strainID_BHI2_24h_Sg",
                                            "20200512_Fabian_14strainID_BHI2_24h_Smi",
                                            "20200512_Fabian_14strainID_BHI2_24h_So",
                                            "20200512_Fabian_14strainID_BHI2_24h_Ssal",
                                            "20200512_Fabian_14strainID_BHI2_24h_Ssan",
                                            "20200512_Fabian_14strainID_BHI2_24h_Vp")

#Model_RF_AnAvSgSmiSoSsalSsanVp_pooled <- Phenoflow::RandomF_FCS(flowData_pooled_AnAvSgSmiSoSsalSsanVp[fcs_names_AnAvSgSmiSoSsalSsanVp_pooled], Sample_Info_pooled_AnAvSgSmiSoSsalSsanVp, sample_col = "name_pooled_AnAvSgSmiSoSsalSsanVp", target_label = "Strain", downsample = 12328, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_AnAvSgSmiSoSsalSsanVp_pooled, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_AnAvSgSmiSoSsalSsanVp_pooled.rds")
Model_RF_AnAvSgSmiSoSsalSsanVp_pooled <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_AnAvSgSmiSoSsalSsanVppooled.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_AnAvSgSmiSoSsalSsanVp <- c("20200512_Fabian_14strainID_BHI2_NA_Mix4_NA_1000_SG_B10.fcs",
                                         "20200512_Fabian_14strainID_BHI2_NA_Mix4_NA_1000_SG_B11.fcs")

flowData_topre_BHI_AnAvSgSmiSoSsalSsanVp <- flowData_transformed_BHI_gated[fcs_topre_BHI_AnAvSgSmiSoSsalSsanVp]
test_pred_BHI_AnAvSgSmiSoSsalSsanVp_pooled <- RandomF_predict(x = Model_RF_AnAvSgSmiSoSsalSsanVp_pooled[[1]], new_data = flowData_topre_BHI_AnAvSgSmiSoSsalSsanVp, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_AnAvSgSmiSoSsalSsanVp_pooled <- left_join(test_pred_BHI_AnAvSgSmiSoSsalSsanVp_pooled, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_AnAvSgSmiSoSsalSsanVp_pooled <- left_join(test_pred_BHI2_AnAvSgSmiSoSsalSsanVp_pooled, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_AnAvSgSmiSoSsalSsanVp_pooled <- test_pred_BHI2_AnAvSgSmiSoSsalSsanVp_pooled %>%
  mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsAnAvSgSmiSoSsalSsanVp_pooled.csv", test_pred_BHI2_AnAvSgSmiSoSsalSsanVp_pooled)


# Mix5
### Model for Sg, Smi, Smu, So, Ssal, Ssan and Ssob in BHI2
# Select the fcs files based on which the model will be trained
flowData_pooled_SgSmiSmuSoSsalSsanSsob <- FCS_pool(x = flowData_transformed_BHI_gated,
                                                   stub = c("20200512_Fabian_14strainID_BHI2_24h_Sg",
                                                            "20200512_Fabian_14strainID_BHI2_24h_Smi",
                                                            "20200512_Fabian_14strainID_BHI2_24h_Smu",
                                                            "20200512_Fabian_14strainID_BHI2_24h_So",
                                                            "20200512_Fabian_14strainID_BHI2_24h_Ssal",
                                                            "20200512_Fabian_14strainID_BHI2_24h_Ssan",
                                                            "20200512_Fabian_14strainID_BHI2_24h_Ssob"))

metadata_pooled_SgSmiSmuSoSsalSsanSsob <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_pooled_SgSmiSmuSoSsalSsanSsob), "_"), rbind)))
colnames(metadata_pooled_SgSmiSmuSoSsalSsanSsob) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain")

name_pooled_SgSmiSmuSoSsalSsanSsob <- flowCore::sampleNames(flowData_pooled_SgSmiSmuSoSsalSsanSsob)
Sample_Info_pooled_SgSmiSmuSoSsalSsanSsob <- cbind(name_pooled_SgSmiSmuSoSsalSsanSsob, metadata_pooled_SgSmiSmuSoSsalSsanSsob)

# Plot not working for some reason
# xyplot(`BL3-A`~`BL1-A`, data = flowData_pooled_SgSmiSmuSoSsalSsanSsob,
#        scales = list(y = list(limits = c(0, 15)),
#                      x = list(limits = c(5, 15))),
#        axis = axis.default, nbin = 125, main = "Test pooled samples", xlab = "BL1-A", ylab = "BL3-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

cells_SgSmiSmuSoSsalSsanSsob_pooled <- flowCore::filter(flowData_pooled_SgSmiSmuSoSsalSsanSsob, polyGate1)
TotalCount_SgSmiSmuSoSsalSsanSsob_pooled <- summary(cells_SgSmiSmuSoSsalSsanSsob_pooled)
TotalCount_SgSmiSmuSoSsalSsanSsob_pooled <- toTable(TotalCount_SgSmiSmuSoSsalSsanSsob_pooled)

fcs_names_SgSmiSmuSoSsalSsanSsob_pooled <- c("20200512_Fabian_14strainID_BHI2_24h_Sg",
                                             "20200512_Fabian_14strainID_BHI2_24h_Smi",
                                             "20200512_Fabian_14strainID_BHI2_24h_Smu",
                                             "20200512_Fabian_14strainID_BHI2_24h_So",
                                             "20200512_Fabian_14strainID_BHI2_24h_Ssal",
                                             "20200512_Fabian_14strainID_BHI2_24h_Ssan",
                                             "20200512_Fabian_14strainID_BHI2_24h_Ssob")

#Model_RF_SgSmiSmuSoSsalSsanSsob_pooled <- Phenoflow::RandomF_FCS(flowData_pooled_SgSmiSmuSoSsalSsanSsob[fcs_names_SgSmiSmuSoSsalSsanSsob_pooled], Sample_Info_pooled_SgSmiSmuSoSsalSsanSsob, sample_col = "name_pooled_SgSmiSmuSoSsalSsanSsob", target_label = "Strain", downsample = 12328, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_SgSmiSmuSoSsalSsanSsob_pooled, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SgSmiSmuSoSsalSsanSsob_pooled.rds")
Model_RF_SgSmiSmuSoSsalSsanSsob_pooled <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_SgSmiSmuSoSsalSsanSsobpooled.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_SgSmiSmuSoSsalSsanSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix5_NA_1000_SG_C8.fcs",
                                          "20200512_Fabian_14strainID_BHI2_NA_Mix5_NA_1000_SG_C9.fcs")

flowData_topre_BHI_SgSmiSmuSoSsalSsanSsob <- flowData_transformed_BHI_gated[fcs_topre_BHI_SgSmiSmuSoSsalSsanSsob]
test_pred_BHI_SgSmiSmuSoSsalSsanSsob_pooled <- RandomF_predict(x = Model_RF_SgSmiSmuSoSsalSsanSsob_pooled[[1]], new_data = flowData_topre_BHI_SgSmiSmuSoSsalSsanSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_SgSmiSmuSoSsalSsanSsob_pooled <- left_join(test_pred_BHI_SgSmiSmuSoSsalSsanSsob_pooled, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_SgSmiSmuSoSsalSsanSsob_pooled <- left_join(test_pred_BHI2_SgSmiSmuSoSsalSsanSsob_pooled, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_SgSmiSmuSoSsalSsanSsob_pooled <- test_pred_BHI2_SgSmiSmuSoSsalSsanSsob_pooled %>%
  mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsSgSmiSmuSoSsalSsanSsob_pooled.csv", test_pred_BHI2_SgSmiSmuSoSsalSsanSsob_pooled)


# Mix6
### Model for Aa, Fn, Pg, Pi, Smu and Ssob in BHI2
# Select the fcs files based on which the model will be trained
flowData_pooled_AaFnPgPiSmuSsob <- FCS_pool(x = flowData_transformed_BHI_gated,
                                            stub = c("20200512_Fabian_14strainID_BHI2_24h_Aa",
                                                     "20200512_Fabian_14strainID_BHI2_24h_Fn",
                                                     "20200512_Fabian_14strainID_BHI2_24h_Pg",
                                                     "20200512_Fabian_14strainID_BHI2_24h_Pi",
                                                     "20200512_Fabian_14strainID_BHI2_24h_Smu",
                                                     "20200512_Fabian_14strainID_BHI2_24h_Ssob"))

metadata_pooled_AaFnPgPiSmuSsob <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_pooled_AaFnPgPiSmuSsob), "_"), rbind)))
colnames(metadata_pooled_AaFnPgPiSmuSsob) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain")

name_pooled_AaFnPgPiSmuSsob <- flowCore::sampleNames(flowData_pooled_AaFnPgPiSmuSsob)
Sample_Info_pooled_AaFnPgPiSmuSsob <- cbind(name_pooled_AaFnPgPiSmuSsob, metadata_pooled_AaFnPgPiSmuSsob)

# Plot not working for some reason
# xyplot(`BL3-A`~`BL1-A`, data = flowData_pooled_AaFnPgPiSmuSsob,
#        scales = list(y = list(limits = c(0, 15)),
#                      x = list(limits = c(5, 15))),
#        axis = axis.default, nbin = 125, main = "Test pooled samples", xlab = "BL1-A", ylab = "BL3-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

cells_AaFnPgPiSmuSsob_pooled <- flowCore::filter(flowData_pooled_AaFnPgPiSmuSsob, polyGate1)
TotalCount_AaFnPgPiSmuSsob_pooled <- summary(cells_AaFnPgPiSmuSsob_pooled)
TotalCount_AaFnPgPiSmuSsob_pooled <- toTable(TotalCount_AaFnPgPiSmuSsob_pooled)

fcs_names_AaFnPgPiSmuSsob_pooled <- c("20200512_Fabian_14strainID_BHI2_24h_Aa",
                                      "20200512_Fabian_14strainID_BHI2_24h_Fn",
                                      "20200512_Fabian_14strainID_BHI2_24h_Pg",
                                      "20200512_Fabian_14strainID_BHI2_24h_Pi",
                                      "20200512_Fabian_14strainID_BHI2_24h_Smu",
                                      "20200512_Fabian_14strainID_BHI2_24h_Ssob")

#Model_RF_AaFnPgPiSmuSsob_pooled <- Phenoflow::RandomF_FCS(flowData_pooled_AaFnPgPiSmuSsob[fcs_names_AaFnPgPiSmuSsob_pooled], Sample_Info_pooled_AaFnPgPiSmuSsob, sample_col = "name_pooled_AaFnPgPiSmuSsob", target_label = "Strain", downsample = 12328, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_AaFnPgPiSmuSsob_pooled, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_AaFnPgPiSmuSsob_pooled.rds")
Model_RF_AaFnPgPiSmuSsob_pooled <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_AaFnPgPiSmuSsobpooled.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_AaFnPgPiSmuSsob <- c("20200512_Fabian_14strainID_BHI2_NA_Mix6_NA_1000_SG_C10.fcs",
                                   "20200512_Fabian_14strainID_BHI2_NA_Mix6_NA_1000_SG_C11.fcs")

flowData_topre_BHI_AaFnPgPiSmuSsob <- flowData_transformed_BHI_gated[fcs_topre_BHI_AaFnPgPiSmuSsob]
test_pred_BHI_AaFnPgPiSmuSsob_pooled <- RandomF_predict(x = Model_RF_AaFnPgPiSmuSsob_pooled[[1]], new_data = flowData_topre_BHI_AaFnPgPiSmuSsob, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_AaFnPgPiSmuSsob_pooled <- left_join(test_pred_BHI_AaFnPgPiSmuSsob_pooled, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_AaFnPgPiSmuSsob_pooled <- left_join(test_pred_BHI2_AaFnPgPiSmuSsob_pooled, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_AaFnPgPiSmuSsob_pooled <- test_pred_BHI2_AaFnPgPiSmuSsob_pooled %>%
  mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsAaFnPgPiSmuSsob_pooled.csv", test_pred_BHI2_AaFnPgPiSmuSsob_pooled)


# Mix7
### Model for Aa, An, Av, Fn, Pg, Pi, Sg, Smi, Smu, So, Ssal, Ssan, Ssob and Vp in BHI2
# Select the fcs files based on which the model will be trained
flowData_pooled_AllStrains <- FCS_pool(x = flowData_transformed_BHI_gated,
                                       stub = c("20200512_Fabian_14strainID_BHI2_24h_Aa",
                                                "20200512_Fabian_14strainID_BHI2_24h_An",
                                                "20200512_Fabian_14strainID_BHI2_24h_Av",
                                                "20200512_Fabian_14strainID_BHI2_24h_Fn",
                                                "20200512_Fabian_14strainID_BHI2_24h_Pg",
                                                "20200512_Fabian_14strainID_BHI2_24h_Pi",
                                                "20200512_Fabian_14strainID_BHI2_24h_Sg",
                                                "20200512_Fabian_14strainID_BHI2_24h_Smi",
                                                "20200512_Fabian_14strainID_BHI2_24h_Smu",
                                                "20200512_Fabian_14strainID_BHI2_24h_So",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssal",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssan",
                                                "20200512_Fabian_14strainID_BHI2_24h_Ssob",
                                                "20200512_Fabian_14strainID_BHI2_24h_Vp"))

metadata_pooled_AllStrains <- data.frame(do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData_pooled_AllStrains), "_"), rbind)))
colnames(metadata_pooled_AllStrains) <- c("Date", "Experimenter", "ExperimentID", "Medium", "Timepoint", "Strain")

name_pooled_AllStrains <- flowCore::sampleNames(flowData_pooled_AllStrains)
Sample_Info_pooled_AllStrains <- cbind(name_pooled_AllStrains, metadata_pooled_AllStrains)

# Plot not working for some reason
# xyplot(`BL3-A`~`BL1-A`, data = flowData_pooled_AllStrains,
#        scales = list(y = list(limits = c(0, 15)),
#                      x = list(limits = c(5, 15))),
#        axis = axis.default, nbin = 125, main = "Test pooled samples", xlab = "BL1-A", ylab = "BL3-A",
#        par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

cells_AllStrains_pooled <- flowCore::filter(flowData_pooled_AllStrains, polyGate1)
TotalCount_AllStrains_pooled <- summary(cells_AllStrains_pooled)
TotalCount_AllStrains_pooled <- toTable(TotalCount_AllStrains_pooled)

fcs_names_AllStrains_pooled <- c("20200512_Fabian_14strainID_BHI2_24h_Aa",
                                 "20200512_Fabian_14strainID_BHI2_24h_An",
                                 "20200512_Fabian_14strainID_BHI2_24h_Av",
                                 "20200512_Fabian_14strainID_BHI2_24h_Fn",
                                 "20200512_Fabian_14strainID_BHI2_24h_Pg",
                                 "20200512_Fabian_14strainID_BHI2_24h_Pi",
                                 "20200512_Fabian_14strainID_BHI2_24h_Sg",
                                 "20200512_Fabian_14strainID_BHI2_24h_Smi",
                                 "20200512_Fabian_14strainID_BHI2_24h_Smu",
                                 "20200512_Fabian_14strainID_BHI2_24h_So",
                                 "20200512_Fabian_14strainID_BHI2_24h_Ssal",
                                 "20200512_Fabian_14strainID_BHI2_24h_Ssan",
                                 "20200512_Fabian_14strainID_BHI2_24h_Ssob",
                                 "20200512_Fabian_14strainID_BHI2_24h_Vp")

#Model_RF_AllStrains_pooled <- Phenoflow::RandomF_FCS(flowData_pooled_AllStrains[fcs_names_AllStrains_pooled], Sample_Info_pooled_AllStrains, sample_col = "name_pooled_AllStrains", target_label = "Strain", downsample = 12328, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)
#saveRDS(object = Model_RF_AllStrains_pooled, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_AllStrains_pooled.rds")
Model_RF_AllStrains_pooled <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RandomForest/RF_BHI_8param_AllStrainspooled.rds")

## Make predictions for relevant mixtures and co-cultures in BHI2
fcs_topre_BHI_AllStrains <- c("20200512_Fabian_14strainID_BHI2_NA_Mix7_NA_1000_SG_D8.fcs",
                              "20200512_Fabian_14strainID_BHI2_NA_Mix7_NA_1000_SG_D9.fcs")

flowData_topre_BHI_AllStrains <- flowData_transformed_BHI_gated[fcs_topre_BHI_AllStrains]
test_pred_BHI_AllStrains_pooled <- RandomF_predict(x = Model_RF_AllStrains_pooled[[1]], new_data = flowData_topre_BHI_AllStrains, cleanFCS = FALSE)

# Export predictions
test_pred_BHI2_AllStrains_pooled <- left_join(test_pred_BHI_AllStrains_pooled, vol2, by = c("Sample" = "Sample_name"))
test_pred_BHI2_AllStrains_pooled <- left_join(test_pred_BHI2_AllStrains_pooled, Sample_Info, by = c("Sample" = "name"))
test_pred_BHI2_AllStrains_pooled <- test_pred_BHI2_AllStrains_pooled %>%
  mutate(Concentration = (Counts*Dilution)/Volume)
# Export as csv file
write.csv2(file = "PredictedCellsAllStrains_pooled.csv", test_pred_BHI2_AllStrains_pooled)


## 8.3. Analysis of performance ----

### 8.3.1. Number of strains vs accuracy ----
# Extract accuracy for each model
modellist <- c("SoFn",
               "SoFnPg",
               "SoFnPgVp",
               "AnAvSgSmiSoSsalSsanVp",
               "SgSmiSmuSoSsalSsanSsob",
               "AaFnPgPiSmuSsob",
               "AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp")

numberstrains <- sapply(modellist, count_capital_letters)

acc_RF <- c(Model_RF_SoFn_pooled[[2]][["overall"]][["Accuracy"]],
            Model_RF_SoFnPg_pooled[[2]][["overall"]][["Accuracy"]],
            Model_RF_SoFnPgVp_pooled[[2]][["overall"]][["Accuracy"]],
            Model_RF_AnAvSgSmiSoSsalSsanVp_pooled[[2]][["overall"]][["Accuracy"]],
            Model_RF_SgSmiSmuSoSsalSsanSsob_pooled[[2]][["overall"]][["Accuracy"]],
            Model_RF_AaFnPgPiSmuSsob_pooled[[2]][["overall"]][["Accuracy"]],
            Model_RF_AllStrains_pooled[[2]][["overall"]][["Accuracy"]])

accuracy_RF <- data.frame(Model = modellist,
                          NumberStrains = numberstrains,
                          Accuracy = acc_RF*100,
                          Method = "RandomForest",
                          row.names = NULL)

random_guessing <- data.frame(NumberStrains = c(2:14),
                              Accuracy = (c(1/(2:14)))*100,
                              Method = "RandomGuessing")

accuracy_RF2 <- dplyr::bind_rows(accuracy_RF, random_guessing)
accuracy_RF2 <- accuracy_RF2[, c(2:4)]

plot_accuracy_RF <- ggplot(data = accuracy_RF2, aes(x = NumberStrains, y = Accuracy, color = Method))+
  geom_point(data = subset(accuracy_RF2, Method == "RandomForest"), size = 7) +
  geom_line(data = subset(accuracy_RF2, Method == "RandomGuessing"), size = 2) +
  labs(x = "Number of strains", y = "Accuracy (%)", color = NULL) +
  scale_color_manual(values = c("RandomForest"= "black", "RandomGuessing"= "darkred"),
                     labels = c("RandomForest"= "Random Forest", "RandomGuessing" = "Random Guessing")) +
  paper_theme_fab +
  scale_y_continuous(limits = c(0, 100)) +
  scale_x_continuous(limits = c(0, 14), breaks = seq(0, 14, by = 1))
print(plot_accuracy_RF)


### 8.3.2. Performance of predictions ----
# RMSE (root mean squared error) for predictions

# Theoretical abundance
theoretical_mocks_rect <- reshape2::dcast(theoretical_mocks, Mix_ID ~ Strain)
theoretical_mocks_rect[is.na(theoretical_mocks_rect)] <- 0

theoretical_mocks_prop <- theoretical_mocks_rect
theoretical_mocks_prop[, -c(1)] <- sweep(as.matrix(theoretical_mocks_rect[, -c(1)]), 1, rowSums(theoretical_mocks_rect[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
colnames(theoretical_mocks_prop)[colnames(theoretical_mocks_prop) == "Mix_ID"] <- "ID"

theoretical_longer <- tidyr::pivot_longer(theoretical_mocks_prop, cols = -ID, names_to = "Species", values_to = "Actual")


# Models trained on 50000 events
FCM_SoFn <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFn_pooled_50kevents.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_SoFnPg <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFnPg_pooled_50kevents.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_SoFnPgVp <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFnPgVp_pooled_50kevents.csv", header = T, sep = ";", stringsAsFactors = F)

# Change , to . in FCM dataframes
FCM_SoFn$Concentration <- gsub(",", ".", FCM_SoFn$Concentration)
FCM_SoFnPg$Concentration <- gsub(",", ".", FCM_SoFnPg$Concentration)
FCM_SoFnPgVp$Concentration <- gsub(",", ".", FCM_SoFnPgVp$Concentration)

FCM_SoFn$Concentration <- as.numeric(FCM_SoFn$Concentration)
FCM_SoFnPg$Concentration <- as.numeric(FCM_SoFnPg$Concentration)
FCM_SoFnPgVp$Concentration <- as.numeric(FCM_SoFnPgVp$Concentration)

FCM_SoFn_mean <- FCM_SoFn
FCM_SoFn_mean$Replicate[is.na(FCM_SoFn_mean$Replicate)] <- "Z"
FCM_SoFn_mean$Timepoint[is.na(FCM_SoFn_mean$Timepoint)] <- "0h"
FCM_SoFn_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFn_mean, FUN = mean)
FCM_SoFn_mean <- subset(FCM_SoFn_mean, Replicate == "Z")
colnames(FCM_SoFn_mean)[colnames(FCM_SoFn_mean) == "Strain"] <- "ID"
FCM_SoFn_mean <- FCM_SoFn_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFn_mean_rect <- reshape2::dcast(FCM_SoFn_mean, ID ~ Predicted_label)
FCM_SoFn_mean_prop <- FCM_SoFn_mean_rect
FCM_SoFn_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_SoFn_mean_prop[, -c(1)]), 1, rowSums(FCM_SoFn_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
FCM_SoFn_mean_prop$Model <- "SoFn"

FCM_SoFnPg_mean <- FCM_SoFnPg
FCM_SoFnPg_mean$Replicate[is.na(FCM_SoFnPg_mean$Replicate)] <- "Z"
FCM_SoFnPg_mean$Timepoint[is.na(FCM_SoFnPg_mean$Timepoint)] <- "0h"
FCM_SoFnPg_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFnPg_mean, FUN = mean)
FCM_SoFnPg_mean <- subset(FCM_SoFnPg_mean, Replicate == "Z")
colnames(FCM_SoFnPg_mean)[colnames(FCM_SoFnPg_mean) == "Strain"] <- "ID"
FCM_SoFnPg_mean <- FCM_SoFnPg_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFnPg_mean_rect <- reshape2::dcast(FCM_SoFnPg_mean, ID ~ Predicted_label)
FCM_SoFnPg_mean_prop <- FCM_SoFnPg_mean_rect
FCM_SoFnPg_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_SoFnPg_mean_prop[, -c(1)]), 1, rowSums(FCM_SoFnPg_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
FCM_SoFnPg_mean_prop$Model <- "SoFnPg"

FCM_SoFnPgVp_mean <- FCM_SoFnPgVp
FCM_SoFnPgVp_mean$Replicate[is.na(FCM_SoFnPgVp_mean$Replicate)] <- "Z"
FCM_SoFnPgVp_mean$Timepoint[is.na(FCM_SoFnPgVp_mean$Timepoint)] <- "0h"
FCM_SoFnPgVp_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFnPgVp_mean, FUN = mean)
FCM_SoFnPgVp_mean <- subset(FCM_SoFnPgVp_mean, Replicate == "Z")
colnames(FCM_SoFnPgVp_mean)[colnames(FCM_SoFnPgVp_mean) == "Strain"] <- "ID"
FCM_SoFnPgVp_mean <- FCM_SoFnPgVp_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFnPgVp_mean_rect <- reshape2::dcast(FCM_SoFnPgVp_mean, ID ~ Predicted_label)
FCM_SoFnPgVp_mean_prop <- FCM_SoFnPgVp_mean_rect
FCM_SoFnPgVp_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_SoFnPgVp_mean_prop[, -c(1)]), 1, rowSums(FCM_SoFnPgVp_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
FCM_SoFnPgVp_mean_prop$Model <- "SoFnPgVp"

FCM_mean_prop_50k <- dplyr::bind_rows(FCM_SoFn_mean_prop, FCM_SoFnPg_mean_prop, FCM_SoFnPgVp_mean_prop)
FCM_mean_prop_50k[is.na(FCM_mean_prop_50k)] <- 0

FCM_SoFn_mean_prop_50k <- subset(FCM_mean_prop_50k, Model == "SoFn")
FCM_SoFn_mean_prop_50k <- FCM_SoFn_mean_prop_50k[, c("ID", "So", "Fn", "Pg", "Vp")]
FCM_SoFnPg_mean_prop_50k <- subset(FCM_mean_prop_50k, Model == "SoFnPg")
FCM_SoFnPg_mean_prop_50k <- FCM_SoFnPg_mean_prop_50k[, c("ID", "So", "Fn", "Pg", "Vp")]
FCM_SoFnPgVp_mean_prop_50k <- subset(FCM_mean_prop_50k, Model == "SoFnPgVp")
FCM_SoFnPgVp_mean_prop_50k <- FCM_SoFnPgVp_mean_prop_50k[, c("ID", "So", "Fn", "Pg", "Vp")]

# Calculation RSME
FCM_SoFn_50k_longer <- tidyr::pivot_longer(FCM_SoFn_mean_prop_50k, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_SoFn_50k <- merge(FCM_SoFn_50k_longer, theoretical_longer, by = c("ID", "Species"))
IDs_SoFn <- unique(merged_SoFn_50k$ID)
RMSE_SoFn_50k <- data.frame(ID = IDs_SoFn,
                        RMSE_SoFn = numeric(length(IDs_SoFn)))

for (i in seq_along(IDs_SoFn)) {
  ID <- IDs_SoFn[i]
  subset_df <- merged_SoFn_50k[merged_SoFn_50k$ID == ID, ]
  RMSE_SoFn_50k$RMSE_SoFn[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}


FCM_SoFnPg_50k_longer <- tidyr::pivot_longer(FCM_SoFnPg_mean_prop_50k, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_SoFnPg_50k <- merge(FCM_SoFnPg_50k_longer, theoretical_longer, by = c("ID", "Species"))
IDs_SoFnPg <- unique(merged_SoFnPg_50k$ID)
RMSE_SoFnPg_50k <- data.frame(ID = IDs_SoFnPg,
                          RMSE_SoFnPg = numeric(length(IDs_SoFnPg)))

for (i in seq_along(IDs_SoFnPg)) {
  ID <- IDs_SoFnPg[i]
  subset_df <- merged_SoFnPg_50k[merged_SoFnPg_50k$ID == ID, ]
  RMSE_SoFnPg_50k$RMSE_SoFnPg[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}

FCM_SoFnPgVp_50k_longer <- tidyr::pivot_longer(FCM_SoFnPgVp_mean_prop_50k, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_SoFnPgVp_50k <- merge(FCM_SoFnPgVp_50k_longer, theoretical_longer, by = c("ID", "Species"))
IDs_SoFnPgVp <- unique(merged_SoFnPgVp_50k$ID)
RMSE_SoFnPgVp_50k <- data.frame(ID = IDs_SoFnPgVp,
                            RMSE_SoFnPgVp = numeric(length(IDs_SoFnPgVp)))

for (i in seq_along(IDs_SoFnPgVp)) {
  ID <- IDs_SoFnPgVp[i]
  subset_df <- merged_SoFnPgVp_50k[merged_SoFnPgVp_50k$ID == ID, ]
  RMSE_SoFnPgVp_50k$RMSE_SoFnPgVp[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}

RMSE_50k <- merge(RMSE_SoFn_50k, RMSE_SoFnPg_50k, by = "ID", all = TRUE)
RMSE_50k <- merge(RMSE_50k, RMSE_SoFnPgVp_50k, by = "ID", all = TRUE)

#saveRDS(object = RMSE_50k, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/RMSE_FCM_50k.rds")


# Models trained on 12328 events
FCM_SoFn <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFn_pooled.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_SoFnPg <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFnPg_pooled.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_SoFnPgVp <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFnPgVp_pooled.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_AnAvSgSmiSoSsalSsanVp <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsAnAvSgSmiSoSsalSsanVp_pooled.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_SgSmiSmuSoSsalSsanSsob <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSgSmiSmuSoSsalSsanSsob_pooled.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_AaFnPgPiSmuSsob <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsAaFnPgPiSmuSsob_pooled.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsAllStrains_pooled.csv", header = T, sep = ";", stringsAsFactors = F)

# Change , to . in FCM dataframes
FCM_SoFn$Concentration <- gsub(",", ".", FCM_SoFn$Concentration)
FCM_SoFnPg$Concentration <- gsub(",", ".", FCM_SoFnPg$Concentration)
FCM_SoFnPgVp$Concentration <- gsub(",", ".", FCM_SoFnPgVp$Concentration)
FCM_AnAvSgSmiSoSsalSsanVp$Concentration <- gsub(",", ".", FCM_AnAvSgSmiSoSsalSsanVp$Concentration)
FCM_SgSmiSmuSoSsalSsanSsob$Concentration <- gsub(",", ".", FCM_SgSmiSmuSoSsalSsanSsob$Concentration)
FCM_AaFnPgPiSmuSsob$Concentration <- gsub(",", ".", FCM_AaFnPgPiSmuSsob$Concentration)
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp$Concentration <- gsub(",", ".", FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp$Concentration)

FCM_SoFn$Concentration <- as.numeric(FCM_SoFn$Concentration)
FCM_SoFnPg$Concentration <- as.numeric(FCM_SoFnPg$Concentration)
FCM_SoFnPgVp$Concentration <- as.numeric(FCM_SoFnPgVp$Concentration)
FCM_AnAvSgSmiSoSsalSsanVp$Concentration <- as.numeric(FCM_AnAvSgSmiSoSsalSsanVp$Concentration)
FCM_SgSmiSmuSoSsalSsanSsob$Concentration <- as.numeric(FCM_SgSmiSmuSoSsalSsanSsob$Concentration)
FCM_AaFnPgPiSmuSsob$Concentration <- as.numeric(FCM_AaFnPgPiSmuSsob$Concentration)
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp$Concentration <- as.numeric(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp$Concentration)

FCM_SoFn_mean <- FCM_SoFn
FCM_SoFn_mean$Replicate[is.na(FCM_SoFn_mean$Replicate)] <- "Z"
FCM_SoFn_mean$Timepoint[is.na(FCM_SoFn_mean$Timepoint)] <- "0h"
FCM_SoFn_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFn_mean, FUN = mean)
FCM_SoFn_mean <- subset(FCM_SoFn_mean, Replicate == "Z")
colnames(FCM_SoFn_mean)[colnames(FCM_SoFn_mean) == "Strain"] <- "ID"
FCM_SoFn_mean <- FCM_SoFn_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFn_mean_rect <- reshape2::dcast(FCM_SoFn_mean, ID ~ Predicted_label)
FCM_SoFn_mean_prop <- FCM_SoFn_mean_rect
FCM_SoFn_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_SoFn_mean_prop[, -c(1)]), 1, rowSums(FCM_SoFn_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
#FCM_SoFn_mean_prop$Model <- "SoFn"

FCM_SoFnPg_mean <- FCM_SoFnPg
FCM_SoFnPg_mean$Replicate[is.na(FCM_SoFnPg_mean$Replicate)] <- "Z"
FCM_SoFnPg_mean$Timepoint[is.na(FCM_SoFnPg_mean$Timepoint)] <- "0h"
FCM_SoFnPg_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFnPg_mean, FUN = mean)
FCM_SoFnPg_mean <- subset(FCM_SoFnPg_mean, Replicate == "Z")
colnames(FCM_SoFnPg_mean)[colnames(FCM_SoFnPg_mean) == "Strain"] <- "ID"
FCM_SoFnPg_mean <- FCM_SoFnPg_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFnPg_mean_rect <- reshape2::dcast(FCM_SoFnPg_mean, ID ~ Predicted_label)
FCM_SoFnPg_mean_prop <- FCM_SoFnPg_mean_rect
FCM_SoFnPg_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_SoFnPg_mean_prop[, -c(1)]), 1, rowSums(FCM_SoFnPg_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
#FCM_SoFnPg_mean_prop$Model <- "SoFnPg"

FCM_SoFnPgVp_mean <- FCM_SoFnPgVp
FCM_SoFnPgVp_mean$Replicate[is.na(FCM_SoFnPgVp_mean$Replicate)] <- "Z"
FCM_SoFnPgVp_mean$Timepoint[is.na(FCM_SoFnPgVp_mean$Timepoint)] <- "0h"
FCM_SoFnPgVp_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFnPgVp_mean, FUN = mean)
FCM_SoFnPgVp_mean <- subset(FCM_SoFnPgVp_mean, Replicate == "Z")
colnames(FCM_SoFnPgVp_mean)[colnames(FCM_SoFnPgVp_mean) == "Strain"] <- "ID"
FCM_SoFnPgVp_mean <- FCM_SoFnPgVp_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFnPgVp_mean_rect <- reshape2::dcast(FCM_SoFnPgVp_mean, ID ~ Predicted_label)
FCM_SoFnPgVp_mean_prop <- FCM_SoFnPgVp_mean_rect
FCM_SoFnPgVp_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_SoFnPgVp_mean_prop[, -c(1)]), 1, rowSums(FCM_SoFnPgVp_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
#FCM_SoFnPgVp_mean_prop$Model <- "SoFnPgVp"

FCM_AnAvSgSmiSoSsalSsanVp_mean <- FCM_AnAvSgSmiSoSsalSsanVp
FCM_AnAvSgSmiSoSsalSsanVp_mean$Replicate[is.na(FCM_AnAvSgSmiSoSsalSsanVp_mean$Replicate)] <- "Z"
FCM_AnAvSgSmiSoSsalSsanVp_mean$Timepoint[is.na(FCM_AnAvSgSmiSoSsalSsanVp_mean$Timepoint)] <- "0h"
FCM_AnAvSgSmiSoSsalSsanVp_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_AnAvSgSmiSoSsalSsanVp_mean, FUN = mean)
FCM_AnAvSgSmiSoSsalSsanVp_mean <- subset(FCM_AnAvSgSmiSoSsalSsanVp_mean, Replicate == "Z")
colnames(FCM_AnAvSgSmiSoSsalSsanVp_mean)[colnames(FCM_AnAvSgSmiSoSsalSsanVp_mean) == "Strain"] <- "ID"
FCM_AnAvSgSmiSoSsalSsanVp_mean <- FCM_AnAvSgSmiSoSsalSsanVp_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_AnAvSgSmiSoSsalSsanVp_mean_rect <- reshape2::dcast(FCM_AnAvSgSmiSoSsalSsanVp_mean, ID ~ Predicted_label)
FCM_AnAvSgSmiSoSsalSsanVp_mean_prop <- FCM_AnAvSgSmiSoSsalSsanVp_mean_rect
FCM_AnAvSgSmiSoSsalSsanVp_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_AnAvSgSmiSoSsalSsanVp_mean_prop[, -c(1)]), 1, rowSums(FCM_AnAvSgSmiSoSsalSsanVp_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
#FCM_AnAvSgSmiSoSsalSsanVp_mean_prop$Model <- "AnAvSgSmiSoSsalSsanVp"

FCM_SgSmiSmuSoSsalSsanSsob_mean <- FCM_SgSmiSmuSoSsalSsanSsob
FCM_SgSmiSmuSoSsalSsanSsob_mean$Replicate[is.na(FCM_SgSmiSmuSoSsalSsanSsob_mean$Replicate)] <- "Z"
FCM_SgSmiSmuSoSsalSsanSsob_mean$Timepoint[is.na(FCM_SgSmiSmuSoSsalSsanSsob_mean$Timepoint)] <- "0h"
FCM_SgSmiSmuSoSsalSsanSsob_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SgSmiSmuSoSsalSsanSsob_mean, FUN = mean)
FCM_SgSmiSmuSoSsalSsanSsob_mean <- subset(FCM_SgSmiSmuSoSsalSsanSsob_mean, Replicate == "Z")
colnames(FCM_SgSmiSmuSoSsalSsanSsob_mean)[colnames(FCM_SgSmiSmuSoSsalSsanSsob_mean) == "Strain"] <- "ID"
FCM_SgSmiSmuSoSsalSsanSsob_mean <- FCM_SgSmiSmuSoSsalSsanSsob_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SgSmiSmuSoSsalSsanSsob_mean_rect <- reshape2::dcast(FCM_SgSmiSmuSoSsalSsanSsob_mean, ID ~ Predicted_label)
FCM_SgSmiSmuSoSsalSsanSsob_mean_prop <- FCM_SgSmiSmuSoSsalSsanSsob_mean_rect
FCM_SgSmiSmuSoSsalSsanSsob_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_SgSmiSmuSoSsalSsanSsob_mean_prop[, -c(1)]), 1, rowSums(FCM_SgSmiSmuSoSsalSsanSsob_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
#FCM_SgSmiSmuSoSsalSsanSsob_mean_prop$Model <- "SgSmiSmuSoSsalSsanSsob"

FCM_AaFnPgPiSmuSsob_mean <- FCM_AaFnPgPiSmuSsob
FCM_AaFnPgPiSmuSsob_mean$Replicate[is.na(FCM_AaFnPgPiSmuSsob_mean$Replicate)] <- "Z"
FCM_AaFnPgPiSmuSsob_mean$Timepoint[is.na(FCM_AaFnPgPiSmuSsob_mean$Timepoint)] <- "0h"
FCM_AaFnPgPiSmuSsob_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_AaFnPgPiSmuSsob_mean, FUN = mean)
FCM_AaFnPgPiSmuSsob_mean <- subset(FCM_AaFnPgPiSmuSsob_mean, Replicate == "Z")
colnames(FCM_AaFnPgPiSmuSsob_mean)[colnames(FCM_AaFnPgPiSmuSsob_mean) == "Strain"] <- "ID"
FCM_AaFnPgPiSmuSsob_mean <- FCM_AaFnPgPiSmuSsob_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_AaFnPgPiSmuSsob_mean_rect <- reshape2::dcast(FCM_AaFnPgPiSmuSsob_mean, ID ~ Predicted_label)
FCM_AaFnPgPiSmuSsob_mean_prop <- FCM_AaFnPgPiSmuSsob_mean_rect
FCM_AaFnPgPiSmuSsob_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_AaFnPgPiSmuSsob_mean_prop[, -c(1)]), 1, rowSums(FCM_AaFnPgPiSmuSsob_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
#FCM_AaFnPgPiSmuSsob_mean_prop$Model <- "AaFnPgPiSmuSsob"

FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean <- FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean$Replicate[is.na(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean$Replicate)] <- "Z"
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean$Timepoint[is.na(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean$Timepoint)] <- "0h"
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean, FUN = mean)
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean <- subset(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean, Replicate == "Z")
colnames(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean)[colnames(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean) == "Strain"] <- "ID"
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean <- FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_rect <- reshape2::dcast(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean, ID ~ Predicted_label)
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_prop <- FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_rect
FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_prop[, -c(1)]), 1, rowSums(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
#FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_prop$Model <- "AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp"

# FCM_mean_prop <- dplyr::bind_rows(FCM_SoFn_mean_prop,
#                                   FCM_SoFnPg_mean_prop,
#                                   FCM_SoFnPgVp_mean_prop,
#                                   FCM_AnAvSgSmiSoSsalSsanVp_mean_prop,
#                                   FCM_SgSmiSmuSoSsalSsanSsob_mean_prop,
#                                   FCM_AaFnPgPiSmuSsob_mean_prop,
#                                   FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_prop)
# FCM_mean_prop[is.na(FCM_mean_prop)] <- 0

#saveRDS(object = FCM_mean_prop, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/FCM_mean_prop_12328cells.rds")

# Calculation RSME
FCM_SoFn_longer <- tidyr::pivot_longer(FCM_SoFn_mean_prop, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_SoFn <- merge(FCM_SoFn_longer, theoretical_longer, by = c("ID", "Species"))
IDs_SoFn <- unique(merged_SoFn$ID)
RMSE_SoFn <- data.frame(ID = IDs_SoFn,
                        RMSE_SoFn = numeric(length(IDs_SoFn)))

for (i in seq_along(IDs_SoFn)) {
  ID <- IDs_SoFn[i]
  subset_df <- merged_SoFn[merged_SoFn$ID == ID, ]
  RMSE_SoFn$RMSE_SoFn[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}


FCM_SoFnPg_longer <- tidyr::pivot_longer(FCM_SoFnPg_mean_prop, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_SoFnPg <- merge(FCM_SoFnPg_longer, theoretical_longer, by = c("ID", "Species"))
IDs_SoFnPg <- unique(merged_SoFnPg$ID)
RMSE_SoFnPg <- data.frame(ID = IDs_SoFnPg,
                          RMSE_SoFnPg = numeric(length(IDs_SoFnPg)))

for (i in seq_along(IDs_SoFnPg)) {
  ID <- IDs_SoFnPg[i]
  subset_df <- merged_SoFnPg[merged_SoFnPg$ID == ID, ]
  RMSE_SoFnPg$RMSE_SoFnPg[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}


FCM_SoFnPgVp_longer <- tidyr::pivot_longer(FCM_SoFnPgVp_mean_prop, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_SoFnPgVp <- merge(FCM_SoFnPgVp_longer, theoretical_longer, by = c("ID", "Species"))
IDs_SoFnPgVp <- unique(merged_SoFnPgVp$ID)
RMSE_SoFnPgVp <- data.frame(ID = IDs_SoFnPgVp,
                            RMSE_SoFnPgVp = numeric(length(IDs_SoFnPgVp)))

for (i in seq_along(IDs_SoFnPgVp)) {
  ID <- IDs_SoFnPgVp[i]
  subset_df <- merged_SoFnPgVp[merged_SoFnPgVp$ID == ID, ]
  RMSE_SoFnPgVp$RMSE_SoFnPgVp[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}


FCM_AnAvSgSmiSoSsalSsanVp_longer <- tidyr::pivot_longer(FCM_AnAvSgSmiSoSsalSsanVp_mean_prop, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_AnAvSgSmiSoSsalSsanVp <- merge(FCM_AnAvSgSmiSoSsalSsanVp_longer, theoretical_longer, by = c("ID", "Species"))
IDs_AnAvSgSmiSoSsalSsanVp <- unique(merged_AnAvSgSmiSoSsalSsanVp$ID)
RMSE_AnAvSgSmiSoSsalSsanVp <- data.frame(ID = IDs_AnAvSgSmiSoSsalSsanVp,
                                         RMSE_AnAvSgSmiSoSsalSsanVp = numeric(length(IDs_AnAvSgSmiSoSsalSsanVp)))

for (i in seq_along(IDs_AnAvSgSmiSoSsalSsanVp)) {
  ID <- IDs_AnAvSgSmiSoSsalSsanVp[i]
  subset_df <- merged_AnAvSgSmiSoSsalSsanVp[merged_AnAvSgSmiSoSsalSsanVp$ID == ID, ]
  RMSE_AnAvSgSmiSoSsalSsanVp$RMSE_AnAvSgSmiSoSsalSsanVp[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}


FCM_SgSmiSmuSoSsalSsanSsob_longer <- tidyr::pivot_longer(FCM_SgSmiSmuSoSsalSsanSsob_mean_prop, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_SgSmiSmuSoSsalSsanSsob <- merge(FCM_SgSmiSmuSoSsalSsanSsob_longer, theoretical_longer, by = c("ID", "Species"))
IDs_SgSmiSmuSoSsalSsanSsob <- unique(merged_SgSmiSmuSoSsalSsanSsob$ID)
RMSE_SgSmiSmuSoSsalSsanSsob <- data.frame(ID = IDs_SgSmiSmuSoSsalSsanSsob,
                                          RMSE_SgSmiSmuSoSsalSsanSsob = numeric(length(IDs_SgSmiSmuSoSsalSsanSsob)))

for (i in seq_along(IDs_SgSmiSmuSoSsalSsanSsob)) {
  ID <- IDs_SgSmiSmuSoSsalSsanSsob[i]
  subset_df <- merged_SgSmiSmuSoSsalSsanSsob[merged_SgSmiSmuSoSsalSsanSsob$ID == ID, ]
  RMSE_SgSmiSmuSoSsalSsanSsob$RMSE_SgSmiSmuSoSsalSsanSsob[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}


FCM_AaFnPgPiSmuSsob_longer <- tidyr::pivot_longer(FCM_AaFnPgPiSmuSsob_mean_prop, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_AaFnPgPiSmuSsob <- merge(FCM_AaFnPgPiSmuSsob_longer, theoretical_longer, by = c("ID", "Species"))
IDs_AaFnPgPiSmuSsob <- unique(merged_AaFnPgPiSmuSsob$ID)
RMSE_AaFnPgPiSmuSsob <- data.frame(ID = IDs_AaFnPgPiSmuSsob,
                                   RMSE_AaFnPgPiSmuSsob = numeric(length(IDs_AaFnPgPiSmuSsob)))

for (i in seq_along(IDs_AaFnPgPiSmuSsob)) {
  ID <- IDs_AaFnPgPiSmuSsob[i]
  subset_df <- merged_AaFnPgPiSmuSsob[merged_AaFnPgPiSmuSsob$ID == ID, ]
  RMSE_AaFnPgPiSmuSsob$RMSE_AaFnPgPiSmuSsob[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}


FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_longer <- tidyr::pivot_longer(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_mean_prop, cols = -ID, names_to = "Species", values_to = "Predicted")
merged_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp <- merge(FCM_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp_longer, theoretical_longer, by = c("ID", "Species"))
IDs_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp <- unique(merged_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp$ID)
RMSE_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp <- data.frame(ID = IDs_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp,
                                                        RMSE_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp = numeric(length(IDs_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp)))

for (i in seq_along(IDs_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp)) {
  ID <- IDs_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp[i]
  subset_df <- merged_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp[merged_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp$ID == ID, ]
  RMSE_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp$RMSE_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp[i] <- sqrt(mean((subset_df$Predicted - subset_df$Actual)^2))
}


RMSE <- merge(RMSE_SoFn, RMSE_SoFnPg, by = "ID", all = TRUE)
RMSE <- merge(RMSE, RMSE_SoFnPgVp, by = "ID", all = TRUE)
RMSE <- merge(RMSE, RMSE_AnAvSgSmiSoSsalSsanVp, by = "ID", all = TRUE)
RMSE <- merge(RMSE, RMSE_SgSmiSmuSoSsalSsanSsob, by = "ID", all = TRUE)
RMSE <- merge(RMSE, RMSE_AaFnPgPiSmuSsob, by = "ID", all = TRUE)
RMSE <- merge(RMSE, RMSE_AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp, by = "ID", all = TRUE)

#saveRDS(object = RMSE, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/RMSE_FCM.rds")





