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
library("ggcyto")

set.seed(777)


#### Loading data ----

Datapath <- "/media/projects1/Fabian/Oral microbiome/StrainRecognitionFC/AttuneNew"
fcsfiles <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = TRUE)
flowData <- flowCore::read.flowSet(files = fcsfiles, transformation = FALSE, emptyValue = F)
metadata <- read.table("/media/projects1/Fabian/Oral microbiome/StrainRecognitionFC/AttuneNew/20201125_MetaData_AttuneNew.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)


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

## Subsetting data in IS Fabian and IS FM
flowData_transformed_Fab <- flowData_transformed[c(1:19)]
flowData_transformed_FM <- flowData_transformed[c(20:39)]

metadata_Fab <- metadata[c(1:19),]
metadata_FM <- metadata[c(20:39),]


#### Exploration of data and gating ----

## Visualization of different parameters
xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed[c(1:39)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "BL1-BL3", xlab = "BL1-A", ylab = "BL3-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL1-A`~`FSC-A`, data = flowData_transformed[c(1:39)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "FSC-BL1", xlab = "FSC-A", ylab = "BL1-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

xyplot(`BL1-A`~`SSC-A`, data = flowData_transformed[c(1:39)],
       scales = list(y = list(limits = c(0, 16)),
                     x = list(limits = c(0, 16))),
       axis = axis.default, nbin = 125, main = "SSC-BL1", xlab = "SSC-A", ylab = "BL1-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)

## Singlets
p_singlets <- xyplot(`SSC-W`~`SSC-H`, data = flowData_transformed[c(1:39)],
                     scales = list(y = list(limits = c(0, 16)),
                                   x = list(limits = c(0, 16))),
                     axis = axis.default, nbin = 125, main = "Singlets SSC (width vs height)", xlab = "SSC-H", ylab = "SSC-W",
                     par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_singlets)

p_density_SSC_W <- autoplot(flowData_transformed[c(1:39)], "SSC-W")+
  labs(title = "Density plot SSC-W",
       x = "SSC-W",
       y = "Density")
print(p_density_SSC_W)

## Gating based on BL1-BL3
# IS Fabian
sqrcut1 <- matrix(c(6.7, 13.3, 14.5, 10, 6.7,
                    2, 5, 13.5, 10.5, 5.3), ncol = 2, nrow = 5)
colnames(sqrcut1) <- c("BL1-A", "BL3-A")
polyGate1 <- polygonGate(.gate = sqrcut1, filterId = "Cells")

p_gating_Fab <- xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_Fab, filter = polyGate1,
                       scales = list(y = list(limits = c(0, 16)),
                                     x = list(limits = c(0, 16))),
                       axis = axis.default, nbin = 125, main = "Gating cells IS Fab (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
                       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_gating_Fab)

# IS FM
sqrcut2 <- matrix(c(9.3, 14.5, 15.5, 10.8, 9.3,
                    4, 6.5, 13.5, 9.7, 7.2), ncol = 2, nrow = 5)
colnames(sqrcut2) <- c("BL1-A", "BL3-A")
polyGate2 <- polygonGate(.gate = sqrcut2, filterId = "Cells")

p_gating_FM <- xyplot(`BL3-A`~`BL1-A`, data = flowData_transformed_FM, filter = polyGate2,
                      scales = list(y = list(limits = c(0, 16)),
                                    x = list(limits = c(0, 16))),
                      axis = axis.default, nbin = 125, main = "Gating cells IS FM (BL1-BL3)", xlab = "BL1-A", ylab = "BL3-A",
                      par.strip.text = list(col = "white", font = 1, cex = 1), smooth = F)
print(p_gating_FM)

## Subsetting cells
flowData_transformed_Fab_gated <- Subset(flowData_transformed_Fab, polyGate1)
flowData_transformed_FM_gated <- Subset(flowData_transformed_FM, polyGate2)

## Singlet analysis after gating
p_density_SSC_W_Fab <- autoplot(flowData_transformed_Fab_gated, "SSC-W")+
  labs(title = "Density plot Fab SSC-W",
       x = "SSC-W",
       y = "Density")
print(p_density_SSC_W_Fab)

p_density_SSC_W_FM <- autoplot(flowData_transformed_FM_gated, "SSC-W")+
  labs(title = "Density plot FM SSC-W",
       x = "SSC-W",
       y = "Density")
print(p_density_SSC_W_FM)


#### Cell concentrations ----

## Cell counts
cells_Fab <- flowCore::filter(flowData_transformed_Fab, polyGate1)
TotalCount_Fab <- summary(cells_Fab)
TotalCount_Fab <- toTable(TotalCount_Fab)

cells_FM <- flowCore::filter(flowData_transformed_FM, polyGate2)
TotalCount_FM <- summary(cells_FM)
TotalCount_FM <- toTable(TotalCount_FM)

## Extracting volumes
# Volumes are in µL
vol_Fab <- as.numeric(flowCore::fsApply(flowData_transformed_Fab, FUN = function(x) x@description$`$VOL`))/1000
vol_FM <- as.numeric(flowCore::fsApply(flowData_transformed_FM, FUN = function(x) x@description$`$VOL`))/1000

## Concentrations
# Calculating concentrations --> Concentrations will be in cells/µL
cell_concentrations_Fab <- data.frame(Sample_name = flowCore::sampleNames(flowData_transformed_Fab),
                                      Strain = metadata_Fab$Strain,
                                      Concentration = (TotalCount_Fab$true*metadata_Fab$Dilution)/vol_Fab)

cell_concentrations_FM <- data.frame(Sample_name = flowCore::sampleNames(flowData_transformed_FM),
                                     Strain = metadata_FM$Strain,
                                     Concentration = (TotalCount_FM$true*metadata_FM$Dilution)/vol_FM)


#### Training Random Forest classifier ----

## Define parameters on which RF will build its model and extract Sample_Info
paramRF = c("FSC-A", "SSC-A", "BL1-A", "BL3-A", "FSC-H", "SSC-H", "BL1-H", "BL3-H", "FSC-W", "SSC-W", "BL1-W", "BL3-W")
name <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = FALSE)
Sample_Info <- cbind(name, metadata)

### IS Fabian
## Select the fcs files based on which the model will be trained --> Streps from Ioanna's stock
fcs_names_streps_Fab <- c("20201125_Fabian_OralOcular_ISFabian_Experiment_OralIoanna_1000_SG_A1.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralIoanna_1000_SG_B1.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralIoanna_1000_SG_C1.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralIoanna_1000_SG_D1.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralIoanna_1000_SG_E1.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralIoanna_1000_SG_F1.fcs")

Sample_Info_streps_Fab <- Sample_Info %>% dplyr::filter(name %in% fcs_names_streps_Fab)

Model_RF_streps_Fab <- Phenoflow::RandomF_FCS(flowData_transformed_Fab_gated[fcs_names_streps_Fab], Sample_Info_streps_Fab, target_label = "Strain", downsample = 22000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for Streps from Fabian's stock
fcs_topre_streps_Fab <- c("20201125_Fabian_OralOcular_ISFabian_Experiment_OralFabian_1000_SG_A2.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralFabian_1000_SG_B2.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralFabian_1000_SG_C2.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralFabian_1000_SG_D2.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralFabian_1000_SG_E2.fcs",
                          "20201125_Fabian_OralOcular_ISFabian_Experiment_OralFabian_1000_SG_F2.fcs")

flowData_topre_streps_Fab <- flowData_transformed_Fab_gated[fcs_topre_streps_Fab]
test_pred_streps_Fab <- RandomF_predict(x = Model_RF_streps_Fab[[1]], new_data =  flowData_topre_streps_Fab, cleanFCS = FALSE)

## Calculate cell concentrations of predictions
vol_Fab2 <- data.frame(Sample_name = flowCore::sampleNames(flowData_transformed_Fab),
                       Volume = vol_Fab)

test_pred_streps_Fab <- left_join(test_pred_streps_Fab, vol_Fab2, by = c("Sample" = "Sample_name"))
test_pred_streps_Fab <- left_join(test_pred_streps_Fab, Sample_Info_streps_Fab, by = c("Sample" = "File_Name"))
test_pred_streps_Fab <- test_pred_streps_Fab %>%
  mutate(Concentration = (Counts*Dilution)/Volume)

## Export as csv file
write.csv2(file = "PredictedCellsStrepsFab.csv", test_pred_streps_Fab)

### IS FM
## Select the fcs files based on which the model will be trained --> Streps from Ioanna's stock
fcs_names_streps_FM <- c("20201125_Fabian_OralOcular_ISFM_Experiment_OralIoanna_1000_SG_A7nezw.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralIoanna_1000_SG_B7.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralIoanna_1000_SG_C7.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralIoanna_1000_SG_D7.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralIoanna_1000_SG_E7.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralIoanna_1000_SG_F7.fcs")

Sample_Info_streps_FM <- Sample_Info %>% dplyr::filter(name %in% fcs_names_streps_FM)

Model_RF_streps_FM <- Phenoflow::RandomF_FCS(flowData_transformed_FM_gated[fcs_names_streps_FM], Sample_Info_streps_FM, target_label = "Strain", downsample = 24000, classification_type = "sample", param = paramRF , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for Streps from Fabian's stock
fcs_topre_streps_FM <- c("20201125_Fabian_OralOcular_ISFM_Experiment_OralFabian_1000_SG_A8.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralFabian_1000_SG_B8.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralFabian_1000_SG_C8.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralFabian_1000_SG_D8.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralFabian_1000_SG_E8.fcs",
                         "20201125_Fabian_OralOcular_ISFM_Experiment_OralFabian_1000_SG_F8.fcs")

flowData_topre_streps_FM <- flowData_transformed_FM_gated[fcs_topre_streps_FM]
test_pred_streps_FM <- RandomF_predict(x = Model_RF_streps_FM[[1]], new_data =  flowData_topre_streps_FM, cleanFCS = FALSE)

## Calculate cell concentrations of predictions
vol_FM2 <- data.frame(Sample_name = flowCore::sampleNames(flowData_transformed_FM),
                       Volume = vol_FM)

test_pred_streps_FM <- left_join(test_pred_streps_FM, vol_FM2, by = c("Sample" = "Sample_name"))
test_pred_streps_FM <- left_join(test_pred_streps_FM, Sample_Info_streps_FM, by = c("Sample" = "File_Name"))
test_pred_streps_FM <- test_pred_streps_FM %>%
  mutate(Concentration = (Counts*Dilution)/Volume)

## Export as csv file
write.csv2(file = "PredictedCellsStrepsFM.csv", test_pred_streps_FM)