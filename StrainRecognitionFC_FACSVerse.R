#### Clearing workspace, loading libraries, setting seed ----

## Clear environment and set working directory
rm(list = ls())
setwd("/media/projects1/Fabian/Oral microbiome/StrainRecognitionFC/FACSVerse")

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

Datapath <- "/media/projects1/Fabian/Oral microbiome/StrainRecognitionFC/FACSVerse/fcsfiles"
fcsfiles <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = TRUE)
flowData <- flowCore::read.flowSet(files = fcsfiles, transformation = FALSE)
metadata <- read.table("/media/projects1/Fabian/Oral microbiome/StrainRecognitionFC/FACSVerse/201910_MetaData3.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)


#### Transformation of data ----

flowData_transformed <- transform(flowData,
                                  `PerCP-Cy5.5-A` = asinh(`PerCP-Cy5.5-A`), 
                                  `SSC-A` = asinh(`SSC-A`), 
                                  `FITC-A` = asinh(`FITC-A`), 
                                  `FSC-A` = asinh(`FSC-A`),
                                  `PerCP-Cy5.5-H` = asinh(`PerCP-Cy5.5-H`), 
                                  `SSC-H` = asinh(`SSC-H`), 
                                  `FITC-H` = asinh(`FITC-H`), 
                                  `FSC-H` = asinh(`FSC-H`))
param = c("PerCP-Cy5.5-A", "FITC-A","SSC-A","FSC-A", "PerCP-Cy5.5-H", "FITC-H","SSC-H","FSC-H")


#### Gating ----

sqrcut1 <- matrix(c(6.2,6.2,14,14,
                    4.5,7.6,15,5), ncol = 2, nrow = 4)
colnames(sqrcut1) <- c("FITC-A", "PerCP-Cy5.5-A")
polyGate1 <- polygonGate(.gate = sqrcut1, filterId = "Total Cells")

##  Gating quality check
xyplot(`PerCP-Cy5.5-A`~`FITC-A`, data = flowData_transformed[c(1, 12, 19, 22, 31, 42)], filter = polyGate1,
       scales = list(y = list(limits = c(0, 16)),
                   x = list(limits = c(4, 16))),
       axis = axis.default, nbin = 125, main = "QC Gating", xlab = "FITC-A", ylab = "PerCP-Cy5.5-A",
       par.strip.text = list(col = "white", font = 1, cex = 1), smooth = FALSE)

flowData_transformed_gated <- Subset(flowData_transformed, polyGate1)


#### Cell concentrations ----

### Cell counts
cells <- flowCore::filter(flowData_transformed, polyGate1)
TotalCount <- summary(cells); TotalCount <- toTable(TotalCount)

### Extracting volumes
# Volumes are in µL
vol <- as.numeric(flowCore::fsApply(flowData_transformed, FUN = function(x) x@description$VOL))/1000

### Concentrations
# Calculating concentrations --> Concentrations will be in cells/µL (dilution = 1000x)
cell_concentrations <- data.frame(Samples = flowCore::sampleNames(flowData_transformed),
                                  Concentration = TotalCount$true/vol)


#### Phenotypic Diversity Analysis ----

## Normalization of data
summary <- fsApply(x = flowData_transformed_gated, FUN = function(x) apply(x, 2, max), use.exprs = TRUE)
max = max(summary[, "FITC-A"])
max2 = max(summary[, "FITC-H"])
mytrans <- function(x) x/max
mytrans2 <- function(x) x/max2

flowData_transformed_2 <- transform(flowData_transformed_gated,
                                    `PerCP-Cy5.5-A` = mytrans(`PerCP-Cy5.5-A`), 
                                    `SSC-A` = mytrans(`SSC-A`), 
                                    `FITC-A` = mytrans(`FITC-A`), 
                                    `FSC-A` = mytrans(`FSC-A`),
                                    `PerCP-Cy5.5-H` = mytrans2(`PerCP-Cy5.5-H`), 
                                    `SSC-H` = mytrans2(`SSC-H`), 
                                    `FITC-H` = mytrans2(`FITC-H`), 
                                    `FSC-H` = mytrans2(`FSC-H`))

## Parameters taken into account for further calculations
param2 <- c("PerCP-Cy5.5-A", "FITC-A","SSC-A","FSC-A")

## Calculate fingerprint with bw = 0.01
fbasis <- flowBasis(flowData_transformed_2, param2, nbin = 128, 
                    bw = 0.01, normalize = function(x) x)

## Calculate ecological parameters from normalized fingerprint 
Diversity.fbasis <- Diversity(fbasis, d = 3, plot = FALSE, R = 999)
#Evenness.fbasis <- Evenness(fbasis, d = 3, plot = FALSE)
#Structural.organization.fbasis <- So(fbasis, d = 3, plot = FALSE)
#Coef.var.fbasis <- CV(fbasis, d = 3, plot = FALSE)

p_alphadiv <- ggplot(data = Diversity.fbasis, aes(x = as.numeric(as.character(metadata$Time)), y = D2, color = metadata$Strains))+
  geom_point(size = 5, alpha = 0.7)+
  geom_line()+
  facet_grid(metadata$State ~ ., scales = "free")+       # Creates different subplots
  theme_bw()+
  labs(color = "Strains", y = "Phenotypic diversity (D2)", x = "Hours", title = "Phenotypic diversity analysis")+
  geom_errorbar(aes(ymin=D2-sd.D2, ymax=D2+sd.D2), width=0.05, color="black")
print(p_alphadiv)

## Export ecological data to .csv file in the chosen directory
#write.csv2(file = "results.metrics.csv",
#           cbind(Diversity.fbasis, Evenness.fbasis,
#                 Structural.organization.fbasis,
#                 Coef.var.fbasis))

## Beta-diversity assessment of fingerprint
beta.div <- beta_div_fcm(fbasis, ord.type = "PCoA")

## Plot ordination
plot_beta_fcm(beta.div, color = metadata$State, shape = as.factor(metadata$Time), labels = list("State", "Timepoint")) + 
  theme_bw() +
  geom_point(size = 8, alpha = 0.5)


#### Training Random Forest classifier ----

name <- list.files(path = Datapath, recursive = TRUE, pattern = ".fcs", full.names = FALSE)
Sample_Info <- cbind(name, metadata)
vol <- data.frame(Sample_names = flowCore::sampleNames(flowData_transformed),
                  Volume = vol)
# Extract the total number of events per fcs file in order to calculate accuracy of model to predict strain
totalTable <- toTable(summary(flowCore::filter(flowData_transformed, polyGate1)))
totalTable <- left_join(totalTable, vol, by = c("sample" = "Sample_names"))
write.csv2(file="totalCells.csv", totalTable)

### Build random forest for So and Fn
fcs_names_SoFn <- c("So_10_3_20191015_160057.fcs",
                    "Fn_10_3_20191015_161647.fcs")
Sample_Info_SoFn <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFn)

Model_RF_SoFn <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFn], Sample_Info_SoFn, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures
fcs_topre_SoFn <- c("Fn_10_3_20191017_145607.fcs",
                    "So_10_3_20191017_143923.fcs",
                    "So_Fn_10_3_20191015_162124.fcs",
                    "So_Fn_Pi_10_3_20191015_162306.fcs",
                    "1_10_3_20191016_141301.fcs",
                    "1_10_3_20191017_143044.fcs")
flowData_topre_SoFn <- flowData_transformed_gated[fcs_topre_SoFn]

test_pred_SoFn <- RandomF_predict(x = Model_RF_SoFn[[1]], new_data = flowData_topre_SoFn, cleanFCS = FALSE)
test_pred_SoFn

## Export predictions
# Convert predictions to concentrations
test_pred_SoFn <- left_join(test_pred_SoFn, vol, by = c("Sample" = "Sample_names"))
test_pred_SoFn <- test_pred_SoFn %>% 
  mutate(Concentration = Counts/Volume)
# Produce csv files with concentrations
write.csv2(file="PredictedCells_SoFn.csv", test_pred_SoFn)


### Build random forest for So, Fn and Pi
fcs_names_SoFnPi <- c("So_10_3_20191015_160057.fcs",
                    "Fn_10_3_20191015_161647.fcs",
                    "Pi_10_3_20191015_161828.fcs")
Sample_Info_SoFnPi <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPi)

Model_RF_SoFnPi <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFnPi], Sample_Info_SoFnPi, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures
fcs_topre_SoFnPi <- c("Fn_10_3_20191017_145607.fcs",
                    "So_10_3_20191017_143923.fcs",
                    "Pi_10_3_20191017_145747.fcs",
                    "So_Fn_10_3_20191015_162124.fcs",
                    "So_Fn_Pi_10_3_20191015_162306.fcs",
                    "So_Fn_Pi_Vp_10_3_20191015_162448.fcs",
                    "1_10_3_20191016_141301.fcs",
                    "1_10_3_20191017_143044.fcs",
                    "2_10_3_20191016_141435.fcs",
                    "2_10_3_20191017_143159.fcs")
flowData_topre_SoFnPi <- flowData_transformed_gated[fcs_topre_SoFnPi]

test_pred_SoFnPi <- RandomF_predict(x = Model_RF_SoFnPi[[1]], new_data = flowData_topre_SoFnPi, cleanFCS = FALSE)
test_pred_SoFnPi

## Export predictions
# Convert predictions to concentrations
test_pred_SoFnPi <- left_join(test_pred_SoFnPi, vol, by = c("Sample" = "Sample_names"))
test_pred_SoFnPi <- test_pred_SoFnPi %>% 
  mutate(Concentration = Counts/Volume)
# Produce csv files with concentrations
write.csv2(file="PredictedCells_SoFnPi.csv", test_pred_SoFnPi)


### Build random forest for So, Fn, Pi and Vp
fcs_names_SoFnPiVp <- c("So_10_3_20191015_160057.fcs",
                      "Fn_10_3_20191015_161647.fcs",
                      "Pi_10_3_20191015_161828.fcs",
                      "Vp_10_3_20191015_161144.fcs")
Sample_Info_SoFnPiVp <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoFnPiVp)

Model_RF_SoFnPiVp <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoFnPiVp], Sample_Info_SoFnPiVp, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures
fcs_topre_SoFnPiVp <- c("Fn_10_3_20191017_145607.fcs",
                      "So_10_3_20191017_143923.fcs",
                      "Pi_10_3_20191017_145747.fcs",
                      "Vp_10_3_20191017_145142.fcs",
                      "So_Fn_10_3_20191015_162124.fcs",
                      "So_Fn_Pi_10_3_20191015_162306.fcs",
                      "So_Fn_Pi_Vp_10_3_20191015_162448.fcs",
                      "1_10_3_20191016_141301.fcs",
                      "1_10_3_20191017_143044.fcs",
                      "2_10_3_20191016_141435.fcs",
                      "2_10_3_20191017_143159.fcs",
                      "3_10_3_20191016_141617.fcs",
                      "3_10_3_20191017_143340.fcs")
flowData_topre_SoFnPiVp <- flowData_transformed_gated[fcs_topre_SoFnPiVp]

test_pred_SoFnPiVp <- RandomF_predict(x = Model_RF_SoFnPiVp[[1]], new_data = flowData_topre_SoFnPiVp, cleanFCS = FALSE)
test_pred_SoFnPiVp

## Export predictions
# Convert predictions to concentrations
test_pred_SoFnPiVp <- left_join(test_pred_SoFnPiVp, vol, by = c("Sample" = "Sample_names"))
test_pred_SoFnPiVp <- test_pred_SoFnPiVp %>% 
  mutate(Concentration = Counts/Volume)
# Produce csv files with concentrations
write.csv2(file="PredictedCells_SoFnPiVp.csv", test_pred_SoFnPiVp)


### Build random forest for So, Ssan, Ssal, Vp, Av, An
fcs_names_SoSsanSsalVpAvAn <- c("So_10_3_20191015_160057.fcs",
                                "Ssan_10_3_20191015_160415.fcs",
                                "Ssal_10_3_20191015_160238.fcs",
                                "Vp_10_3_20191015_161144.fcs",
                                "Av_10_3_20191015_161506.fcs",
                                "An_10_3_20191015_161325.fcs")
Sample_Info_SoSsanSsalVpAvAn <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoSsanSsalVpAvAn)

Model_RF_SoSsanSsalVpAvAn <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoSsanSsalVpAvAn], Sample_Info_SoSsanSsalVpAvAn, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures
fcs_topre_SoSsanSsalVpAvAn <- c("So_10_3_20191017_143923.fcs",
                                "Ssan_10_3_20191017_144300.fcs",
                                "Ssal_10_3_20191017_144126.fcs",
                                "Vp_10_3_20191017_145142.fcs",
                                "Av_10_3_20191017_145406.fcs",
                                "4_10_3_20191016_141808.fcs",
                                "4_10_3_20191017_143513.fcs")
flowData_topre_SoSsanSsalVpAvAn <- flowData_transformed_gated[fcs_topre_SoSsanSsalVpAvAn]

test_pred_SoSsanSsalVpAvAn <- RandomF_predict(x = Model_RF_SoSsanSsalVpAvAn[[1]], new_data = flowData_topre_SoSsanSsalVpAvAn, cleanFCS = FALSE)
test_pred_SoSsanSsalVpAvAn

## Export predictions
# Convert predictions to concentrations
test_pred_SoSsanSsalVpAvAn <- left_join(test_pred_SoSsanSsalVpAvAn, vol, by = c("Sample" = "Sample_names"))
test_pred_SoSsanSsalVpAvAn <- test_pred_SoSsanSsalVpAvAn %>% 
  mutate(Concentration = Counts/Volume)
# Produce csv files with concentrations
write.csv2(file="PredictedCells_SoSsanSsalVpAvAn.csv", test_pred_SoSsanSsalVpAvAn)


### Build random forest for So, Ssan, Ssal, Vp, Av, An, Smi, Sg
fcs_names_SoSsanSsalVpAvAnSmiSg <- c("So_10_3_20191015_160057.fcs",
                                "Ssan_10_3_20191015_160415.fcs",
                                "Ssal_10_3_20191015_160238.fcs",
                                "Vp_10_3_20191015_161144.fcs",
                                "Av_10_3_20191015_161506.fcs",
                                "An_10_3_20191015_161325.fcs",
                                "Smi_10_3_20191015_160639.fcs",
                                "Sg_10_3_20191015_160531.fcs")
Sample_Info_SoSsanSsalVpAvAnSmiSg <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SoSsanSsalVpAvAnSmiSg)

Model_RF_SoSsanSsalVpAvAnSmiSg <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SoSsanSsalVpAvAnSmiSg], Sample_Info_SoSsanSsalVpAvAnSmiSg, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures
fcs_topre_SoSsanSsalVpAvAnSmiSg <- c("So_10_3_20191017_143923.fcs",
                                     "Ssan_10_3_20191017_144300.fcs",
                                     "Ssal_10_3_20191017_144126.fcs",
                                     "Vp_10_3_20191017_145142.fcs",
                                     "Av_10_3_20191017_145406.fcs",
                                     "Smi_10_3_20191017_144619.fcs",
                                     "Sg_10_3_20191017_144441.fcs",
                                     "4_10_3_20191016_141808.fcs",
                                     "4_10_3_20191017_143513.fcs")
flowData_topre_SoSsanSsalVpAvAnSmiSg <- flowData_transformed_gated[fcs_topre_SoSsanSsalVpAvAnSmiSg]

test_pred_SoSsanSsalVpAvAnSmiSg <- RandomF_predict(x = Model_RF_SoSsanSsalVpAvAnSmiSg[[1]], new_data = flowData_topre_SoSsanSsalVpAvAnSmiSg, cleanFCS = FALSE)
test_pred_SoSsanSsalVpAvAnSmiSg

## Export predictions
# Convert predictions to concentrations
test_pred_SoSsanSsalVpAvAnSmiSg <- left_join(test_pred_SoSsanSsalVpAvAnSmiSg, vol, by = c("Sample" = "Sample_names"))
test_pred_SoSsanSsalVpAvAnSmiSg <- test_pred_SoSsanSsalVpAvAnSmiSg %>% 
  mutate(Concentration = Counts/Volume)
# Produce csv files with concentrations
write.csv2(file="PredictedCells_SoSsanSsalVpAvAn.csv", test_pred_SoSsanSsalVpAvAnSmiSg)


### Build random forest for Smu, Ssob, Fn, Pi, Aa
fcs_names_SmuSsobFnPiAa <- c("Smu_10_3_20191015_160820.fcs",
                             "Ssob_10_3_20191015_161002.fcs",
                             "Fn_10_3_20191015_161647.fcs",
                             "Pi_10_3_20191015_161828.fcs",
                             "Aa_10_3_20191015_161944.fcs")
Sample_Info_SmuSsobFnPiAa <- Sample_Info %>% dplyr::filter(name %in% fcs_names_SmuSsobFnPiAa)

Model_RF_SmuSsobFnPiAa <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_SmuSsobFnPiAa], Sample_Info_SmuSsobFnPiAa, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures
fcs_topre_SmuSsobFnPiAa <- c("Smu_10_3_20191017_144759.fcs",
                             "Ssob_10_3_20191017_144939.fcs",
                             "Fn_10_3_20191017_145607.fcs",
                             "Pi_10_3_20191017_145747.fcs",
                             "Aa_10_3_20191017_145853.fcs",
                             "5_10_3_20191016_141949.fcs",
                             "5_10_3_20191017_143645.fcs")
flowData_topre_SmuSsobFnPiAa <- flowData_transformed_gated[fcs_topre_SmuSsobFnPiAa]

test_pred_SmuSsobFnPiAa <- RandomF_predict(x = Model_RF_SmuSsobFnPiAa[[1]], new_data = flowData_topre_SmuSsobFnPiAa, cleanFCS = FALSE)
test_pred_SmuSsobFnPiAa

## Export predictions
# Convert predictions to concentrations
test_pred_SmuSsobFnPiAa <- left_join(test_pred_SmuSsobFnPiAa, vol, by = c("Sample" = "Sample_names"))
test_pred_SmuSsobFnPiAa <- test_pred_SmuSsobFnPiAa %>% 
  mutate(Concentration = Counts/Volume)
# Produce csv files with concentrations
write.csv2(file="PredictedCells_SmuSsobFnPiAa.csv", test_pred_SmuSsobFnPiAa)


### Build random forest for all strains
fcs_names_all <- c("Aa_10_3_20191015_161944.fcs",
                   "Av_10_3_20191015_161506.fcs",
                   "Fn_10_3_20191015_161647.fcs",
                   "Pi_10_3_20191015_161828.fcs",
                   "Sg_10_3_20191015_160531.fcs",
                   "Smi_10_3_20191015_160639.fcs",
                   "Smu_10_3_20191015_160820.fcs",
                   "So_10_3_20191015_160057.fcs",
                   "Ssal_10_3_20191015_160238.fcs",
                   "Ssan_10_3_20191015_160415.fcs",
                   "Ssob_10_3_20191015_161002.fcs",
                   "Vp_10_3_20191015_161144.fcs")
Sample_Info_all <- Sample_Info %>% dplyr::filter(name %in% fcs_names_all)

Model_RF_all <- Phenoflow::RandomF_FCS(flowData_transformed_gated[fcs_names_all], Sample_Info_all, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

## Make predictions for relevant mixtures and co-cultures
fcs_topre_all <- c("So_10_3_20191017_143923.fcs",
                   "Ssan_10_3_20191017_144300.fcs",
                   "Ssal_10_3_20191017_144126.fcs",
                   "Vp_10_3_20191017_145142.fcs",
                   "Av_10_3_20191017_145406.fcs",
                   "Smi_10_3_20191017_144619.fcs",
                   "Sg_10_3_20191017_144441.fcs",
                   "Smu_10_3_20191017_144759.fcs",
                   "Ssob_10_3_20191017_144939.fcs",
                   "Fn_10_3_20191017_145607.fcs",
                   "Pi_10_3_20191017_145747.fcs",
                   "Aa_10_3_20191017_145853.fcs",
                   "1_10_3_20191016_141301.fcs",
                   "1_10_3_20191017_143044.fcs",
                   "2_10_3_20191016_141435.fcs",
                   "2_10_3_20191017_143159.fcs",
                   "3_10_3_20191016_141617.fcs",
                   "3_10_3_20191017_143340.fcs",
                   "4_10_3_20191016_141808.fcs",
                   "4_10_3_20191017_143513.fcs",
                   "5_10_3_20191016_141949.fcs",
                   "5_10_3_20191017_143645.fcs",
                   "6_10_3_20191016_142112.fcs",
                   "6_10_3_20191017_143749.fcs")
flowData_topre_all <- flowData_transformed_gated[fcs_topre_all]

test_pred_all <- RandomF_predict(x = Model_RF_all[[1]], new_data = flowData_topre_all, cleanFCS = FALSE)
test_pred_all

## Export predictions
# Convert predictions to concentrations
test_pred_all <- left_join(test_pred_all, vol, by = c("Sample" = "Sample_names"))
test_pred_all <- test_pred_all %>% 
  mutate(Concentration = Counts/Volume)
# Produce csv files with concentrations
write.csv2(file="PredictedCells_all.csv", test_pred_all)




################################
# Some extra models #
################################

# Run random forest for all Streptococcus

# Select fcs files
fcs_names_streps <- c("Sg_10_3_20191015_160531.fcs",
                      "Smi_10_3_20191015_160639.fcs",
                      "Smu_10_3_20191015_160820.fcs",
                      "So_10_3_20191015_160057.fcs",
                      "Ssal_10_3_20191015_160238.fcs",
                      "Ssan_10_3_20191015_160415.fcs",
                      "Ssob_10_3_20191015_161002.fcs")
# Sample info has to contain a column called 'name' which matches the samplenames of the fcs files
Sample_Info_streps <- Sample_Info %>% dplyr::filter(name %in% fcs_names_streps)
# Random forest model
Model_RF_streps <- Phenoflow::RandomF_FCS(flowData_transformed[fcs_names_streps], Sample_Info_streps, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param, p_train = 0.75, seed = 777, cleanFCS = F, timesplit = 0.1, TimeChannel = "Time", plot_fig = T)

# Run random forest for all strains except Streptococcus

# Select fcs files
fcs_names_no_streps <- c("Aa_10_3_20191015_161944.fcs",
                         "Av_10_3_20191015_161506.fcs",
                         "Fn_10_3_20191015_161647.fcs",
                         "Pg_10_3_20191017_145248.fcs",
                         "Pi_10_3_20191015_161828.fcs",
                         "Vp_10_3_20191015_161144.fcs")
# Sample info has to contain a column called 'name' which matches the samplenames of the fcs files
Sample_Info_no_streps <- Sample_Info %>% dplyr::filter(name %in% fcs_names_no_streps)
# Random forest model
Model_RF_no_streps <- Phenoflow::RandomF_FCS(flowData_transformed[fcs_names_no_streps], Sample_Info_no_streps, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param, p_train = 0.75, seed = 777, cleanFCS = F, timesplit = 0.1, TimeChannel = "Time", plot_fig = T)

#### Can one get high performance when training based on genus level?

#Select the fcs files based on which the model will be trained- all strains expect An
fcs_names_genus <- c("Aa_10_3_20191015_161944.fcs",
                     "Av_10_3_20191015_161506.fcs",
                     "Fn_10_3_20191015_161647.fcs",
                     "Pg_10_3_20191017_145248.fcs",
                     "Pi_10_3_20191015_161828.fcs",
                     "Sg_10_3_20191015_160531.fcs",
                     "Smi_10_3_20191015_160639.fcs",
                     "Smu_10_3_20191015_160820.fcs",
                     "So_10_3_20191015_160057.fcs",
                     "Ssal_10_3_20191015_160238.fcs",
                     "Ssan_10_3_20191015_160415.fcs",
                     "Ssob_10_3_20191015_161002.fcs",
                     "Vp_10_3_20191015_161144.fcs")
# Sample info has to contain a column called 'name' which matches the samplenames of the fcs files
Sample_Info_genus <- Sample_Info %>% dplyr::filter(name %in% fcs_names_genus)
# Random forest model
Model_RF_13strains_genus <- Phenoflow::RandomF_FCS(flowData_transformed[fcs_names_genus], Sample_Info_genus, target_label = "Genus", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

### Compare previous model to model on species level, only including 2 streps

#Select the fcs files based on which the model will be trained- all strains expect An
fcs_names_8strains <- c("Aa_10_3_20191015_161944.fcs",
                        "Av_10_3_20191015_161506.fcs",
                        "Fn_10_3_20191015_161647.fcs",
                        "Pg_10_3_20191017_145248.fcs",
                        "Pi_10_3_20191015_161828.fcs",
                        "Smu_10_3_20191015_160820.fcs",
                        "So_10_3_20191015_160057.fcs",
                        "Vp_10_3_20191015_161144.fcs")
# Sample info has to contain a column called 'name' which matches the samplenames of the fcs files
Sample_Info_8strains <- Sample_Info %>% dplyr::filter(name %in% fcs_names_8strains)
# Random forest model
Model_RF_8strains <- Phenoflow::RandomF_FCS(flowData_transformed[fcs_names_8strains], Sample_Info_8strains, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)

### Compare previous model to model on species level, only including 1 strep

#Select the fcs files based on which the model will be trained- all strains expect An
fcs_names_7strains <- c("Aa_10_3_20191015_161944.fcs",
                        "Av_10_3_20191015_161506.fcs",
                        "Fn_10_3_20191015_161647.fcs",
                        "Pg_10_3_20191017_145248.fcs",
                        "Pi_10_3_20191015_161828.fcs",
                        "So_10_3_20191015_160057.fcs",
                        "Vp_10_3_20191015_161144.fcs")
# Sample info has to contain a column called 'name' which matches the samplenames of the fcs files
Sample_Info_7strains <- Sample_Info %>% dplyr::filter(name %in% fcs_names_7strains)
# Random forest model
Model_RF_7strains <- Phenoflow::RandomF_FCS(flowData_transformed[fcs_names_7strains], Sample_Info_7strains, target_label = "Strains", downsample = 10000, classification_type = "sample", param = param , p_train = 0.75, seed = 777, cleanFCS = FALSE, timesplit = 0.1, TimeChannel = "Time", plot_fig = TRUE)


