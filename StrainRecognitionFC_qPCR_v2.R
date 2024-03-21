# 1. Clearing workspace, loading libraries, setting seed ----

# Clear environment and set working directory
rm(list = ls())
setwd("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM")

# Load libraries
#library("Phenoflow")
#library("flowViz")
#library("flowFDA")
#library("flowAI")
#library("vegan")
library("ggplot2")
library("RColorBrewer")
#library("ggrepel")
#library("ape")
#library("gridExtra")
#library("grid")
#library("scales")
#library("cowplot")
library("reshape2")
library("dplyr")
library("tidyverse")

seed <- 777
set.seed(seed)


# 2. Load data ----

qPCR_data <- read.csv(file = "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/qPCR/data/qPCR_DataForR.csv", header = T, sep = ";", stringsAsFactors = F)

csv_mocks <- read.csv(file = "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/qPCR/data/qPCR_mocks.csv", header = T, sep = ";", stringsAsFactors = F)
csv_counts <- read.csv(file = "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/qPCR/data/Counts_mocks.csv", header = T, sep = ";", stringsAsFactors = F)

qPCR_mocks <- qPCR_data[c(37:41), ]

qPCR_mocks <- qPCR_mocks[, c(1:13)]
rownames(qPCR_mocks) <- NULL

qPCR_cocult <- qPCR_data[c(1:36), ]

# Smu and Ssob were used in co-culture 4 instead of Smi and Ssan --> did not run qPCR for Smu and Ssob for these samples, so this co-culture is useless
# Pi was actually not Pi, so co-culture 5 and 6 are useless
qPCR_cocult <- qPCR_cocult[!(qPCR_cocult$Sample %in% c("Co4", "Co5", "Co6")), c(1:13)]
qPCR_cocult[is.na(qPCR_cocult)] <- 0
rownames(qPCR_cocult) <- NULL


# 3. Calculations qPCR ----

# Gene copy number per species
CopyNumber_So <- 1
CopyNumber_Fn <- 5
CopyNumber_Pg <- 4
CopyNumber_Vp <- 4

# Parameters relevant for calculations
Volume_Mocks_BHI <- 0.75       # mL
Volume_CoCultures_BHI <- 1     # mL
DilutionDNA <- 10              # DNA was diluted 10x before PCR reaction
VolumeDNA <- 25                # µL

# Calculating cell concentrations based on copy numbers from qPCR
mocks_qPCR <- data.frame(Sample_Name = qPCR_mocks$Sample,
                         Cells_So = qPCR_mocks$So/CopyNumber_So,
                         CellsSD_So = qPCR_mocks$SD_So/CopyNumber_So,
                         Cells_Fn = qPCR_mocks$Fn/CopyNumber_Fn,
                         CellsSD_Fn = qPCR_mocks$SD_Fn/CopyNumber_Fn,
                         Cells_Pg = qPCR_mocks$Pg/CopyNumber_Pg,
                         CellsSD_Pg = qPCR_mocks$SD_Pg/CopyNumber_Pg,
                         Cells_Vp = qPCR_mocks$Vp/CopyNumber_Vp,
                         CellsSD_Vp = qPCR_mocks$SD_Vp/CopyNumber_Vp)

cocult_qPCR <- data.frame(Sample_Name = qPCR_cocult$Sample,
                          Cells_So = qPCR_cocult$So/CopyNumber_So,
                          CellsSD_So = qPCR_cocult$SD_So/CopyNumber_So,
                          Cells_Fn = qPCR_cocult$Fn/CopyNumber_Fn,
                          CellsSD_Fn = qPCR_cocult$SD_Fn/CopyNumber_Fn,
                          Cells_Pg = qPCR_cocult$Pg/CopyNumber_Pg,
                          CellsSD_Pg = qPCR_cocult$SD_Pg/CopyNumber_Pg,
                          Cells_Vp = qPCR_cocult$Vp/CopyNumber_Vp,
                          CellsSD_Vp = qPCR_cocult$SD_Vp/CopyNumber_Vp)

mocks_concentration <- data.frame(Sample_Name = mocks_qPCR$Sample_Name,
                                  Concentration_So = (mocks_qPCR$Cells_So*DilutionDNA*VolumeDNA)/Volume_Mocks_BHI,
                                  ConcentrationSD_So = (mocks_qPCR$CellsSD_So*DilutionDNA*VolumeDNA)/Volume_Mocks_BHI,
                                  Concentration_Fn = (mocks_qPCR$Cells_Fn*DilutionDNA*VolumeDNA)/Volume_Mocks_BHI,
                                  ConcentrationSD_Fn = (mocks_qPCR$CellsSD_Fn*DilutionDNA*VolumeDNA)/Volume_Mocks_BHI,
                                  Concentration_Pg = (mocks_qPCR$Cells_Pg*DilutionDNA*VolumeDNA)/Volume_Mocks_BHI,
                                  ConcentrationSD_Pg = (mocks_qPCR$CellsSD_Pg*DilutionDNA*VolumeDNA)/Volume_Mocks_BHI,
                                  Concentration_Vp = (mocks_qPCR$Cells_Vp*DilutionDNA*VolumeDNA)/Volume_Mocks_BHI,
                                  ConcentrationSD_Vp = (mocks_qPCR$CellsSD_Vp*DilutionDNA*VolumeDNA)/Volume_Mocks_BHI)

cocult_concentration <- data.frame(Sample_Name = cocult_qPCR$Sample_Name,
                                   Concentration_So = (cocult_qPCR$Cells_So*DilutionDNA*VolumeDNA)/Volume_CoCultures_BHI,
                                   ConcentrationSD_So = (cocult_qPCR$CellsSD_So*DilutionDNA*VolumeDNA)/Volume_CoCultures_BHI,
                                   Concentration_Fn = (cocult_qPCR$Cells_Fn*DilutionDNA*VolumeDNA)/Volume_CoCultures_BHI,
                                   ConcentrationSD_Fn = (cocult_qPCR$CellsSD_Fn*DilutionDNA*VolumeDNA)/Volume_CoCultures_BHI,
                                   Concentration_Pg = (cocult_qPCR$Cells_Pg*DilutionDNA*VolumeDNA)/Volume_CoCultures_BHI,
                                   ConcentrationSD_Pg = (cocult_qPCR$CellsSD_Pg*DilutionDNA*VolumeDNA)/Volume_CoCultures_BHI,
                                   Concentration_Vp = (cocult_qPCR$Cells_Vp*DilutionDNA*VolumeDNA)/Volume_CoCultures_BHI,
                                   ConcentrationSD_Vp = (cocult_qPCR$CellsSD_Vp*DilutionDNA*VolumeDNA)/Volume_CoCultures_BHI)

conc_mocks_qPCR <- data.frame(Sample_Name = mocks_concentration$Sample_Name,
                              qPCR_So = mocks_concentration$Concentration_So,
                              qPCR_So_SD = mocks_concentration$ConcentrationSD_So,
                              qPCR_Fn = mocks_concentration$Concentration_Fn,
                              qPCR_Fn_SD = mocks_concentration$ConcentrationSD_Fn,
                              qPCR_Pg = mocks_concentration$Concentration_Pg,
                              qPCR_Pg_SD = mocks_concentration$ConcentrationSD_Pg,
                              qPCR_Vp = mocks_concentration$Concentration_Vp,
                              qPCR_Vp_SD = mocks_concentration$ConcentrationSD_Vp)

conc_cocult_qPCR <- data.frame(Sample_Name = cocult_concentration$Sample_Name,
                               qPCR_So = cocult_concentration$Concentration_So,
                               qPCR_So_SD = cocult_concentration$ConcentrationSD_So,
                               qPCR_Fn = cocult_concentration$Concentration_Fn,
                               qPCR_Fn_SD = cocult_concentration$ConcentrationSD_Fn,
                               qPCR_Pg = cocult_concentration$Concentration_Pg,
                               qPCR_Pg_SD = cocult_concentration$ConcentrationSD_Pg,
                               qPCR_Vp = cocult_concentration$Concentration_Vp,
                               qPCR_Vp_SD = cocult_concentration$ConcentrationSD_Vp)

saveRDS(object = conc_mocks_qPCR, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/qPCR_mocks.rds")
saveRDS(object = conc_cocult_qPCR, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/qPCR_cocult.rds")


