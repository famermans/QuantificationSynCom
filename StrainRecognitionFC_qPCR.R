#### Clearing workspace, loading libraries, setting seed ----

## Clear environment and set working directory
rm(list = ls())
setwd("/media/projects1/Fabian/Oral microbiome/StrainRecognitionFC/qPCR")

## Load libraries
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

set.seed(777)


#### Loading data ----

csv_mocks <- read.csv(file = "data/qPCR_mocks.csv", header = T, sep = ";", stringsAsFactors = F)
csv_counts <- read.csv(file = "data/Counts_mocks.csv", header = T, sep = ";", stringsAsFactors = F)

#### Calculations ----

## Gene copy number per species
CopyNumber_So <- 1
CopyNumber_Fn <- 5
CopyNumber_Pg <- 4
CopyNumber_Vp <- 4

## Parameters relevant for calculations
Volume <- 1               # mL
Dilution <- 10
VolumeDNA <- 25           # µL

### Calculating cell concentrations based on copy numbers from qPCR
mocks <- data.frame(Sample_Name = csv_mocks$Sample_Name,
                    Cells_So = csv_mocks$Quantity_Mean_So/CopyNumber_So,
                    CellsSD_So = csv_mocks$Quantity_SD_So/CopyNumber_So,
                    Cells_Fn = csv_mocks$Quantity_Mean_Fn/CopyNumber_Fn,
                    CellsSD_Fn = csv_mocks$Quantity_SD_Fn/CopyNumber_Fn,
                    Cells_Pg = csv_mocks$Quantity_Mean_Pg/CopyNumber_Pg,
                    CellsSD_Pg = csv_mocks$Quantity_SD_Pg/CopyNumber_Pg,
                    Cells_Vp = csv_mocks$Quantity_Mean_Vp/CopyNumber_Vp,
                    CellsSD_Vp = csv_mocks$Quantity_SD_Vp/CopyNumber_Vp)

mocks_concentration <- data.frame(Sample_Name = mocks$Sample_Name,
                                  Concentration_So = mocks$Cells_So*(Dilution*VolumeDNA/Volume),
                                  ConcentrationSD_So = mocks$CellsSD_So*(Dilution*VolumeDNA/Volume),
                                  Concentration_Fn = mocks$Cells_Fn*(Dilution*VolumeDNA/Volume),
                                  ConcentrationSD_Fn = mocks$CellsSD_Fn*(Dilution*VolumeDNA/Volume),
                                  Concentration_Pg = mocks$Cells_Pg*(Dilution*VolumeDNA/Volume),
                                  ConcentrationSD_Pg = mocks$CellsSD_Pg*(Dilution*VolumeDNA/Volume),
                                  Concentration_Vp = mocks$Cells_Vp*(Dilution*VolumeDNA/Volume),
                                  ConcentrationSD_Vp = mocks$CellsSD_Vp*(Dilution*VolumeDNA/Volume))

mocks_qPCR <- data.frame(Sample_Name = mocks_concentration$Sample_Name,
                         qPCR_So = mocks_concentration$Concentration_So,
                         qPCR_Fn = mocks_concentration$Concentration_Fn,
                         qPCR_Pg = mocks_concentration$Concentration_Pg,
                         qPCR_Vp = mocks_concentration$Concentration_Vp)

### Reshaping data for visualization purposes
## Mering data frames
data <- merge(mocks_qPCR, csv_counts, by = "Sample_Name")

## Melting data frame
melted <- melt(data, "Sample_Name")

## Create new column in melted data frame with category (= qPCR or flow cytometry)
melted$cat <- ''

melted[melted$variable == 'qPCR_So', ]$cat <- "qPCR"
melted[melted$variable == 'qPCR_Fn', ]$cat <- "qPCR"
melted[melted$variable == 'qPCR_Pg', ]$cat <- "qPCR"
melted[melted$variable == 'qPCR_Vp', ]$cat <- "qPCR"

melted[melted$variable == 'FCM_So', ]$cat <- "FCM"
melted[melted$variable == 'FCM_Fn', ]$cat <- "FCM"
melted[melted$variable == 'FCM_Pg', ]$cat <- "FCM"
melted[melted$variable == 'FCM_Vp', ]$cat <- "FCM"

## Create new column in melted dataframe with species name
melted$species <- ''

melted[melted$variable == 'qPCR_So', ]$species <- "So"
melted[melted$variable == 'FCM_So', ]$species <- "So"

melted[melted$variable == 'qPCR_Fn', ]$species <- "Fn"
melted[melted$variable == 'FCM_Fn', ]$species <- "Fn"

melted[melted$variable == 'qPCR_Pg', ]$species <- "Pg"
melted[melted$variable == 'FCM_Pg', ]$species <- "Pg"

melted[melted$variable == 'qPCR_Vp', ]$species <- "Vp"
melted[melted$variable == 'FCM_Vp', ]$species <- "Vp"

## Change name of 'value' column in data frame to 'concentration'
melted <- melted %>% 
  dplyr::rename(concentration = value)

### Visualization
plot1 <- ggplot(data = melted, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'stack')+
  theme_bw()+
  facet_grid(~ Sample_Name)+
  labs(title = "FCM vs qPCR mocks",
       x = "Technique",
       y = "Cell concentration (cells/mL)")

print(plot1)











