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
Volume_BHI2 <- 0.75       # mL
Volume_MM <- 1            # mL
Dilution <- 10
VolumeDNA <- 25           # µL
# I'm in doubt whether division by VolumeqPCR should be performed, as it is also done for the standard...
#VolumeqPCR <- 5           # µL

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

mocks_BHI2 <- mocks[1:5, ]
mocks_MM <- mocks[6:10, ]

mocks_concentration_BHI2 <- data.frame(Sample_Name = mocks_BHI2$Sample_Name,
                                   Concentration_So = mocks_BHI2$Cells_So*(Dilution*VolumeDNA/Volume_BHI2),
                                   ConcentrationSD_So = mocks_BHI2$CellsSD_So*(Dilution*VolumeDNA/Volume_BHI2),
                                   Concentration_Fn = mocks_BHI2$Cells_Fn*(Dilution*VolumeDNA/Volume_BHI2),
                                   ConcentrationSD_Fn = mocks_BHI2$CellsSD_Fn*(Dilution*VolumeDNA/Volume_BHI2),
                                   Concentration_Pg = mocks_BHI2$Cells_Pg*(Dilution*VolumeDNA/Volume_BHI2),
                                   ConcentrationSD_Pg = mocks_BHI2$CellsSD_Pg*(Dilution*VolumeDNA/Volume_BHI2),
                                   Concentration_Vp = mocks_BHI2$Cells_Vp*(Dilution*VolumeDNA/Volume_BHI2),
                                   ConcentrationSD_Vp = mocks_BHI2$CellsSD_Vp*(Dilution*VolumeDNA/Volume_BHI2))

mocks_concentration_MM <- data.frame(Sample_Name = mocks_MM$Sample_Name,
                                  Concentration_So = mocks_MM$Cells_So*(Dilution*VolumeDNA/Volume_MM),
                                  ConcentrationSD_So = mocks_MM$CellsSD_So*(Dilution*VolumeDNA/Volume_MM),
                                  Concentration_Fn = mocks_MM$Cells_Fn*(Dilution*VolumeDNA/Volume_MM),
                                  ConcentrationSD_Fn = mocks_MM$CellsSD_Fn*(Dilution*VolumeDNA/Volume_MM),
                                  Concentration_Pg = mocks_MM$Cells_Pg*(Dilution*VolumeDNA/Volume_MM),
                                  ConcentrationSD_Pg = mocks_MM$CellsSD_Pg*(Dilution*VolumeDNA/Volume_MM),
                                  Concentration_Vp = mocks_MM$Cells_Vp*(Dilution*VolumeDNA/Volume_MM),
                                  ConcentrationSD_Vp = mocks_MM$CellsSD_Vp*(Dilution*VolumeDNA/Volume_MM))

mocks_concentration <- rbind(mocks_concentration_BHI2, mocks_concentration_MM)

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

## Create new column in melted data frame with category (= qPCR or calculation based on FCM)
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

## Export to csv
#write.csv2(file = "qPCR_RF_mocks.csv", melted)

### Visualization
plot1 <- ggplot(data = melted, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'stack')+
  theme_bw()+
  facet_grid(~ Sample_Name)+
  labs(title = "Theoretical vs qPCR mocks",
       x = "Technique",
       y = "Cell concentration (cells/mL)",
       fill = "Species")

print(plot1)

plot2 <- ggplot(data = melted, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'fill')+
  theme_bw()+
  facet_grid(~ Sample_Name)+
  labs(title = "Theoretical vs qPCR mocks",
       x = "Technique",
       y = "Relative abundance",
       fill = "Species")

print(plot2)