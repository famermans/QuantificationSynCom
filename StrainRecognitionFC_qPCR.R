#### Clearing workspace, loading libraries, setting seed ----

## Clear environment and set working directory
rm(list = ls())
setwd("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM")

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

csv_mocks <- read.csv(file = "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/qPCR/data/qPCR_mocks.csv", header = T, sep = ";", stringsAsFactors = F)
csv_counts <- read.csv(file = "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/qPCR/data/Counts_mocks.csv", header = T, sep = ";", stringsAsFactors = F)


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

#### Visualization ----
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

plot3 <- ggplot(data = melted, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'dodge')+
  theme_bw()+
  facet_grid(~ Sample_Name)+
  labs(title = "Theoretical vs qPCR mocks",
       x = "Technique",
       y = "Cell concentration (cells/mL)",
       fill = "Species")+
#  scale_y_log10(breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000))
  scale_y_continuous(trans='log10', breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000))+
  coord_cartesian(ylim = c(10000, 10000000000))

print(plot3)

## Plot ony part of the mocks (mock 1, 2 and 8)
melted_sel <- filter(melted, Sample_Name == "M1"| Sample_Name == "M2"| Sample_Name == "M8")

plot4 <- ggplot(data = melted_sel, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'fill')+
  theme_bw()+
  facet_grid(~ Sample_Name)+
  labs(title = "Theoretical vs qPCR mocks",
       x = "Technique",
       y = "Relative abundance",
       fill = "Species")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 24, hjust = 0.5))

print(plot4)

plot5 <- ggplot(data = melted_sel, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'dodge')+
  theme_bw()+
  facet_grid(~ Sample_Name)+
  labs(title = "Theoretical vs qPCR mocks",
       x = "Technique",
       y = "Cell concentration (cells/mL)",
       fill = "Species")+
  scale_y_continuous(trans='log10', breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000))+
  coord_cartesian(ylim = c(10000, 10000000000))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 24, hjust = 0.5))

print(plot5)

## Other selection of samples (M1, M2, M3, M8)
melted_sel2 <- filter(melted, Sample_Name == "M1"| Sample_Name == "M2"| Sample_Name == "M8" | Sample_Name == "M3" & cat == "FCM")

plot_sel <- ggplot(data = melted_sel2, aes(x = Sample_Name, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'fill')+
  theme_bw()+
  labs(title = "Prediction mocks",
       x = "Mock",
       y = "Relative abundance",
       fill = "Species")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 24, hjust = 0.5))

print(plot_sel)

#### Visualization mocks with Illumina data ----
csv_mocks2 <- read.csv(file = "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/qPCR/data/qPCR_RF_mocks2.csv", header = T, sep = ";", stringsAsFactors = F)
csv_mocks2 <- csv_mocks2[, -1]

plot6 <- ggplot(data = csv_mocks2, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'dodge')+
  theme_bw()+
  facet_grid(~ Sample_Name)+
  labs(title = "Mock Communities",
       x = "Technique",
       y = "Cell concentration (cells/mL)",
       fill = "Species")+
  scale_y_continuous(trans='log10', breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000))+
  coord_cartesian(ylim = c(100000, 10000000000))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 24, hjust = 0.5))

print(plot6)

plot7 <- ggplot(data = csv_mocks2, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'fill')+
  theme_bw()+
  facet_grid(~ Sample_Name)+
  labs(title = "Mock Communities",
       x = "Technique",
       y = "Relative abundance",
       fill = "Species")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 90),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 24, hjust = 0.5))

print(plot7)


#### Visualization co-cultures ----

csv_cocultures <- read.csv(file = "/media/cmetfcm/Fabian/Oral_Microbiology/Strain_Recognition_FCM/qPCR/data/qPCR_cocultures2.csv", header = T, sep = ";", stringsAsFactors = F)

plot8 <- ggplot(data = csv_cocultures, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'dodge')+
  theme_bw()+
  facet_grid(~ Replicate)+
  labs(title = "Co-culture 2 (So, Fn, Pg)",
       x = "Technique",
       y = "Cell concentration (cells/mL)",
       fill = "Species")+
  scale_y_continuous(trans='log10', breaks = c(10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000, 10000000000))+
  coord_cartesian(ylim = c(100000, 10000000000))+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 24, hjust = 0.5))

print(plot8)

plot9 <- ggplot(data = csv_cocultures, aes(x = cat, y = concentration, fill = species))+
  geom_bar(stat = 'identity', position = 'fill')+
  theme_bw()+
  facet_grid(~ Replicate)+
  labs(title = "Co-culture 2 (So, Fn, Pg)",
       x = "Technique",
       y = "Relative abundance",
       fill = "Species")+
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        plot.title = element_text(size = 24, hjust = 0.5))

print(plot9)
