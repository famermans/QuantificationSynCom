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
library("cowplot")
library("reshape2")
library("dplyr")
library("tidyverse")
library("ggpubr")

seed <- 777
set.seed(seed)

source(file = "/Projects1/Fabian/paper_theme_fab.R")


# 2. Loading data ----

qPCR_mocks <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/qPCR_mocks.rds")
qPCR_cocult <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/qPCR_cocult.rds")

illumina <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/genustable_top13.rds")
metadata_illumina <- readxl::read_excel("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/DADA2/Metadata_LGC_NGS3818.xlsx", sheet = "Sheet1") %>% 
  as.data.frame()

FCM_SoFn <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFn_pooled_50kevents.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_SoFnPg <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFnPg_pooled_50kevents.csv", header = T, sep = ";", stringsAsFactors = F)
FCM_SoFnPgVp <- read.csv(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/PredictionsRF/PredictedCellsSoFnPgVp_pooled_50kevents.csv", header = T, sep = ";", stringsAsFactors = F)

theoretical_mocks <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/theoretical_mocks.rds")

# Change , to . in FCM dataframes
FCM_SoFn$Concentration <- gsub(",", ".", FCM_SoFn$Concentration)
FCM_SoFnPg$Concentration <- gsub(",", ".", FCM_SoFnPg$Concentration)
FCM_SoFnPgVp$Concentration <- gsub(",", ".", FCM_SoFnPgVp$Concentration)

FCM_SoFn$Concentration <- as.numeric(FCM_SoFn$Concentration)
FCM_SoFnPg$Concentration <- as.numeric(FCM_SoFnPg$Concentration)
FCM_SoFnPgVp$Concentration <- as.numeric(FCM_SoFnPgVp$Concentration)

# Concentrations are in cells/µL -> change to cells/mL
FCM_SoFn$Concentration <- FCM_SoFn$Concentration*1000
FCM_SoFnPg$Concentration <- FCM_SoFnPg$Concentration*1000
FCM_SoFnPgVp$Concentration <- FCM_SoFnPgVp$Concentration*1000


# 3. Calculate means, relative abundances and format data into one data frame ----

## 3.1. 50k cells ----
theoretical_mocks_rect <- reshape2::dcast(theoretical_mocks, Mix_ID ~ Strain)
theoretical_mocks_rect$Technique <- "Theoretical"
theoretical_mocks_rect[is.na(theoretical_mocks_rect)] <- 0

theoretical_mocks_prop <- sweep(as.matrix(theoretical_mocks_rect[, -c(1, 16)]), 1, rowSums(theoretical_mocks_rect[, -c(1, 16)]), FUN = "/") %>% 
  as.data.frame()

theoretical_mocks_prop$Type <- "Mock"
theoretical_mocks_prop$ID <- theoretical_mocks_rect$Mix_ID
theoretical_mocks_prop$Technique <- "Theoretical"


qPCR_mocks_means <- qPCR_mocks[, c(1, 2, 4, 6, 8)]
colnames(qPCR_mocks_means) <- c("ID", "So", "Fn", "Pg", "Vp")

qPCR_mocks_means_prop <- sweep(as.matrix(qPCR_mocks_means[, -c(1)]), 1, rowSums(qPCR_mocks_means[, -c(1)]), FUN = "/") %>% 
  as.data.frame()

qPCR_mocks_means_prop$ID <- qPCR_mocks_means$ID
qPCR_mocks_means_prop$Type <- "Mock"
qPCR_mocks_means_prop$Technique <- "qPCR"

qPCR_cocult_means <- qPCR_cocult[, c(1, 2, 4, 6, 8)]
colnames(qPCR_cocult_means) <- c("ID", "So", "Fn", "Pg", "Vp")

qPCR_cocult_means_prop <- qPCR_cocult_means
qPCR_cocult_means_prop[, -c(1)] <- sweep(as.matrix(qPCR_cocult_means_prop[, -c(1)]), 1, rowSums(qPCR_cocult_means_prop[, -c(1)]), FUN = "/") %>% 
  as.data.frame()
qPCR_cocult_means_prop$Type <- "Co-culture"
qPCR_cocult_means_prop$Technique <- "qPCR"

qPCR_means_prop <- rbind(qPCR_mocks_means_prop, qPCR_cocult_means_prop)


FCM_SoFn_mean <- FCM_SoFn
FCM_SoFn_mean$Replicate[is.na(FCM_SoFn_mean$Replicate)] <- "Z"
FCM_SoFn_mean$Timepoint[is.na(FCM_SoFn_mean$Timepoint)] <- "0h"
FCM_SoFn_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFn_mean, FUN = mean)
FCM_SoFn_mean$ID <- paste(FCM_SoFn_mean$Strain, FCM_SoFn_mean$Replicate, FCM_SoFn_mean$Timepoint, sep = "_")
FCM_SoFn_mean$ID <- gsub("_Z_0h", "", FCM_SoFn_mean$ID)
#FCM_SoFn_mean$Predicted_label <- gsub("\\*", "", FCM_SoFn_mean$Predicted_label)
FCM_SoFn_mean <- FCM_SoFn_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFn_mean_rect <- reshape2::dcast(FCM_SoFn_mean, ID ~ Predicted_label)
FCM_SoFn_mean_rect$Type <- c(rep("Co-culture", 18), rep("Mock", 5))
FCM_SoFn_mean_prop <- FCM_SoFn_mean_rect
FCM_SoFn_mean_prop[, -c(1, 4)] <- sweep(as.matrix(FCM_SoFn_mean_prop[, -c(1, 4)]), 1, rowSums(FCM_SoFn_mean_prop[, -c(1, 4)]), FUN = "/") %>% 
  as.data.frame()
FCM_SoFn_mean_prop$Technique <- "FCM_Model_SoFn"

FCM_SoFnPg_mean <- FCM_SoFnPg
FCM_SoFnPg_mean$Replicate[is.na(FCM_SoFnPg_mean$Replicate)] <- "Z"
FCM_SoFnPg_mean$Timepoint[is.na(FCM_SoFnPg_mean$Timepoint)] <- "0h"
FCM_SoFnPg_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFnPg_mean, FUN = mean)
FCM_SoFnPg_mean$ID <- paste(FCM_SoFnPg_mean$Strain, FCM_SoFnPg_mean$Replicate, FCM_SoFnPg_mean$Timepoint, sep = "_")
FCM_SoFnPg_mean$ID <- gsub("_Z_0h", "", FCM_SoFnPg_mean$ID)
#FCM_SoFnPg_mean$Predicted_label <- gsub("\\*", "", FCM_SoFnPg_mean$Predicted_label)
FCM_SoFnPg_mean <- FCM_SoFnPg_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFnPg_mean_rect <- reshape2::dcast(FCM_SoFnPg_mean, ID ~ Predicted_label)
FCM_SoFnPg_mean_rect$Type <- c(rep("Co-culture", 18), rep("Mock", 5))
FCM_SoFnPg_mean_prop <- FCM_SoFnPg_mean_rect
FCM_SoFnPg_mean_prop[, -c(1, 5)] <- sweep(as.matrix(FCM_SoFnPg_mean_prop[, -c(1, 5)]), 1, rowSums(FCM_SoFnPg_mean_prop[, -c(1, 5)]), FUN = "/") %>% 
  as.data.frame()
FCM_SoFnPg_mean_prop$Technique <- "FCM_Model_SoFnPg"

FCM_SoFnPgVp_mean <- FCM_SoFnPgVp
FCM_SoFnPgVp_mean$Replicate[is.na(FCM_SoFnPgVp_mean$Replicate)] <- "Z"
FCM_SoFnPgVp_mean$Timepoint[is.na(FCM_SoFnPgVp_mean$Timepoint)] <- "0h"
FCM_SoFnPgVp_mean <- aggregate(Concentration ~ Predicted_label + Strain + Replicate + Timepoint, data = FCM_SoFnPgVp_mean, FUN = mean)
FCM_SoFnPgVp_mean$ID <- paste(FCM_SoFnPgVp_mean$Strain, FCM_SoFnPgVp_mean$Replicate, FCM_SoFnPgVp_mean$Timepoint, sep = "_")
FCM_SoFnPgVp_mean$ID <- gsub("_Z_0h", "", FCM_SoFnPgVp_mean$ID)
#FCM_SoFnPgVp_mean$Predicted_label <- gsub("\\*", "", FCM_SoFnPgVp_mean$Predicted_label)
FCM_SoFnPgVp_mean <- FCM_SoFnPgVp_mean[, c("ID", "Predicted_label", "Concentration")]
FCM_SoFnPgVp_mean_rect <- reshape2::dcast(FCM_SoFnPgVp_mean, ID ~ Predicted_label)
FCM_SoFnPgVp_mean_rect$Type <- c(rep("Co-culture", 18), rep("Mock", 5))
FCM_SoFnPgVp_mean_prop <- FCM_SoFnPgVp_mean_rect
FCM_SoFnPgVp_mean_prop[, -c(1, 6)] <- sweep(as.matrix(FCM_SoFnPgVp_mean_prop[, -c(1, 6)]), 1, rowSums(FCM_SoFnPgVp_mean_prop[, -c(1, 6)]), FUN = "/") %>% 
  as.data.frame()
FCM_SoFnPgVp_mean_prop$Technique <- "FCM_Model_SoFnPgVp"

FCM_mean_prop <- dplyr::bind_rows(FCM_SoFn_mean_prop, FCM_SoFnPg_mean_prop, FCM_SoFnPgVp_mean_prop)
FCM_mean_prop[is.na(FCM_mean_prop)] <- 0

FCM_SoFn_mean_rect$Technique <- "FCM_Model_SoFn"
FCM_SoFnPg_mean_rect$Technique <- "FCM_Model_SoFnPg"
FCM_SoFnPgVp_mean_rect$Technique <- "FCM_Model_SoFnPgVp"

FCM_mean <- dplyr::bind_rows(FCM_SoFn_mean_rect, FCM_SoFnPg_mean_rect, FCM_SoFnPgVp_mean_rect)
FCM_mean[is.na(FCM_mean)] <- 0


illumina_FM <- illumina[, !names(illumina) %in% c("MOCK01", "MOCK02", "FM388_empty")] # Remove unwanted samples
illumina_FM <- illumina_FM[rowSums(illumina_FM != 0) > 0, ] # Remove genera that do not contain any observations
illumina_FM <- as.data.frame(t(illumina_FM))
illumina_FM[illumina_FM < 0] <- 0 # Replace negative values with zero
illumina_FM$Other <- illumina_FM$Staphylococcus + illumina_FM$Other # Add Staphylococcus to Other column as this is not a strain of interest
illumina_FM <- illumina_FM[, !names(illumina_FM) %in% c("Staphylococcus")]
colnames(illumina_FM) <- c("Fn", "Pg", "So", "Vp", "Other")
illumina_FM <- illumina_FM/100
illumina_prop <- illumina_FM
illumina_prop$Technique <- "Illumina"
illumina_prop$Name_LGC <- row.names(illumina_prop)
row.names(illumina_prop) <- NULL
metadata_illumina_merge <- metadata_illumina[-c(1, 2, 20), c("Name_LGC", "Sample", "Replicate", "Timepoint", "Type")]
row.names(metadata_illumina_merge) <- NULL
metadata_illumina_merge$ID <- paste(metadata_illumina_merge$Sample, metadata_illumina_merge$Replicate, metadata_illumina_merge$Timepoint, sep = "_")
metadata_illumina_merge <- metadata_illumina_merge[, c("Name_LGC", "ID", "Type")]
metadata_illumina_merge[c(1:6, 17), "ID"] <- c("Mix1_A1", "Mix2", "Mix3", "Mix8", "Mix9", "Mix1_MM", "Mix1_A2")
illumina_prop <- merge(illumina_prop, metadata_illumina_merge, by = "Name_LGC")
illumina_prop <- illumina_prop[, -c(1)]


all_data <- dplyr::bind_rows(theoretical_mocks_prop, qPCR_means_prop, FCM_mean_prop, illumina_prop)
all_data[is.na(all_data)] <- 0

all_data2 <- subset(all_data, ID != "Mix1_MM" & ID != "Mix1_A1" & ID != "Mix1_A2" & ID != "Mix4" & ID != "Mix5" & ID != "Mix6" & ID != "Mix7")
all_data2_illumina_mix1 <- subset(all_data, ID == "Mix1_A1" | ID == "Mix1_A2")
all_data2_illumina_mix1 <- colMeans(all_data2_illumina_mix1[, c(1:14, 18)])
all_data2_illumina_mix1 <- data.frame(t(all_data2_illumina_mix1))
all_data2_illumina_mix1$Type <- "Mock"
all_data2_illumina_mix1$Technique <- "Illumina"
all_data2_illumina_mix1$ID <- "Mix1"
all_data2 <- dplyr::bind_rows(all_data2, all_data2_illumina_mix1)
row.names(all_data2) <- NULL
all_data2 <- all_data2[, colSums(all_data2 != 0) > 0]

# Perform copy number correction for So, Fn, Pg and Vp
copyNumber_So <- 4
copyNumber_Fn <- 5
copyNumber_Pg <- 4
copyNumber_Vp <- 4

illumina_prop_copycorrected <- subset(all_data2, Technique == "Illumina")
illumina_prop_copycorrected$So <- illumina_prop_copycorrected$So/copyNumber_So
illumina_prop_copycorrected$Fn <- illumina_prop_copycorrected$Fn/copyNumber_Fn
illumina_prop_copycorrected$Pg <- illumina_prop_copycorrected$Pg/copyNumber_Pg
illumina_prop_copycorrected$Vp <- illumina_prop_copycorrected$Vp/copyNumber_Vp

illumina_prop_copycorrected[, -c(5:7)] <- sweep(as.matrix(illumina_prop_copycorrected[, -c(5:7)]), 1, rowSums(illumina_prop_copycorrected[, -c(5:7)]), FUN = "/") %>% 
  as.data.frame()

all_data2_copycorrected <- subset(all_data2, Technique != "Illumina")
all_data2_copycorrected <- dplyr::bind_rows(all_data2_copycorrected, illumina_prop_copycorrected)

## 3.2. Pooled models FCM ----
FCM_mean_prop_12k <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/FCM_mean_prop_12328cells.rds")

FCM_mean_prop_12k_2 <- FCM_mean_prop_12k[(FCM_mean_prop_12k$ID == "Mix1" & FCM_mean_prop_12k$Model == "SoFn") |
                                           (FCM_mean_prop_12k$ID == "Mix2" & FCM_mean_prop_12k$Model == "SoFnPg") |
                                           (FCM_mean_prop_12k$ID == "Mix3" & FCM_mean_prop_12k$Model == "SoFnPgVp") |
                                           (FCM_mean_prop_12k$ID == "Mix8" & FCM_mean_prop_12k$Model == "SoFn") |
                                           (FCM_mean_prop_12k$ID == "Mix9" & FCM_mean_prop_12k$Model == "SoFn") |
                                           FCM_mean_prop_12k$ID == "Mix4" | FCM_mean_prop_12k$ID == "Mix5" | FCM_mean_prop_12k$ID == "Mix6" | FCM_mean_prop_12k$ID == "Mix7", ]
FCM_mean_prop_12k_2 <- FCM_mean_prop_12k_2[, c(1:3, 5:16)]
FCM_mean_prop_12k_2$Technique <- "FCM"

theoretical_mocks_prop2 <- theoretical_mocks_prop[, c(1:14, 16:17)]

mocks_12k <- rbind(FCM_mean_prop_12k_2, theoretical_mocks_prop2)
row.names(mocks_12k) <- NULL

# Predictions in silico mocks
FCM_silico <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/predictions_mocks_silico.rds")
FCM_silico_mean <- FCM_silico[, c("ID", "Aa", "An", "Av", "Fn", "Pg", "Pi", "Sg", "Smi", "Smu", "So", "Ssal", "Ssan", "Ssob", "Vp")]
FCM_silico_mean_prop <- FCM_silico_mean
FCM_silico_mean_prop[, -c(1)] <- sweep(as.matrix(FCM_silico_mean_prop[, -c(1)]), 1, rowSums(FCM_silico_mean_prop[, -c(1)]), FUN = "/")

FCM_mean_prop_12k_2$Technique <- "FCM - In vitro"
FCM_silico_mean_prop$Technique <- "FCM - In silico"
mocks_12k_2 <- rbind(FCM_mean_prop_12k_2, FCM_silico_mean_prop, theoretical_mocks_prop2)
row.names(mocks_12k_2) <- NULL
mocks_12k_2$ID <- gsub("Mix", "Mock ", mocks_12k_2$ID)


# 4. Calculate performance mocks compared to theoretical ----

# Calculate RMSE
# Format data theoretical
theoretical_longer <- tidyr::pivot_longer(theoretical_mocks_prop[, c(1:14, 16)], cols = -ID, names_to = "Species", values_to = "Actual")

# FCM
RMSE_FCM <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/RMSE_FCM_50k.rds")

# qPCR
qPCR_longer <- tidyr::pivot_longer(qPCR_mocks_means_prop[, c(1:5)], cols = -ID, names_to = "Species", values_to = "qPCR")

merged_qPCR <- merge(qPCR_longer, theoretical_longer, by = c("ID", "Species"))

IDs_qPCR <- unique(merged_qPCR$ID)
RMSE_qPCR <- data.frame(ID = IDs_qPCR,
                        RMSE_qPCR = numeric(length(IDs_qPCR)))

for (i in seq_along(IDs_qPCR)) {
  ID <- IDs_qPCR[i]
  subset_df <- merged_qPCR[merged_qPCR$ID == ID, ]
  RMSE_qPCR$RMSE_qPCR[i] <- sqrt(mean((subset_df$qPCR - subset_df$Actual)^2))
}


# Illumina
illumina_mock_prop <- subset(illumina_prop_copycorrected, Type == "Mock")
illumina_longer <- tidyr::pivot_longer(illumina_mock_prop[, c(1:4, 6)], cols = -ID, names_to = "Species", values_to = "illumina")

merged_illumina <- merge(illumina_longer, theoretical_longer, by = c("ID", "Species"))

IDs_illumina <- unique(merged_illumina$ID)
RMSE_illumina <- data.frame(ID = IDs_illumina,
                            RMSE_illumina = numeric(length(IDs_illumina)))

for (i in seq_along(IDs_illumina)) {
  ID <- IDs_illumina[i]
  subset_df <- merged_illumina[merged_illumina$ID == ID, ]
  RMSE_illumina$RMSE_illumina[i] <- sqrt(mean((subset_df$illumina - subset_df$Actual)^2))
}


RMSE <- merge(RMSE_FCM, RMSE_qPCR, by = "ID", all = TRUE)
RMSE <- merge(RMSE, RMSE_illumina, by = "ID", all = TRUE)

# In silico mocks
RMSE_FCM_invitro <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/RMSE_FCM.rds")
RMSE_FCM_invitro <- data.frame(ID = RMSE_FCM_invitro$ID,
                               RMSE_invitro = c(RMSE_FCM_invitro[1, 2],
                                                RMSE_FCM_invitro[2, 3],
                                                RMSE_FCM_invitro[3, 4],
                                                RMSE_FCM_invitro[4, 5],
                                                RMSE_FCM_invitro[5, 6],
                                                RMSE_FCM_invitro[6, 7],
                                                RMSE_FCM_invitro[7, 8],
                                                RMSE_FCM_invitro[8, 2],
                                                RMSE_FCM_invitro[9, 2]))

FCM_insilico <- mocks_12k_2[mocks_12k_2$Technique == "FCM - In silico", ]
row.names(FCM_insilico) <- NULL
FCM_insilico_longer <- tidyr::pivot_longer(FCM_insilico[, c(1:15)], cols = -ID, names_to = "Species", values_to = "FCM_Insilico")
FCM_insilico_longer$ID <- gsub("Mock ", "Mix", FCM_insilico_longer$ID)

merged_FCM_insilico <- merge(FCM_insilico_longer, theoretical_longer, by = c("ID", "Species"))

IDs_FCM_insilico <- unique(merged_FCM_insilico$ID)
RMSE_FCM_insilico <- data.frame(ID = IDs_FCM_insilico,
                                RMSE_insilico = numeric(length(IDs_FCM_insilico)))

for (i in seq_along(IDs_FCM_insilico)) {
  ID <- IDs_FCM_insilico[i]
  subset_df <- merged_FCM_insilico[merged_FCM_insilico$ID == ID, ]
  RMSE_FCM_insilico$RMSE_insilico[i] <- sqrt(mean((subset_df$FCM_Insilico - subset_df$Actual)^2))
}

# RMSE 12k cells (in vitro vs in silico)
RMSE_FCM_12k <- merge(RMSE_FCM_invitro, RMSE_FCM_insilico, by = "ID", all = TRUE)
RMSE_FCM_12k$ID <- gsub("Mix", "Mock ", RMSE_FCM_12k$ID)


# 5. Visualization ----

## 5.1. Relative abundance ----
all_data_melted <- reshape2::melt(all_data, id.vars = c("ID", "Type", "Technique"), variable.name = c("Strain"), value.name = c("Relative_abundance"))
all_data2_melted <- reshape2::melt(all_data2, id.vars = c("ID", "Type", "Technique"), variable.name = c("Strain"), value.name = c("Relative_abundance"))
all_data2_copycorrected_melted <- reshape2::melt(all_data2_copycorrected, id.vars = c("ID", "Type", "Technique"), variable.name = c("Strain"), value.name = c("Relative_abundance"))

all_data_melted$Relative_abundance <- all_data_melted$Relative_abundance*100
all_data2_melted$Relative_abundance <- all_data2_melted$Relative_abundance*100
all_data2_copycorrected_melted$Relative_abundance <- all_data2_copycorrected_melted$Relative_abundance*100

plot_relative_all <- ggplot(data = all_data_melted, aes(x = ID, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_grid(Type ~ Technique) +
  theme_bw() +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(plot_relative_all)

plot_relative_all_sample <- ggplot(data = all_data_melted, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap( ~ ID) +
  theme_bw() +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(plot_relative_all_sample)

plot_relative_all2 <- ggplot(data = all_data2_melted, aes(x = ID, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_grid(Type ~ Technique) +
  theme_bw() +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(plot_relative_all2)

plot_relative_all2_sample <- ggplot(data = all_data2_melted, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap( ~ ID) +
  theme_bw() +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(plot_relative_all2_sample)

plot_relative_all2_copycorrected <- ggplot(data = all_data2_copycorrected_melted, aes(x = ID, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_grid(Type ~ Technique) +
  theme_bw() +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(plot_relative_all2_copycorrected)

plot_relative_all2_copycorrected_sample <- ggplot(data = all_data2_copycorrected_melted, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap( ~ ID) +
  theme_bw() +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
print(plot_relative_all2_copycorrected_sample)


### 5.1.1. Mocks ----
all_data2_mocks <- subset(all_data2_copycorrected_melted, Type == "Mock")
all_data2_mocks$ID <- gsub("Mix1", "Mock 1", all_data2_mocks$ID)
all_data2_mocks$ID <- gsub("Mix2", "Mock 2", all_data2_mocks$ID)
all_data2_mocks$ID <- gsub("Mix3", "Mock 3", all_data2_mocks$ID)
all_data2_mocks$ID <- gsub("Mix8", "Mock 8", all_data2_mocks$ID)
all_data2_mocks$ID <- gsub("Mix9", "Mock 9", all_data2_mocks$ID)
all_data2_mocks$Technique <- gsub("Theoretical", "Actual", all_data2_mocks$Technique)

plot_mocks <- ggplot(data = all_data2_mocks, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, nrow = length(unique(ID))) +
  theme_bw() +
  labs(x = "Technique", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#A3A500", "So" = "#00BF7D", "Vp" = "#00B0F6", "Other" = "#E76BF3"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")), "Other" = "Other")) +
  scale_x_discrete(labels = c("Actual" = "Theoretical", "FCM_Model_SoFn" = "FCM SoFn", "FCM_Model_SoFnPg" = "FCM SoFnPg", "FCM_Model_SoFnPgVp" = "FCM SoFnPgVp", "Illumina" = "Illumina", "qPCR" = "qPCR"))
print(plot_mocks)

plot_mocks <- ggplot(data = all_data2_mocks, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, nrow = length(unique(ID))) +
  theme_bw() +
  labs(x = "Technique", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800", "Other" = "#E76BF3"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")), "Other" = "Other")) +
  scale_x_discrete(labels = c("Actual" = "Theoretical", "FCM_Model_SoFn" = "FCM SoFn", "FCM_Model_SoFnPg" = "FCM SoFnPg", "FCM_Model_SoFnPgVp" = "FCM SoFnPgVp", "Illumina" = "Illumina", "qPCR" = "qPCR"))
print(plot_mocks)

mocks_12k_melted <- reshape2::melt(mocks_12k, id.vars = c("ID", "Technique"), variable.name = c("Strain"), value.name = c("Relative_abundance"))
mocks_12k_melted$ID <- gsub("Mix", "Mock ", mocks_12k_melted$ID)
mocks_12k_melted$Technique <- gsub("Theoretical", "Actual", mocks_12k_melted$Technique)
mocks_12k_melted$Relative_abundance <- mocks_12k_melted$Relative_abundance*100

plot_mocks_12k <- ggplot(data = mocks_12k_melted, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, nrow = length(unique(ID))) +
  theme_bw() +
  labs(x = "Technique", y = "Relative abundance (%)") +
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
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800", "An" = "#53B400", "Av" = "#00BC56", "Aa" = "#00C094", "Sg" = "#00BFC4", "Smi" = "#00B6EB", "Smu" = "#06A4FF", "Pi" = "#A58AFF", "Ssal" = "#DF70F8", "Ssan" = "#FB61D7", "Ssob" = "#FF66A8"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")), "An" = expression(italic("A. naeslundii")), "Av" = expression(italic("A. viscosus")), "Aa" = expression(italic("A. actinomycetemcomitans")), "Sg" = expression(italic("S. gordonii")), "Smi" = expression(italic("S. mitis")), "Smu" = expression(italic("S. mutans")), "Pi" = expression(italic("P. intermedia")), "Ssal" = expression(italic("S. salivarius")), "Ssan" = expression(italic("S. sanguinis")), "Ssob" = expression(italic("S. sobrinus")))) +
  scale_x_discrete(labels = c("Actual" = "Actual", "FCM" = "FCM"))
print(plot_mocks_12k)

mocks_12k_2_melted <- reshape2::melt(mocks_12k_2, id.vars = c("ID", "Technique"), variable.name = c("Strain"), value.name = c("Relative_abundance"))
mocks_12k_2_melted$Technique <- gsub("Theoretical", "Actual", mocks_12k_2_melted$Technique)
mocks_12k_2_melted$Relative_abundance <- mocks_12k_2_melted$Relative_abundance*100

plot_mocks_12k_2 <- ggplot(data = mocks_12k_2_melted, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, nrow = length(unique(ID))) +
  theme_bw() +
  labs(x = "Technique", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800", "An" = "#53B400", "Av" = "#00BC56", "Aa" = "#00C094", "Sg" = "#00BFC4", "Smi" = "#00B6EB", "Smu" = "#06A4FF", "Pi" = "#A58AFF", "Ssal" = "#DF70F8", "Ssan" = "#FB61D7", "Ssob" = "#FF66A8"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")), "An" = expression(italic("A. naeslundii")), "Av" = expression(italic("A. viscosus")), "Aa" = expression(italic("A. actinomycetemcomitans")), "Sg" = expression(italic("S. gordonii")), "Smi" = expression(italic("S. mitis")), "Smu" = expression(italic("S. mutans")), "Pi" = expression(italic("P. intermedia")), "Ssal" = expression(italic("S. salivarius")), "Ssan" = expression(italic("S. sanguinis")), "Ssob" = expression(italic("S. sobrinus")))) +
  scale_x_discrete(labels = c("Actual" = "Theoretical", "FCM - In vitro" = expression(paste("FCM - ", italic("In vitro"))), "FCM - In silico" = expression(paste("FCM - ", italic("In silico")))))
print(plot_mocks_12k_2)


### 5.1.2. Co-cultures ----
all_data2_cocult <- subset(all_data2_copycorrected_melted, Type == "Co-culture")
all_data2_cocult$ID <- as.factor(all_data2_cocult$ID)
order_IDs <- c("Co1_A_24h", "Co1_B_24h", "Co1_C_24h",
               "Co1_A_48h", "Co1_B_48h", "Co1_C_48h",
               "Co2_A_24h", "Co2_B_24h", "Co2_C_24h",
               "Co2_A_48h", "Co2_B_48h", "Co2_C_48h",
               "Co3_A_24h", "Co3_B_24h", "Co3_C_24h",
               "Co3_A_48h", "Co3_B_48h", "Co3_C_48h")
all_data2_cocult$ID <- factor(all_data2_cocult$ID, levels = order_IDs)

plot_cocult <- ggplot(data = all_data2_cocult, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, ncol = 3,
             labeller = labeller(ID = c(Co1_A_24h = "Co-culture 1A 24h", Co1_B_24h = "Co-culture 1B 24h", Co1_C_24h = "Co-culture 1C 24h",
                                        Co1_A_48h = "Co-culture 1A 48h", Co1_B_48h = "Co-culture 1B 48h", Co1_C_48h = "Co-culture 1C 48h",
                                        Co2_A_24h = "Co-culture 2A 24h", Co2_B_24h = "Co-culture 2B 24h", Co2_C_24h = "Co-culture 2C 24h",
                                        Co2_A_48h = "Co-culture 2A 48h", Co2_B_48h = "Co-culture 2B 48h", Co2_C_48h = "Co-culture 2C 48h",
                                        Co3_A_24h = "Co-culture 3A 24h", Co3_B_24h = "Co-culture 3B 24h", Co3_C_24h = "Co-culture 3C 24h",
                                        Co3_A_48h = "Co-culture 3A 48h", Co3_B_48h = "Co-culture 3B 48h", Co3_C_48h = "Co-culture 3C 48h"))) +
  theme_bw() +
  labs(x = "Technique", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800", "Other" = "#E76BF3"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")), "Other" = "Other")) +
  scale_x_discrete(labels = c("FCM_Model_SoFn" = "FCM SoFn", "FCM_Model_SoFnPg" = "FCM SoFnPg", "FCM_Model_SoFnPgVp" = "FCM SoFnPgVp", "Illumina" = "Illumina", "qPCR" = "qPCR"))
print(plot_cocult)

all_data2_cocult_box <- all_data2_cocult
all_data2_cocult_box$Culture <- substr(all_data2_cocult_box$ID, 1, 3)
all_data2_cocult_box$Replicate <- substr(all_data2_cocult_box$ID, 5, 5)
all_data2_cocult_box$Timepoint <- substr(all_data2_cocult_box$ID, 7, 9)
all_data2_cocult_box <- all_data2_cocult_box[, c(1, 3:8)]
all_data2_cocult_box_clean <- all_data2_cocult_box[all_data2_cocult_box$Relative_abundance != 0, ]

plot_cocult_box <- ggplot(data = all_data2_cocult_box_clean, aes(x = Timepoint, y = Relative_abundance, fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Replicate, y = Relative_abundance), position = position_dodge(width = 0.8), size = 6) +
  facet_grid(Technique ~ Culture,
             labeller = labeller(Culture = c(Co1 = "Co-culture 1", Co2 = "Co-culture 2", Co3 = "Co-culture 3"),
                                 Technique = c(FCM_Model_SoFn = "FCM SoFn", FCM_Model_SoFnPg = "FCM SoFnPg", FCM_Model_SoFnPgVp = "FCM SoFnPgVp", qPCR = "qPCR", Illumina = "Illumina"))) +
  theme_bw() +
  labs(y = "Relative abudnance (%)", fill = "Strain") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800", "Other" = "#E76BF3"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula"))))
print(plot_cocult_box)

plot_cocult_box2 <- ggplot(data = all_data2_cocult_box_clean, aes(x = Timepoint, y = Relative_abundance, fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  #geom_text(aes(label = Replicate, y = Relative_abundance), position = position_dodge(width = 0.8), size = 5) +
  ggrepel::geom_text_repel(aes(label = Replicate, y = Relative_abundance), position = position_dodge(width = 0.8), size = 5) +
  facet_grid(Culture ~ Technique,
             labeller = labeller(Culture = c(Co1 = "Co-culture 1", Co2 = "Co-culture 2", Co3 = "Co-culture 3"),
                                 Technique = c(FCM_Model_SoFn = "FCM SoFn", FCM_Model_SoFnPg = "FCM SoFnPg", FCM_Model_SoFnPgVp = "FCM SoFnPgVp", qPCR = "qPCR", Illumina = "Illumina"))) +
  theme_bw() +
  labs(y = "Relative abundance (%)", fill = "Strain") +
  theme(axis.text.x = element_text(vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800", "Other" = "#E76BF3"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula"))))
print(plot_cocult_box2)

all_data2_cocult_box_clean2 <- all_data2_cocult_box[all_data2_cocult_box$Strain != "Other", ]
all_data2_cocult_box_clean2 <- all_data2_cocult_box_clean2 %>% 
  dplyr::filter(!(Strain == "Pg" & Culture == "Co1" & Technique %in% c("FCM_Model_SoFn", "Illumina", "qPCR"))) %>% 
  dplyr::filter(!(Strain == "Vp" & Culture == "Co1" & Technique %in% c("FCM_Model_SoFn", "FCM_Model_SoFnPg", "Illumina", "qPCR"))) %>% 
  dplyr::filter(!(Strain == "Vp" & Culture == "Co2" & Technique %in% c("FCM_Model_SoFn", "FCM_Model_SoFnPg", "Illumina", "qPCR"))) %>% 
  dplyr::filter(!(Strain %in% c("Pg", "Vp") & Technique == "FCM_Model_SoFn")) %>% 
  dplyr::filter(!(Strain %in% c("Vp") & Technique == "FCM_Model_SoFnPg"))

plot_cocult_box3 <- ggplot(data = all_data2_cocult_box_clean2, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_grid(Culture ~ Timepoint,
             labeller = labeller(Culture = c(Co1 = "Co-culture 1", Co2 = "Co-culture 2", Co3 = "Co-culture 3"))) +
  theme_bw() +
  labs(y = "Relative abundance (%)", fill = "Strain") +
  theme(axis.text.x = element_text(vjust = 0.5, angle = 90),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")))) +
  scale_x_discrete(labels = c("FCM_Model_SoFn" = "FCM SoFn", "FCM_Model_SoFnPg" = "FCM SoFnPg", "FCM_Model_SoFnPgVp" = "FCM SoFnPgVp", "Illumina" = "Illumina", "qPCR" = "qPCR"))
print(plot_cocult_box3)

plot_cocult_box4 <- ggplot(data = all_data2_cocult_box_clean2, aes(x = Timepoint, y = Relative_abundance, fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  ggrepel::geom_text_repel(aes(label = Replicate, y = Relative_abundance), position = position_dodge(width = 0.8), size = 5) +
  facet_grid(Culture ~ Technique,
             labeller = labeller(Culture = c(Co1 = "Co-culture 1", Co2 = "Co-culture 2", Co3 = "Co-culture 3"),
                                 Technique = c(FCM_Model_SoFn = "FCM SoFn", FCM_Model_SoFnPg = "FCM SoFnPg", FCM_Model_SoFnPgVp = "FCM SoFnPgVp", qPCR = "qPCR", Illumina = "Illumina"))) +
  theme_bw() +
  labs(y = "Relative abundance (%)", fill = "Strain") +
  theme(axis.text.x = element_text(vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula"))))
print(plot_cocult_box4)


### 5.1.3. FCM wrongly predicted strains ----
# Mocks
rel_mocks_FCM <- subset(all_data2_mocks, Technique == "FCM_Model_SoFn" | Technique == "FCM_Model_SoFnPg" | Technique == "FCM_Model_SoFnPgVp")
rel_mocks_FCM <- rel_mocks_FCM[rel_mocks_FCM$Strain != "Other", ]

rel_mocks_FCM_wider <- rel_mocks_FCM %>% 
  pivot_wider(names_from = Strain, values_from = Relative_abundance)

rel_mocks_FCM_M1M8M9 <- subset(rel_mocks_FCM_wider, ID == "Mock 1" | ID == "Mock 8" | ID == "Mock 9")
rel_mocks_FCM_M1M8M9$wrong <- rel_mocks_FCM_M1M8M9$Pg + rel_mocks_FCM_M1M8M9$Vp

rel_mocks_FCM_M2 <- subset(rel_mocks_FCM_wider, ID == "Mock 2")
rel_mocks_FCM_M2$wrong <- rel_mocks_FCM_M2$Vp

rel_mocks_FCM_M3 <- subset(rel_mocks_FCM_wider, ID == "Mock 3")
rel_mocks_FCM_M3$wrong <- 0

rel_mocks_FCM_wider2 <- rbind(rel_mocks_FCM_M1M8M9, rel_mocks_FCM_M2, rel_mocks_FCM_M3)

# Co-cultures
rel_cocult_FCM <- subset(all_data2_cocult_box, Technique == "FCM_Model_SoFn" | Technique == "FCM_Model_SoFnPg" | Technique == "FCM_Model_SoFnPgVp")
rel_cocult_FCM <- rel_cocult_FCM[rel_cocult_FCM$Strain != "Other", ]

rel_cocult_FCM_wider <- rel_cocult_FCM %>% 
  pivot_wider(names_from = Strain, values_from = Relative_abundance)

rel_cocult_FCM_Co1 <- subset(rel_cocult_FCM_wider, Culture == "Co1")
rel_cocult_FCM_Co1$wrong <- rel_cocult_FCM_Co1$Pg + rel_cocult_FCM_Co1$Vp

rel_cocult_FCM_Co2 <- subset(rel_cocult_FCM_wider, Culture == "Co2")
rel_cocult_FCM_Co2$wrong <- rel_cocult_FCM_Co2$Vp

rel_cocult_FCM_Co3 <- subset(rel_cocult_FCM_wider, Culture == "Co3")
rel_cocult_FCM_Co3$wrong <- 0

rel_cocult_FCM_wider2 <- rbind(rel_cocult_FCM_Co1, rel_cocult_FCM_Co2, rel_cocult_FCM_Co3)

rel_cocult_FCM_mean <- rel_cocult_FCM_wider2 %>% 
  dplyr::group_by(Culture, Technique, Timepoint) %>% 
  summarize(Fn_mean = mean(Fn, na.rm = TRUE),
            So_mean = mean(So, na.rm = TRUE),
            Pg_mean = mean(Pg, na.rm = TRUE),
            Vp_mean = mean(Vp, na.rm = TRUE),
            wrong_mean = mean(wrong, na.rm = TRUE),
            Fn_sd = sd(Fn, na.rm = TRUE),
            So_sd = sd(So, na.rm = TRUE),
            Pg_sd = sd(Pg, na.rm = TRUE),
            Vp_sd = sd(Vp, na.rm = TRUE),
            wrong_sd = sd(wrong, na.rm = TRUE))


## 5.2. Absolute abundance ----
### 5.2.1. Mocks ----
FCM_mocks <- FCM_mean[(FCM_mean$ID == "Mix1" & FCM_mean$Technique == "FCM_Model_SoFn") |
                        (FCM_mean$ID == "Mix2" & FCM_mean$Technique == "FCM_Model_SoFnPg") |
                        (FCM_mean$ID == "Mix3" & FCM_mean$Technique == "FCM_Model_SoFnPgVp") |
                        (FCM_mean$ID == "Mix8" & FCM_mean$Technique == "FCM_Model_SoFn") |
                        (FCM_mean$ID == "Mix9" & FCM_mean$Technique == "FCM_Model_SoFn"), ]
FCM_mocks <- FCM_mocks[, c("ID", "So", "Fn", "Pg", "Vp")]
FCM_mocks$Technique <- "FCM"

qPCR_mocks_means$Technique <- "qPCR"

theoretical_mocks2 <- theoretical_mocks_rect[theoretical_mocks_rect$Mix_ID == "Mix1" | theoretical_mocks_rect$Mix_ID == "Mix2" | theoretical_mocks_rect$Mix_ID == "Mix3" | theoretical_mocks_rect$Mix_ID == "Mix8" | theoretical_mocks_rect$Mix_ID == "Mix9", ]
theoretical_mocks2$Technique <- "Actual"
colnames(theoretical_mocks2)[colnames(theoretical_mocks2) == "Mix_ID"] <- "ID"
theoretical_mocks2 <- theoretical_mocks2[, colSums(theoretical_mocks2 != 0) > 0]

absolute_mocks <- dplyr::bind_rows(FCM_mocks, qPCR_mocks_means, theoretical_mocks2)
absolute_mocks$ID <- gsub("Mix", "Mock ", absolute_mocks$ID)
absolute_mocks_melted <- reshape2::melt(absolute_mocks, id.vars = c("ID", "Technique"), variable.name = c("Strain"), value.name = c("Concentration"))

plot_mocks_absolute <- ggplot(data = absolute_mocks_melted, aes(x = Technique, y = Concentration, fill = Strain))+
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, nrow = length(unique(ID))) +
  theme_bw() +
  labs(x = "Technique", y = "Concentration (cells/mL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")))) +
  scale_y_continuous(breaks = c(0, 500000000, 1000000000, 2000000000, 3000000000, 4000000000, 5000000000, 6000000000)) +
  scale_x_discrete(labels = c("Actual" = "Theoretical", "FCM" = "FCM", "qPCR" = "qPCR"))
print(plot_mocks_absolute)

absolute_mocks_FCM <- absolute_mocks_melted[absolute_mocks_melted$Technique == "FCM" | absolute_mocks_melted$Technique == "Actual", ]
absolute_mocks_FCM$Technique <- gsub("FCM", "Predicted FCM", absolute_mocks_FCM$Technique)
absolute_mocks_FCM$Technique <- gsub("Actual", "Calculated", absolute_mocks_FCM$Technique)

plot_mocks_absolute_FCM <- ggplot(data = absolute_mocks_FCM, aes(x = Technique, y = Concentration, fill = Strain))+
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, nrow = length(unique(ID))) +
  theme_bw() +
  labs(x = "Technique", y = "Concentration (cells/mL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")))) +
  scale_y_continuous(breaks = seq(0, 300000000, by = 50000000)) +
  scale_x_discrete(labels = c("Calculated" = "Theoretical", "Predicted FCM" = "FCM"))
print(plot_mocks_absolute_FCM)

plot_mocks_absolute_grid <- plot_grid(plot_mocks_absolute, plot_mocks_absolute_FCM, labels = c("A", "B"), label_size = 34, hjust = -0.15, ncol = 2, nrow = 1)
legend_absolute_mocks <- get_legend(plot_mocks_absolute +
                                      theme(legend.position = "bottom"))
plot_mocks_absolute_combined <- plot_grid(plot_mocks_absolute_grid, legend_absolute_mocks, nrow = 2, rel_heights = c(1, 0.05))
print(plot_mocks_absolute_combined)


### 5.2.2. Co-cultures ----
FCM_cocult <- FCM_mean[(FCM_mean$ID == "Co1_A_24h" & FCM_mean$Technique == "FCM_Model_SoFn") |
                         (FCM_mean$ID == "Co1_B_24h" & FCM_mean$Technique == "FCM_Model_SoFn") |
                         (FCM_mean$ID == "Co1_C_24h" & FCM_mean$Technique == "FCM_Model_SoFn") |
                         (FCM_mean$ID == "Co1_A_48h" & FCM_mean$Technique == "FCM_Model_SoFn") |
                         (FCM_mean$ID == "Co1_B_48h" & FCM_mean$Technique == "FCM_Model_SoFn") |
                         (FCM_mean$ID == "Co1_C_48h" & FCM_mean$Technique == "FCM_Model_SoFn") |
                         (FCM_mean$ID == "Co2_A_24h" & FCM_mean$Technique == "FCM_Model_SoFnPg") |
                         (FCM_mean$ID == "Co2_B_24h" & FCM_mean$Technique == "FCM_Model_SoFnPg") |
                         (FCM_mean$ID == "Co2_C_24h" & FCM_mean$Technique == "FCM_Model_SoFnPg") |
                         (FCM_mean$ID == "Co2_A_48h" & FCM_mean$Technique == "FCM_Model_SoFnPg") |
                         (FCM_mean$ID == "Co2_B_48h" & FCM_mean$Technique == "FCM_Model_SoFnPg") |
                         (FCM_mean$ID == "Co2_C_48h" & FCM_mean$Technique == "FCM_Model_SoFnPg") |
                         (FCM_mean$ID == "Co3_A_24h" & FCM_mean$Technique == "FCM_Model_SoFnPgVp") |
                         (FCM_mean$ID == "Co3_B_24h" & FCM_mean$Technique == "FCM_Model_SoFnPgVp") |
                         (FCM_mean$ID == "Co3_C_24h" & FCM_mean$Technique == "FCM_Model_SoFnPgVp") |
                         (FCM_mean$ID == "Co3_A_48h" & FCM_mean$Technique == "FCM_Model_SoFnPgVp") |
                         (FCM_mean$ID == "Co3_B_48h" & FCM_mean$Technique == "FCM_Model_SoFnPgVp") |
                         (FCM_mean$ID == "Co3_C_48h" & FCM_mean$Technique == "FCM_Model_SoFnPgVp"), ]
FCM_cocult <- FCM_cocult[, c("ID", "So", "Fn", "Pg", "Vp")]
FCM_cocult$Technique <- "FCM"

qPCR_cocult_means$Technique <- "qPCR"

absolute_cocult <- dplyr::bind_rows(FCM_cocult, qPCR_cocult_means)

absolute_cocult$ID <- as.factor(absolute_cocult$ID)
order_IDs <- c("Co1_A_24h", "Co1_B_24h", "Co1_C_24h",
               "Co1_A_48h", "Co1_B_48h", "Co1_C_48h",
               "Co2_A_24h", "Co2_B_24h", "Co2_C_24h",
               "Co2_A_48h", "Co2_B_48h", "Co2_C_48h",
               "Co3_A_24h", "Co3_B_24h", "Co3_C_24h",
               "Co3_A_48h", "Co3_B_48h", "Co3_C_48h")
absolute_cocult$ID <- factor(absolute_cocult$ID, levels = order_IDs)

absolute_cocult_melted <- reshape2::melt(absolute_cocult, id.vars = c("ID", "Technique"), variable.name = c("Strain"), value.name = c("Concentration"))

plot_cocult_absolute <- ggplot(data = absolute_cocult_melted, aes(x = Technique, y = Concentration, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, ncol = 6,
             labeller = labeller(ID = c(Co1_A_24h = "Co-culture 1A 24h", Co1_B_24h = "Co-culture 1B 24h", Co1_C_24h = "Co-culture 1C 24h",
                                        Co1_A_48h = "Co-culture 1A 48h", Co1_B_48h = "Co-culture 1B 48h", Co1_C_48h = "Co-culture 1C 48h",
                                        Co2_A_24h = "Co-culture 2A 24h", Co2_B_24h = "Co-culture 2B 24h", Co2_C_24h = "Co-culture 2C 24h",
                                        Co2_A_48h = "Co-culture 2A 48h", Co2_B_48h = "Co-culture 2B 48h", Co2_C_48h = "Co-culture 2C 48h",
                                        Co3_A_24h = "Co-culture 3A 24h", Co3_B_24h = "Co-culture 3B 24h", Co3_C_24h = "Co-culture 3C 24h",
                                        Co3_A_48h = "Co-culture 3A 48h", Co3_B_48h = "Co-culture 3B 48h", Co3_C_48h = "Co-culture 3C 48h"))) +
  theme_bw() +
  labs(x = "Technique", y = "Concentration (cells/mL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")))) +
  scale_y_continuous(breaks = c(0, 500000000, 1000000000, 2000000000, 3000000000, 4000000000, 5000000000, 6000000000))
print(plot_cocult_absolute)


absolute_cocult_box <- absolute_cocult
absolute_cocult_box$Culture <- substr(absolute_cocult_box$ID, 1, 3)
absolute_cocult_box$Replicate <- substr(absolute_cocult_box$ID, 5, 5)
absolute_cocult_box$Timepoint <- substr(absolute_cocult_box$ID, 7, 9)

absolute_cocult_box_melted <- reshape2::melt(absolute_cocult_box, id.vars = c("Replicate", "Timepoint", "Culture", "Technique", "ID"), variable.name = "Strain", value.name = "Concentration")
absolute_cocult_box_melted_clean <- absolute_cocult_box_melted[absolute_cocult_box_melted$Concentration != 0, ]

plot_cocult_absolute_box <- ggplot(data = absolute_cocult_box_melted_clean, aes(x = Timepoint, y = Concentration, fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_text(aes(label = Replicate, y = Concentration), position = position_dodge(width = 0.8), size = 6) +
  facet_grid(Technique ~ Culture, scales = "free_y",
             labeller = labeller(Culture = c(Co1 = "Co-culture 1", Co2 = "Co-culture 2", Co3 = "Co-culture 3"),
                                 Technique = c(FCM = "Flow cytometry", qPCR = "qPCR"))) +
  theme_bw() +
  labs(y = "Concentration (cells/mL)", fill = "Strain") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula"))))
print(plot_cocult_absolute_box)

plot_cocult_absolute_box2 <- ggplot(data = absolute_cocult_box_melted_clean, aes(x = Timepoint, y = Concentration, fill = Strain)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  #geom_text(aes(label = Replicate, y = Concentration), position = position_dodge(width = 0.8), size = 5) +
  ggrepel::geom_text_repel(aes(label = Replicate, y = Concentration), position = position_dodge(width = 0.8), size = 5) +
  scale_y_log10() +
  facet_grid(Culture ~ Technique,
             labeller = labeller(Culture = c(Co1 = "Co-culture 1", Co2 = "Co-culture 2", Co3 = "Co-culture 3"),
                                 Technique = c(FCM = "Flow cytometry", qPCR = "qPCR"))) +
  theme_bw() +
  labs(y = "Concentration (cells/mL)", fill = "Strain") +
  theme(axis.text.x = element_text(vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula"))))
print(plot_cocult_absolute_box2)


# FCM only --> look at growth (cell numbers)
absolute_cocult_FCM <- absolute_cocult_melted[absolute_cocult_melted$Technique == "FCM", ]
absolute_cocult_FCM$Timepoint <- as.factor(substr(as.character(absolute_cocult_FCM$ID), start = nchar(as.character(absolute_cocult_FCM$ID))-2, stop = nchar(as.character(absolute_cocult_FCM$ID))))
absolute_cocult_FCM$ID <- as.factor(substr(as.character(absolute_cocult_FCM$ID), start = 1, stop = nchar(as.character(absolute_cocult_FCM$ID))-4))
absolute_cocult_FCM$ID <- gsub("Co1_", "Co-culture 1", absolute_cocult_FCM$ID)
absolute_cocult_FCM$ID <- gsub("Co2_", "Co-culture 2", absolute_cocult_FCM$ID)
absolute_cocult_FCM$ID <- gsub("Co3_", "Co-culture 3", absolute_cocult_FCM$ID)

plot_cocult_absolute_FCM <- ggplot(data = absolute_cocult_FCM, aes(x = Timepoint, y = Concentration, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, ncol = 3) +
  theme_bw() +
  labs(x = "Timepoint", y = "Concentration (cells/mL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")))) +
  scale_y_continuous(breaks = seq(0, 500000000, by = 100000000))
print(plot_cocult_absolute_FCM)

# qPCR only --> look at growth (cell numbers)
absolute_cocult_qPCR <- absolute_cocult_melted[absolute_cocult_melted$Technique == "qPCR", ]
absolute_cocult_qPCR$Timepoint <- as.factor(substr(as.character(absolute_cocult_qPCR$ID), start = nchar(as.character(absolute_cocult_qPCR$ID))-2, stop = nchar(as.character(absolute_cocult_qPCR$ID))))
absolute_cocult_qPCR$ID <- as.factor(substr(as.character(absolute_cocult_qPCR$ID), start = 1, stop = nchar(as.character(absolute_cocult_qPCR$ID))-4))
absolute_cocult_qPCR$ID <- gsub("Co1_", "Co-culture 1", absolute_cocult_qPCR$ID)
absolute_cocult_qPCR$ID <- gsub("Co2_", "Co-culture 2", absolute_cocult_qPCR$ID)
absolute_cocult_qPCR$ID <- gsub("Co3_", "Co-culture 3", absolute_cocult_qPCR$ID)

plot_cocult_absolute_qPCR <- ggplot(data = absolute_cocult_qPCR, aes(x = Timepoint, y = Concentration, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, ncol = 3) +
  theme_bw() +
  labs(x = "Timepoint", y = "Concentration (cells/mL)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "none",
        axis.text = element_text(size = 18),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 22),
        strip.text = element_text(size = 18),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")))) +
  scale_y_continuous(breaks = seq(0, 6000000000, by = 1000000000))
print(plot_cocult_absolute_qPCR)

# Combined plot co-cultures
plot_cocult_absolute_grid <- plot_grid(plot_cocult_absolute, plot_cocult_absolute_FCM, labels = c("A", "B"), label_size = 34, hjust = -0.15, ncol = 2, nrow = 1, rel_widths = c(1, 0.6))
legend_absolute_cocult <- get_legend(plot_cocult_absolute +
                                      theme(legend.position = "bottom"))
plot_cocult_absolute_combined <- plot_grid(plot_cocult_absolute_grid, legend_absolute_cocult, nrow = 2, rel_heights = c(1, 0.05))
print(plot_cocult_absolute_combined)

plot_cocult_absolute_grid2 <- plot_grid(plot_cocult_absolute_qPCR, plot_cocult_absolute_FCM, labels = c("A", "B"), label_size = 34, hjust = -0.15, ncol = 2, nrow = 1)
legend_absolute_cocult2 <- get_legend(plot_cocult_absolute_qPCR +
                                       theme(legend.position = "bottom"))
plot_cocult_absolute_combined2 <- plot_grid(plot_cocult_absolute_grid2, legend_absolute_cocult2, nrow = 2, rel_heights = c(1, 0.05))
print(plot_cocult_absolute_combined2)

# qPCR for Pg --> might there be some necrotrophy happening?
qPCR_cocult_Pg <- qPCR_cocult_means[, c(1, 4)]
qPCR_cocult_Pg <- dplyr::filter(qPCR_cocult_Pg, Pg > 0)
qPCR_cocult_Pg$Cocult <- as.factor(substr(as.character(qPCR_cocult_Pg$ID), start = 1, stop = 3))
qPCR_cocult_Pg$Replicate <- as.factor(substr(as.character(qPCR_cocult_Pg$ID), start = 5, stop = 5))
qPCR_cocult_Pg$Timepoint <- as.factor(substr(as.character(qPCR_cocult_Pg$ID), start = nchar(as.character(qPCR_cocult_Pg$ID))-2, stop = nchar(as.character(qPCR_cocult_Pg$ID))))
qPCR_cocult_Pg$Cocult <- gsub("Co2", "Co-culture 2", qPCR_cocult_Pg$Cocult)
qPCR_cocult_Pg$Cocult <- gsub("Co3", "Co-culture 3", qPCR_cocult_Pg$Cocult)

plot_cocult_absolute_Pg <- ggplot(data = qPCR_cocult_Pg, aes(x = Timepoint, y = Pg, color = Replicate, group = Replicate)) +
  geom_point(size = 7, alpha = 0.6) +
  geom_line(linewidth = 2, alpha = 0.6) +
  labs(y = "Concentration (cells/mL)", x = NULL) +
  facet_wrap(~ Cocult, ncol = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text.align = 0,
        panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 5000000), breaks = c(0, 500000, 1000000, 2000000, 3000000, 4000000, 5000000)) +
  scale_color_manual(values = c("A" = "darkred", "B" = "blue3", "C" = "#A3A500"))
print(plot_cocult_absolute_Pg)


## 5.3. RMSE relative abundance ----
# RMSE 50k cells
RMSE_melted <- reshape2::melt(RMSE, id.vars = c("ID"), variable.name = c("Technique"), value.name = c("RMSE"))

plot_RMSE <- ggplot(data = RMSE_melted, aes(x = ID, y = RMSE, color = Technique)) +
  geom_point(size = 7, alpha = 0.6) +
  labs(x = "Mock", y = "RMSE", color = "Technique") +
  scale_color_manual(values = c("RMSE_SoFn"= "#FF6C00", "RMSE_SoFnPg"= "darkred", "RMSE_SoFnPgVp" = "red3", "RMSE_qPCR" = "blue3", "RMSE_illumina" = "#008080"),
                     labels = c("RMSE_SoFn"= "FCM SoFn", "RMSE_SoFnPg"= "FCM SoFnPg", "RMSE_SoFnPgVp" = "FCM SoFnPgVp", "RMSE_qPCR" = "qPCR", "RMSE_illumina" = "Illumina")) +
  paper_theme_fab +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.text.align = 0) +
  scale_x_discrete(labels = c("Mix1" = "Mock 1", "Mix2" = "Mock 2", "Mix3" = "Mock 3", "Mix8" = "Mock 8", "Mix9" = "Mock 9")) +
  guides(color = guide_legend(nrow = 3))
print(plot_RMSE)

# RMSE FCM 12k cells
RMSE_FCM_12k_melted <- reshape2::melt(RMSE_FCM_12k, id.vars = c("ID"), variable.name = c("Technique"), value.name = c("RMSE"))

plot_RMSE_FCM_12k <- ggplot(data = RMSE_FCM_12k_melted, aes(x = ID, y = RMSE, color = Technique)) +
  geom_point(size = 7, alpha = 0.6) +
  labs(x = "Mock", y = "RMSE", color = NULL) +
  scale_color_manual(values = c("RMSE_invitro"= "darkred", "RMSE_insilico" = "blue3"),
                     labels = c("RMSE_invitro"= expression(paste("FCM - ", italic("In vitro"))), "RMSE_insilico" = expression(paste("FCM - ", italic("In silico"))))) +
  paper_theme_fab +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.text.align = 0) +
  guides(color = guide_legend(nrow = 2))
print(plot_RMSE_FCM_12k)

# Combined figure rel abundance and RMSE FCM 12k cells
plot_FCM_12k_grid <- plot_grid(plot_mocks_12k_2, plot_RMSE_FCM_12k, labels = c("A", "B"), label_size = 34, hjust = -0.15, ncol = 2, nrow = 1, rel_widths =c(1, 0.42))
print(plot_FCM_12k_grid)

# Combined figure rel abundance and RMSE 50k cells
plot_50k_grid <- plot_grid(plot_mocks, plot_RMSE, labels = c("A", "B"), label_size = 34, hjust = -0.15, ncol = 2, nrow = 1, rel_widths =c(1, 0.6))
print(plot_50k_grid)


# 6. Statistical analysis cell concentrations co-cultures ----

## 6.1. qPCR ----
total_concentration_qPCR_cocult <- aggregate(Concentration ~ ID + Timepoint, data = absolute_cocult_qPCR, FUN = sum)
total_concentration_qPCR_cocult$ID <- substr(total_concentration_qPCR_cocult$ID, 1, nchar(total_concentration_qPCR_cocult$ID)-1)

# REMARK: smallest possible p-value in Wilcox test with 3 replicates is 0.1 --> do not use Wilcox test!
sample_names_cocult <- unique(total_concentration_qPCR_cocult$ID)
wilcox_qPCR_cocult <- lapply(sample_names_cocult, function(sample_name) {
  df_sample <- total_concentration_qPCR_cocult[total_concentration_qPCR_cocult$ID == sample_name, ]
  test_result <- wilcox.test(Concentration ~ Timepoint, data = df_sample, paired = TRUE, alternative = "two.sided")
  return(test_result)
})

sample_names_cocult <- unique(total_concentration_qPCR_cocult$ID)
ttest_qPCR_cocult <- lapply(sample_names_cocult, function(sample_name) {
  df_sample <- total_concentration_qPCR_cocult[total_concentration_qPCR_cocult$ID == sample_name, ]
  test_result <- t.test(Concentration ~ Timepoint, data = df_sample, paired = TRUE, alternative = "two.sided")
  return(test_result)
})

stat_qPCR_cocult_1 <- total_concentration_qPCR_cocult[total_concentration_qPCR_cocult$ID == "Co-culture 1", ]
stat_qPCR_cocult_1$ID <- rep(c("A", "B", "C"), 2)
stat_qPCR_cocult_1 <- pivot_wider(stat_qPCR_cocult_1, names_from = Timepoint, values_from = Concentration)
colnames(stat_qPCR_cocult_1) <- c("ID", "h24", "h48")
stat_qPCR_cocult_1$diff <- stat_qPCR_cocult_1$h48-stat_qPCR_cocult_1$h24

shapiro_qPCR_cocult_1 <- shapiro.test(stat_qPCR_cocult_1$diff)
ggqqplot(stat_qPCR_cocult_1$diff)
ggdensity(stat_qPCR_cocult_1$diff,
          main = "Density plot of cell concentration",
          xlab = "Difference in concentration")

stat_qPCR_cocult_2 <- total_concentration_qPCR_cocult[total_concentration_qPCR_cocult$ID == "Co-culture 2", ]
stat_qPCR_cocult_2$ID <- rep(c("A", "B", "C"), 2)
stat_qPCR_cocult_2 <- pivot_wider(stat_qPCR_cocult_2, names_from = Timepoint, values_from = Concentration)
colnames(stat_qPCR_cocult_2) <- c("ID", "h24", "h48")
stat_qPCR_cocult_2$diff <- stat_qPCR_cocult_2$h48-stat_qPCR_cocult_2$h24

shapiro_qPCR_cocult_2 <- shapiro.test(stat_qPCR_cocult_2$diff)
ggqqplot(stat_qPCR_cocult_2$diff)
ggdensity(stat_qPCR_cocult_2$diff,
          main = "Density plot of cell concentration",
          xlab = "Difference in concentration")

stat_qPCR_cocult_3 <- total_concentration_qPCR_cocult[total_concentration_qPCR_cocult$ID == "Co-culture 3", ]
stat_qPCR_cocult_3$ID <- rep(c("A", "B", "C"), 2)
stat_qPCR_cocult_3 <- pivot_wider(stat_qPCR_cocult_3, names_from = Timepoint, values_from = Concentration)
colnames(stat_qPCR_cocult_3) <- c("ID", "h24", "h48")
stat_qPCR_cocult_3$diff <- stat_qPCR_cocult_3$h48-stat_qPCR_cocult_3$h24

shapiro_qPCR_cocult_3 <- shapiro.test(stat_qPCR_cocult_3$diff)
ggqqplot(stat_qPCR_cocult_3$diff)
ggdensity(stat_qPCR_cocult_3$diff,
          main = "Density plot of cell concentration",
          xlab = "Difference in concentration")

## 6.2. FCM ----
total_concentration_FCM_cocult <- aggregate(Concentration ~ ID + Timepoint, data = absolute_cocult_FCM, FUN = sum)
total_concentration_FCM_cocult$ID <- substr(total_concentration_FCM_cocult$ID, 1, nchar(total_concentration_FCM_cocult$ID)-1)

# REMARK: smallest possible p-value in Wilcox test with 3 replicates is 0.1 --> do not use Wilcox test!
sample_names_cocult <- unique(total_concentration_FCM_cocult$ID)
wilcox_FCM_cocult <- lapply(sample_names_cocult, function(sample_name) {
  df_sample <- total_concentration_FCM_cocult[total_concentration_FCM_cocult$ID == sample_name, ]
  test_result <- wilcox.test(Concentration ~ Timepoint, data = df_sample, paired = TRUE, alternative = "two.sided")
  return(test_result)
})

sample_names_cocult <- unique(total_concentration_FCM_cocult$ID)
ttest_FCM_cocult <- lapply(sample_names_cocult, function(sample_name) {
  df_sample <- total_concentration_FCM_cocult[total_concentration_FCM_cocult$ID == sample_name, ]
  test_result <- t.test(Concentration ~ Timepoint, data = df_sample, paired = TRUE, alternative = "two.sided")
  return(test_result)
})

stat_FCM_cocult_1 <- total_concentration_FCM_cocult[total_concentration_FCM_cocult$ID == "Co-culture 1", ]
stat_FCM_cocult_1$ID <- rep(c("A", "B", "C"), 2)
stat_FCM_cocult_1 <- pivot_wider(stat_FCM_cocult_1, names_from = Timepoint, values_from = Concentration)
colnames(stat_FCM_cocult_1) <- c("ID", "h24", "h48")
stat_FCM_cocult_1$diff <- stat_FCM_cocult_1$h48-stat_FCM_cocult_1$h24

shapiro_FCM_cocult_1 <- shapiro.test(stat_FCM_cocult_1$diff)
ggqqplot(stat_FCM_cocult_1$diff)
ggdensity(stat_FCM_cocult_1$diff,
          main = "Density plot of cell concentration",
          xlab = "Difference in concentration")

stat_FCM_cocult_2 <- total_concentration_FCM_cocult[total_concentration_FCM_cocult$ID == "Co-culture 2", ]
stat_FCM_cocult_2$ID <- rep(c("A", "B", "C"), 2)
stat_FCM_cocult_2 <- pivot_wider(stat_FCM_cocult_2, names_from = Timepoint, values_from = Concentration)
colnames(stat_FCM_cocult_2) <- c("ID", "h24", "h48")
stat_FCM_cocult_2$diff <- stat_FCM_cocult_2$h48-stat_FCM_cocult_2$h24

shapiro_FCM_cocult_2 <- shapiro.test(stat_FCM_cocult_2$diff)
ggqqplot(stat_FCM_cocult_2$diff)
ggdensity(stat_FCM_cocult_2$diff,
          main = "Density plot of cell concentration",
          xlab = "Difference in concentration")

stat_FCM_cocult_3 <- total_concentration_FCM_cocult[total_concentration_FCM_cocult$ID == "Co-culture 3", ]
stat_FCM_cocult_3$ID <- rep(c("A", "B", "C"), 2)
stat_FCM_cocult_3 <- pivot_wider(stat_FCM_cocult_3, names_from = Timepoint, values_from = Concentration)
colnames(stat_FCM_cocult_3) <- c("ID", "h24", "h48")
stat_FCM_cocult_3$diff <- stat_FCM_cocult_3$h48-stat_FCM_cocult_3$h24

shapiro_FCM_cocult_3 <- shapiro.test(stat_FCM_cocult_3$diff)
ggqqplot(stat_FCM_cocult_3$diff)
ggdensity(stat_FCM_cocult_3$diff,
          main = "Density plot of cell concentration",
          xlab = "Difference in concentration")
