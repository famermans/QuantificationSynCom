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

qPCR_cocult_means_prop <- sweep(as.matrix(qPCR_cocult_means[, -c(1)]), 1, rowSums(qPCR_cocult_means[, -c(1)]), FUN = "/") %>% 
  as.data.frame()

qPCR_cocult_means_prop$ID <- qPCR_cocult_means$ID
qPCR_cocult_means_prop$Type <- "Co-culture"
qPCR_cocult_means_prop$Technique <- "qPCR"
qPCR_cocult_means_prop$Replicate <- rep(c("A", "B", "C"), 6)
qPCR_cocult_means_prop$Timepoint <- c(rep("24h", 9), rep("48h", 9))
qPCR_cocult_means_prop$ID <- paste(qPCR_cocult_means_prop$ID, qPCR_cocult_means_prop$Replicate, qPCR_cocult_means_prop$Timepoint, sep = "_")
qPCR_cocult_means_prop <- qPCR_cocult_means_prop[, c(1:7)]

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
FCM_SoFn_mean_rect$Type <- c(rep("Co-culture", 12), rep("Mock", 5))
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

ID_models <- data.frame(ID = c("Mix1", "Mix2", "Mix3", "Mix4", "Mix5", "Mix6", "Mix7", "Mix8", "Mix9"),
                        Model = c("SoFn", "SoFnPg", "SoFnPgVp", "AnAvSgSmiSoSsalSsanVp", "SgSmiSmuSoSsalSsanSsob", "AaFnPgPiSmuSsob", "AaAnAvFnPgPiSgSmiSmuSoSsalSsanSsobVp", "SoFn", "SoFn"))

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
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text.align = 0) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#A3A500", "So" = "#00BF7D", "Vp" = "#00B0F6", "Other" = "#E76BF3"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")), "Other" = "Other")) +
  scale_x_discrete(labels = c("Actual" = "Actual", "FCM_Model_SoFn" = "FCM SoFn", "FCM_Model_SoFnPg" = "FCM SoFnPg", "FCM_Model_SoFnPgVp" = "FCM SoFnPgVp", "Illumina" = "Illumina", "qPCR" = "qPCR"))
print(plot_mocks)

mocks_12k_melted <- reshape2::melt(mocks_12k, id.vars = c("ID", "Technique"), variable.name = c("Strain"), value.name = c("Relative_abundance"))
mocks_12k_melted$ID <- gsub("Mix1", "Mock 1", mocks_12k_melted$ID)
mocks_12k_melted$ID <- gsub("Mix2", "Mock 2", mocks_12k_melted$ID)
mocks_12k_melted$ID <- gsub("Mix3", "Mock 3", mocks_12k_melted$ID)
mocks_12k_melted$ID <- gsub("Mix4", "Mock 4", mocks_12k_melted$ID)
mocks_12k_melted$ID <- gsub("Mix5", "Mock 5", mocks_12k_melted$ID)
mocks_12k_melted$ID <- gsub("Mix6", "Mock 6", mocks_12k_melted$ID)
mocks_12k_melted$ID <- gsub("Mix7", "Mock 7", mocks_12k_melted$ID)
mocks_12k_melted$ID <- gsub("Mix8", "Mock 8", mocks_12k_melted$ID)
mocks_12k_melted$ID <- gsub("Mix9", "Mock 9", mocks_12k_melted$ID)
mocks_12k_melted$Technique <- gsub("Theoretical", "Actual", mocks_12k_melted$Technique)

plot_mocks_12k <- ggplot(data = mocks_12k_melted, aes(x = Technique, y = Relative_abundance, fill = Strain)) +
  geom_bar(position = "stack", stat = "identity", color = "black") +
  facet_wrap(~ ID, nrow = length(unique(ID))) +
  theme_bw() +
  labs(x = "Sample", y = "Relative abundance (%)") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.position = "bottom",
        axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        strip.text = element_text(size = 14),
        legend.text.align = 0) +
  scale_fill_manual(values = c("Fn" = "#F8766D", "Pg" = "#E38900", "So" = "#C49A00", "Vp" = "#99A800", "An" = "#53B400", "Av" = "#00BC56", "Aa" = "#00C094", "Sg" = "#00BFC4", "Smi" = "#00B6EB", "Smu" = "#06A4FF", "Pi" = "#A58AFF", "Ssal" = "#DF70F8", "Ssan" = "#FB61D7", "Ssob" = "#FF66A8"),
                    labels = c("Fn" = expression(italic("F. nucleatum")), "Pg" = expression(italic("P. gingivalis")), "So" = expression(italic("S. oralis")), "Vp" = expression(italic("V. parvula")), "An" = expression(italic("A. naeslundii")), "Av" = expression(italic("A. viscosus")), "Aa" = expression(italic("A. actinomycetemcomitans")), "Sg" = expression(italic("S. gordonii")), "Smi" = expression(italic("S. mitis")), "Smu" = expression(italic("S. mutans")), "Pi" = expression(italic("P. intermedia")), "Ssal" = expression(italic("S. salivarius")), "Ssan" = expression(italic("S. sanguinis")), "Ssob" = expression(italic("S. sobrinus")))) +
  scale_x_discrete(labels = c("Actual" = "Actual", "FCM" = "FCM"))
print(plot_mocks_12k)


## 5.2. RMSE ----
RMSE_melted <- reshape2::melt(RMSE, id.vars = c("ID"), variable.name = c("Technique"), value.name = c("RMSE"))

plot_RMSE <- ggplot(data = RMSE_melted, aes(x = ID, y = RMSE, color = Technique)) +
  geom_point(size = 5, alpha = 0.6) +
  labs(x = "Mock", y = "RMSE", color = "Technique") +
  scale_color_manual(values = c("RMSE_SoFn"= "#FF6C00", "RMSE_SoFnPg"= "darkred", "RMSE_SoFnPgVp" = "red3", "RMSE_qPCR" = "blue3", "RMSE_illumina" = "#008080"),
                     labels = c("RMSE_SoFn"= "FCM SoFn", "RMSE_SoFnPg"= "FCM SoFnPg", "RMSE_SoFnPgVp" = "FCM SoFnPgVp", "RMSE_qPCR" = "qPCR", "RMSE_illumina" = "Illumina")) +
  paper_theme_fab +
  theme(legend.position = "right",
        axis.title.x = element_blank()) +
  scale_x_discrete(labels = c("Mix1" = "Mock 1", "Mix2" = "Mock 2", "Mix3" = "Mock 3", "Mix8" = "Mock 8", "Mix9" = "Mock 9"))
print(plot_RMSE)
