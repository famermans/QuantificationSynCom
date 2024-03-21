# Step 1: Loading packages ----
# Set working directory, remove objects from workspace
setwd("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/DADA2")
rm(list = ls())

# r packages 
library(reshape2)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(vegan)
library(tibble)
library(dplyr)
library(DESeq2)
library(scales)
library(cluster)
library(ComplexHeatmap)
library(ggrepel)
library(circlize)
library(plyr)
library(Hmisc)
library(corrplot)
library(iNEXT)
library(readxl)
library(RColorBrewer)
#library(pairwiseAdonis)
#library(mixOmics)
#library(ComplexHeatmap)

# Specifying data location 
input_path <- "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/DADA2/DADA2/"
output_path <- "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/DADA2/Results/"
#functions needed
#source("KDP_tax_summarize.R")
#source("filter_mixomics.R")
#source("KDP_univar_statistics.R")
#source("volcanokim.R")
# Color vectors 
sevencolors <- c("#9BD0E3", "#DFDA5F", "#9FDF9D", "#E3B1D2", "#E3B573", "#6EDFC4", "#B2E078")
eightcolors <- c("#a24f7e", "#004692", "#5e4fa2", "#3682BA", "#98D5A4", "#FEF0A7", "#EE6445", "#E5E4E2")
elevencolors <- c("#000000", "#CD6600", "#36648B", "#008B00", "#7A378B", "#838B83", "#8B0000", "#323C4D", "#A03D44", "#98B736",
                  "#C8ADA4")
elevencolors2 <- c("#a24f7e", "#004692", "#5e4fa2", "#3682BA", "#5CB7A9", "#98D5A4", "#FEF0A7", "#E5E4E2", "#FA9C58", "#EE6445",
                   "#9E0142")
fourteencolors <- c("#a24f7e", "#004692", "#5e4fa2", "#3682BA", "#5CB7A9", "#98D5A4", "#D0EC9C", "#F3FAAD", "#FEF0A7", "#FA9C58",
                    "#EE6445", "#9E0142", "#fabebe", "#E5E4E2")
thirtycolors <- c("#323F24", "#CB51D7", "#72E245", "#DE4F2D", "#81DDC7", "#6584C6", "#CFAE3C", "#914261", "#BFDC86", "#CDC6BE",
                  "#D3469A", "#3C315A", "#607B30", "#DD4469", "#5B9072", "#CA9FC7", "#7F592D", "#5A656D", "#D1867F", "#4C2426",
                  "#79B6CE", "#6E67D1", "#CBDC3F", "#63D98B", "#97362B", "#CBB37D", "#7E3F85", "#CF8338", "#5DA73A", "#E5E4E2")
thirtycolors2 <- c("#a24f7e", "#004692", "#5e4fa2", "#3682BA", "#5CB7A9", "#98D5A4", "#D0EC9C", "#F3FAAD", "#FEF0A7", "#FA9C58",
                   "#EE6445", "#9E0142", "#fabebe", "#bcf60c", "#f032e6", "#46f0f0", "#911eb4", "#f58231", "#4363d8", "#ffe119",
                   "#3cb44b", "#e6194b", "#CBDC3F", "#63D98B", "#97362B", "#CBB37D", "#7E3F85", "#CF8338", "#5DA73A", "#E5E4E2")
fourtycolors <- c("#D3469A", "#3C315A", "#607B30", "#DD4469", "#5B9072", "#CA9FC7", "#7F592D", "#5A656D", "#D1867F", "#4C2426",
                  "#79B6CE", "#6E67D1", "#CBDC3F", "#63D98B", "#97362B", "#CBB37D", "#7E3F85", "#CF8338", "#5DA73A", "#CB80D0",
                  "#7eb96b", "#9366e9", "#6eb729", "#b55edb", "#46c353", "#a92fa4", "#afc52e", "#3468e5", "#b4a823", "#5a4fc4",
                  "#82cf62", "#e357bc", "#4b9925", "#cf6bd6", "#369039", "#e53486", "#39c685", "#bc2c84", "#8eab32", "#6b4bb1")

rainbowcolors <- rainbow(n = 134)

donor_colors <- c("#a24f7e", "#004692", "#5e4fa2", "#3682BA", "#5CB7A9", "#98D5A4", "#FEF0A7", "#EE6445", "#9E0142")
type_colors <- c("seagreen3", "#9E0142", "steelblue1")
population_colors <- c("yellow4", "#004692", "steelblue1", "seagreen3", "#9E0142")

source(file = "/Projects1/Fabian/paper_theme_fab.R")


# Step 2: Loading data ----
# Loading seq data 
ASVtable <- read.csv2(paste(input_path, "/ASVtable.csv", sep = ""), sep = ",", row.names = 1)
ASVtable <- t(ASVtable)
Taxtable <- read.csv2(paste(input_path, "/ASVtax138.species.csv", sep = ""), sep = ",", row.names = 1)
ASVsampleinfo <- readxl::read_excel("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/DADA2/Metadata_LGC_NGS3818.xlsx", sheet = "Sheet1") %>% 
  as.data.frame()
rownames(ASVsampleinfo) <- ASVsampleinfo$Name_LGC

# Order factor levels
ASVsampleinfo$Sample <- factor(ASVsampleinfo$Sample, levels = c(unique(ASVsampleinfo$Sample)))
ASVsampleinfo$Replicate <- factor(ASVsampleinfo$Replicate, levels = c(unique(ASVsampleinfo$Replicate)))
ASVsampleinfo$Timepoint <- factor(ASVsampleinfo$Timepoint, levels = c(unique(ASVsampleinfo$Timepoint)))
ASVsampleinfo$Date <- factor(ASVsampleinfo$Date, levels = c(unique(ASVsampleinfo$Date)))

metatable <- ASVsampleinfo


# Step 3: Preprocessing data ----
# Combining ASV table with taxtable
ASVtax <- dplyr::select(Taxtable, -Species)
ASVtax <- cbind(ASVtax[colnames(ASVtable), ], t(ASVtable))

ASVtable_prop <- sweep(ASVtable, 1, rowSums(ASVtable), '/')
rowSums(ASVtable_prop) # Check if previous operation was executed correctly


# Samples      
## Aggregate at genus level 
genusASVtable <- aggregate(. ~ Genus, data = ASVtax[ , -c(1:5)], FUN = sum)
rownames(genusASVtable) <- genusASVtable[ , 1]
genusASVtable <- genusASVtable[ , -1]

genusASVtable_prop <- sweep(genusASVtable, 2, colSums(genusASVtable), '/')
colSums(genusASVtable_prop) # Check if previous operation was executed correctly

## Aggregate at family level 
familyASVtable <- aggregate(. ~ Family, data = ASVtax[ , -c(1:4, 6)], FUN = sum)
rownames(familyASVtable) <- familyASVtable[ , 1]
familyASVtable <- familyASVtable[ , -1]

familyASVtable_prop <- sweep(familyASVtable, 2, colSums(familyASVtable), '/')
colSums(familyASVtable_prop) # Check if previous operation was executed correctly

## Aggregate at Phylum level 
phylumASVtable <- aggregate(. ~ Phylum, data = ASVtax[ , -c(1,3:6)], FUN = sum)
rownames(phylumASVtable) <- phylumASVtable[ , 1]
phylumASVtable <- phylumASVtable[ , -1]

phylumASVtable_prop <- sweep(phylumASVtable, 2, colSums(phylumASVtable), '/')
colSums(phylumASVtable_prop) # Check if previous operation was executed correctly


# Step 4: Quality control -> Rarefaction curves! ----
# Quality Control samples 
# 
#   ## Taxa distribution
#     readnr <- rowSums(ASVtable_QC)
#     ASVnr <- rowSums(ASVtable_QC!=0)
#     genusnr <- colSums(genusASVtable_QC!=0)
#     phylumnr <- colSums(phylumASVtable_QC!=0)
#     
#     taxoverview_QC <- data.frame(Sample=rownames(ASVtable_QC),ASV=ASVnr,Genus=genusnr,Phylum=phylumnr,Reads=readnr)
#     taxoverview <- melt(taxoverview_QC,id.vars=c("Sample"))
#     
#     plottaxoverview <- ggplot(taxoverview,aes(x=Sample,y=as.numeric(value))) + geom_bar(aes(fill=factor(variable)),stat='identity',position="identity") + facet_grid(.~variable,scales="free",drop=TRUE) + coord_flip()
#     plottaxoverview <- plottaxoverview + scale_fill_manual(values=sevencolors)
#     plottaxoverview <- plottaxoverview + ylab("Number of unique taxa") + xlab("Sample") + theme(legend.title=element_blank(),legend.position="none")
#     plottaxoverview
#     
#  ## Composition at genus level
#   genustable_QC_bar <- decostand(subset(genusASVtable_QC_prop,select=-c(NP1, NP2)),MARGIN=2,"total")*100
#   selecttop20par <- rowSums(genustable_QC_bar)
#   genustable_QC_bar <- rownames_to_column(genustable_QC_bar)
#   genustable_QC_bar <- as.data.frame(genustable_QC_bar %>% top_n(29,selecttop20par))
#   other <- 100-colSums(genustable_QC_bar[,-1])
#   rownames(genustable_QC_bar) <- genustable_QC_bar$rowname
#   genustable_QC_bar <- genustable_QC_bar[,-1]
#   genustable_QC_bar <- rbind(genustable_QC_bar,other)
#   rownames(genustable_QC_bar)[nrow(genustable_QC_bar)] <- "Other"
# 
#   genustable_QC_bar <- as.data.frame(t(genustable_QC_bar))
#   genustable_QC_bar <- cbind(genustable_QC_bar,rownames=rownames(genustable_QC_bar))
#   genustable_QC_bar <- melt(genustable_QC_bar,id.vars=c('rownames'))
# 
#   p <- ggplot(data=genustable_QC_bar,aes(x=rownames,y=as.numeric(as.character(value))))
#   p <- p + geom_bar(stat="identity",aes(fill=variable))
#   p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + scale_fill_manual(values=fourtycolors,name="genus")
#   p <- p + coord_flip() + ylab("Relative abundance (%)") + xlab(NULL) + theme(legend.text=element_text(face='italic')) + guides(fill=guide_legend(ncol=4)) + theme(legend.justification = "center")
#   p
# 
# # Replicates 
#   
#   # Composition at genus level
#     genustable_rep_bar <- decostand(genusASVtable_rep_prop,MARGIN=2,"total")*100
#     selecttop20par <- rowSums(genustable_rep_bar)
#     genustable_rep_bar <- rownames_to_column(genustable_rep_bar)
#     genustable_rep_bar <- as.data.frame(genustable_rep_bar %>% top_n(29,selecttop20par))
#     other <- 100-colSums(genustable_rep_bar[,-1])
#     rownames(genustable_rep_bar) <- genustable_rep_bar$rowname
#     genustable_rep_bar <- genustable_rep_bar[,-1]
#     genustable_rep_bar <- rbind(genustable_rep_bar,other)
#     rownames(genustable_rep_bar)[nrow(genustable_rep_bar)] <- "Other"
#     
#     genustable_rep_bar <- as.data.frame(t(genustable_rep_bar))
#     genustable_rep_bar <- cbind(genustable_rep_bar,rownames=rownames(genustable_rep_bar))
#     genustable_rep_bar <- melt(genustable_rep_bar,id.vars=c('rownames'))
# 
#     p <- ggplot(data=genustable_rep_bar,aes(x=rownames,y=as.numeric(as.character(value))))
#     p <- p + geom_bar(stat="identity",aes(fill=variable)) 
#     p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + scale_fill_manual(values=fourtycolors,name="genus")
#     p <- p + coord_flip() + ylab("Relative abundance (%)") + xlab(NULL) + theme(legend.text=element_text(face='italic')) + guides(fill=guide_legend(ncol=3)) + theme(legend.justification = "center")
#     p
#     
#   # Composition at ASV level
#     ASVtable_rep_bar <- decostand(data.frame(t(ASVtable_rep_prop)),MARGIN=2,"total")*100
#     selecttop20par <- rowSums((ASVtable_rep_bar))
#    # ASVtable_rep_bar <- data.frame(t(ASVtable_rep_bar))
#     ASVtable_rep_bar <- rownames_to_column(ASVtable_rep_bar)
#     ASVtable_rep_bar <- as.data.frame(ASVtable_rep_bar %>% top_n(29,selecttop20par))
#     other <- 100-colSums(ASVtable_rep_bar[,-1])
#     rownames(ASVtable_rep_bar) <- ASVtable_rep_bar$rowname
#     ASVtable_rep_bar <- ASVtable_rep_bar[,-1]
#     ASVtable_rep_bar <- rbind(ASVtable_rep_bar,other)
#     rownames(ASVtable_rep_bar)[nrow(ASVtable_rep_bar)] <- "Other"
#     
#     ASVtable_rep_bar <- as.data.frame(t(ASVtable_rep_bar))
#     ASVtable_rep_bar <- cbind(ASVtable_rep_bar,rownames=rownames(ASVtable_rep_bar))
#     ASVtable_rep_bar <- melt(ASVtable_rep_bar,id.vars=c('rownames'))
# 
#     p <- ggplot(data=ASVtable_rep_bar,aes(x=rownames,y=as.numeric(as.character(value))))
#     p <- p + geom_bar(stat="identity",aes(fill=variable)) 
#     p <- p + theme(axis.text.x = element_text(angle = 0, hjust = 1)) + scale_fill_manual(values=fourtycolors,name="ASV")
#     p <- p + coord_flip() + ylab("Relative abundance (%)") + xlab(NULL) + theme(legend.text=element_text(face='italic')) + guides(fill=guide_legend(ncol=4)) + theme(legend.justification = "center")
#     p
# 
# # Samples 
# 
#   # Taxa distribution
#     readnr <- rowSums(ASVtable)
#     ASVnr <- rowSums(ASVtable!=0)
#     genusnr <- colSums(genusASVtable!=0)
#     #familynr <- colSums(tax_summarize(ASVtable,Taxtable,"Family")!=0)
#     #ordernr <- colSums(tax_summarize(ASVtable,Taxtable,"Order")!=0)
#     #classnr <- colSums(tax_summarize(ASVtable,Taxtable,"Class")!=0)
#     phylumnr <- colSums(phylumASVtable!=0)
#     
#     taxoverview_QC <- data.frame(Sample=rownames(ASVtable),ASV=ASVnr,Genus=genusnr,Phylum=phylumnr,Reads=readnr)
#     
#     taxoverview_QC  <- cbind(ASVsampleinfo,taxoverview_QC)
#     taxoverview <- melt(taxoverview_QC,id.vars=c(colnames(ASVsampleinfo), 'Sample'))
#     
#     plottaxoverview <- ggplot(taxoverview,aes(x=Sample,y=as.numeric(value))) + geom_bar(aes(fill=factor(variable)),stat='identity',position="identity") + facet_grid(Nutrient_condition~variable,scales="free",drop=TRUE) + coord_flip()
#     plottaxoverview <- plottaxoverview + scale_fill_manual(values=sevencolors)
#     plottaxoverview <- plottaxoverview + ylab("Number of unique taxa") + xlab("Sample") + theme(legend.title=element_blank(),legend.position="none")
#     plottaxoverview

# Rarefaction
set.seed <- 777
datarare <- t(ASVtable)
# Find minimal total number of reads of all samples
raremax <- min(colSums(datarare))
# Draw rarefaction curves
rarecurve(round(ASVtable), step = 20, sample = raremax, col = c("blue"), cex = 0.6)
rarecurvedata <- rarecurve(round(ASVtable), step = 20, sample = raremax, cex = 0.6, tidy = TRUE) # format data for use in ggplot

metatable2 <- metatable
colnames(metatable2)[colnames(metatable2) == "Name_LGC"] <- "Site"
colnames(metatable2)[colnames(metatable2) == "Sample"] <- "Sample_name"
rarecurvedata_merged <- merge(rarecurvedata, metatable2, by = "Site", all.x = TRUE)
colnames(rarecurvedata_merged)[colnames(rarecurvedata_merged) == "Sample"] <- "Sample_size"

p_rarecurve <- ggplot(rarecurvedata_merged, aes(x = Sample_size, y = Species))+
  geom_line(aes(color = Timepoint, group = Site))+
  theme_bw()+
  theme(legend.position = "right")

data_ends <- rarecurvedata_merged %>%
  group_by(Site) %>%
  dplyr::summarize(Sample = max(Sample_size), Species = max(Species))
colnames(data_ends)[colnames(data_ends) == "Sample"] <- "Sample_size"

data_ends2 <- dplyr::left_join(data_ends, metatable2, by = "Site")

p_rarecurve + geom_text_repel(aes(label = Sample_name), data = data_ends2,
  fontface ="plain", color = "black", size = 3) +
  labs(x = "Number of reads", y = "Number of ASVs")

# Select data Co2 as example
rarecurvedata_select <- subset(rarecurvedata_merged, Sample_name %in% c('Co2'))

p_rarecurve_select <- ggplot(rarecurvedata_select, aes(x = Sample_size, y = Species))+
  geom_line(aes(color = Timepoint, group = Site))+
  paper_theme_fab+
  theme(legend.position = "right")

data_ends_select <- rarecurvedata_select %>% 
  group_by(Site) %>% 
  dplyr::summarize(Sample = max(Sample_size), Species = max(Species))
data_ends_select <- dplyr::left_join(data_ends_select, metatable2, by = "Site")
colnames(data_ends_select)[colnames(data_ends_select) == "Sample"] <- "Sample_size"

p_rarecurve_select + geom_text_repel(aes(label = Replicate), data = data_ends_select,
                                     fontface = "plain", color = "black", size = 6)+
  labs(x = "Number of reads", y = "Number of ASVs")


# Quality reads
quality_scores_Fs <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/quality_scores_Fs.rds")
quality_scores_Rs <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/quality_scores_Rs.rds")
quality_scores_Fs_filt <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/quality_scores_Fs_filt.rds")
quality_scores_Rs_filt <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/quality_scores_Rs_filt.rds")
mean_quality_Fs <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/mean_quality_Fs.rds")
mean_quality_Rs <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/mean_quality_Rs.rds")
mean_quality_Fs_filt <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/mean_quality_Fs_filt.rds")
mean_quality_Rs_filt <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/mean_quality_Rs_filt.rds")
area_quality_filt <- readRDS(file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/area_quality_filt.rds")

rownames(area_quality_filt) <- area_quality_filt$Sample
area_quality_meta <- merge(area_quality_filt, ASVsampleinfo, by = "row.names")
area_quality_meta <- area_quality_meta[, -c(1:2)]

area_quality_meta_melted <- area_quality_meta[, c(1:3, 5:7, 9)] %>% 
  reshape2::melt(id.vars = c("Name_LGC", "Sample.y", "Replicate", "Timepoint", "Type"), value.name = "Area")

plot_quality_reads <- ggplot(area_quality_meta_melted, aes(x = Sample.y, y = Area, color = variable, shape = Replicate, size = Timepoint))+
  geom_point()+
  labs(y = "Area under mean quality", x = NULL, color = "Direction")+
  paper_theme_fab+
  scale_color_manual(values = c("darkblue", "darkorange"), labels = c("Forward", "Reverse"))+
  scale_size_manual(values = c(7, 5, 3)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
print(plot_quality_reads)



# Step 5: Normalisation ----
# Proportions
ASVtable_proportions <- decostand(ASVtable, MARGIN = 1, "total")*100
rowSums(ASVtable_proportions)

# rlog 
ASVtable_rlog <- rlogTransformation(as.matrix(t(ASVtable) + 1), blind = TRUE)

# Variance Stabilizing Transformation
ASVtable_vst <- varianceStabilizingTransformation(as.matrix(t(ASVtable) + 1), blind = TRUE, fitType = "parametric")

# size factor normalisation
# deseqdata <- DESeqDataSetFromMatrix(t(ASVtable), colData = ASVsampleinfo, design = ~ NC_PD)
# deseqdata <- estimateSizeFactors(deseqdata)
# ASVtable_sf <- counts(deseqdata, normalized = TRUE)



# Step 6: Visualizing proportions ----
# Phylum level
phylumtable_bar <- decostand(phylumASVtable_prop, MARGIN = 2, "total")*100
selecttop20bar <- rowSums(phylumtable_bar)
phylumtable_bar <- rownames_to_column(phylumtable_bar)
phylumtable_bar <- as.data.frame(phylumtable_bar %>% top_n(7, selecttop20bar))
other <- 100-colSums(phylumtable_bar[, -1])
rownames(phylumtable_bar) <- phylumtable_bar$rowname
phylumtable_bar <- phylumtable_bar[, -1]
phylumtable_bar <- rbind(phylumtable_bar, other)
rownames(phylumtable_bar)[nrow(phylumtable_bar)] <- "Other"

phylumtable_bar <- cbind(ASVsampleinfo, as.data.frame(t(phylumtable_bar)))
#phylumtable_bar <- subset(phylumtable_bar, Type == "Mock")
phylumtable_bar <- melt(phylumtable_bar, id.vars = colnames(ASVsampleinfo))
phylumtable_bar$visual <- paste(phylumtable_bar$Sample, phylumtable_bar$Replicate, phylumtable_bar$Timepoint, sep = "_")

p <- ggplot(data = phylumtable_bar, aes(x = visual, y = value, fill = variable)) 
p <- p + geom_bar(position = "stack", stat = "identity", color = "black") + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0, hjust = 1)) + facet_grid(~ Type) + scale_fill_manual(name = 'Phylum', values = eightcolors)
p <- p + ylab("Relative abundance (%)") + xlab("Sample") + theme(legend.text = element_text(face = 'italic')) + guides(fill = guide_legend(ncol = 1))
p


# Genus level
genustable_bar <- decostand(genusASVtable_prop, MARGIN = 2, "total")*100
selecttop20bar <- rowSums(genustable_bar)
genustable_bar <- rownames_to_column(genustable_bar)
genustable_bar <- as.data.frame(genustable_bar %>% top_n(29, selecttop20bar))
other <- 100-colSums(genustable_bar[, -1])
rownames(genustable_bar) <- genustable_bar$rowname
genustable_bar <- genustable_bar[, -1]
genustable_bar <- rbind(genustable_bar, other)
rownames(genustable_bar)[nrow(genustable_bar)] <- "Other"

genustable_bar  <- cbind(ASVsampleinfo, as.data.frame(t(genustable_bar)))
#genustable_bar  <- subset(genustable_bar, Type == "Co-culture" | Type == "Mock")
genustable_bar  <- melt(genustable_bar, id.vars = colnames(ASVsampleinfo))
genustable_bar$visual <- paste(genustable_bar$Sample, genustable_bar$Replicate, genustable_bar$Timepoint, sep = "_")

p <- ggplot(data = genustable_bar, aes(x = visual, y = value, fill = variable))
p <- p + geom_bar(position = "stack", stat = "identity", color = "black") + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 90), strip.text.y = element_text(angle = 0, hjust = 1)) + scale_fill_manual(name = 'Genus', values = thirtycolors)
p <- p + ylab("Relative abundance (%)") + xlab("Population") + theme(legend.text = element_text(face = 'italic')) + guides(fill = guide_legend(ncol = 2)) + coord_flip()
p


genustable_bar_top13 <- decostand(genusASVtable_prop, MARGIN = 2, "total")*100
selecttop13bar <- rowSums(genustable_bar_top13)
genustable_bar_top13 <- rownames_to_column(genustable_bar_top13)
genustable_bar_top13 <- as.data.frame(genustable_bar_top13 %>% top_n(13, selecttop13bar))
other_top13 <- 100-colSums(genustable_bar_top13[, -1])
rownames(genustable_bar_top13) <- genustable_bar_top13$rowname
genustable_bar_top13 <- genustable_bar_top13[, -1]
genustable_bar_top13 <- rbind(genustable_bar_top13, other_top13)
rownames(genustable_bar_top13)[nrow(genustable_bar_top13)] <- "Other"

saveRDS(object = genustable_bar_top13, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/genustable_top13.rds")

genustable_bar_top13 <- t(genustable_bar_top13)
genustable_bar_top13 <- merge(genustable_bar_top13, ASVsampleinfo, by = "row.names")
genustable_bar_top13 <- genustable_bar_top13[, -1]
#genustable_bar_top13 <- subset(genustable_bar_top13, Type == "Co-culture" | Type == "Mock")
genustable_bar_top13  <- melt(genustable_bar_top13, id.vars = colnames(ASVsampleinfo))
genustable_bar_top13$visual <- paste(genustable_bar_top13$Sample, genustable_bar_top13$Replicate, genustable_bar_top13$Timepoint, sep = "_")

saveRDS(object = genustable_bar_top13, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/genustable_top13_melted.rds")

p <- ggplot(data = genustable_bar_top13, aes(y = visual, x = value, fill = variable)) 
p <- p + geom_bar(position = "stack", stat = "identity", color = "black") + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 0), strip.text.y = element_text(angle = 0, hjust = 1)) + scale_fill_manual(name = 'Genus', values = fourteencolors)
p <- p + xlab("Relative abundance (%)") + ylab("Population") + theme(legend.text = element_text(face = 'italic')) + guides(fill = guide_legend(ncol = 2))
p


# ASV level 
ASVtable_bar <- decostand(data.frame(t(ASVtable)), MARGIN = 2, "total")*100
selecttop20bar <- rowSums(ASVtable_bar)
ASVtable_bar <- rownames_to_column(ASVtable_bar)
ASVtable_bar <- as.data.frame(ASVtable_bar %>% top_n(29, selecttop20bar))
other <- 100-colSums(ASVtable_bar[, -1])
rownames(ASVtable_bar) <- ASVtable_bar$rowname
ASVtable_bar <- ASVtable_bar[, -1]
ASVtable_bar <- rbind(ASVtable_bar, other)
rownames(ASVtable_bar)[nrow(ASVtable_bar)] <- "Other"

ASVtable_bar <- t(ASVtable_bar)
ASVtable_bar <- merge(ASVtable_bar, ASVsampleinfo, by = "row.names")
ASVtable_bar <- ASVtable_bar[, -1]
ASVtable_bar <- melt(ASVtable_bar, id.vars = colnames(ASVsampleinfo))
ASVtable_bar$visual <- paste(ASVtable_bar$Sample, ASVtable_bar$Replicate, ASVtable_bar$Timepoint, sep = "_")

p <- ggplot(data = ASVtable_bar, aes(x = visual, y = value, fill = variable)) 
p <- p + geom_bar(position = 'stack', stat = 'identity', color = "black") + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 0), strip.text.y = element_text(angle = 0, hjust = 1)) + scale_fill_manual(name = 'ASV',values = thirtycolors2)
p <- p + ylab("Relative abundance (%)") + xlab("Population") + theme(legend.text = element_text(face = 'italic')) + guides(fill = guide_legend(ncol = 2)) + coord_flip()
p

