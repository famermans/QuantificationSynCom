# Step 1: Load libraries, set settings ----
# Clear environment and set working directory
rm(list = ls())
setwd("/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/DADA2")

# Load libraries
library("dada2")
library("ggplot2")
library("lattice")
library("stringr")
library("knitr")
library("dplyr")

if(!dir.exists(file.path("DADA2"))){
  dir.create(file.path("DADA2"), recursive = TRUE)}
outputloc <- file.path("DADA2")

species <- "yes" # Need species assignment?
database <- "silva" # Options: silva, unite or midas
#setDadaOpt(OMEGA_A = 10^-40, OMEGA_C = 10^-40, OMEGA_P = 10^-4, DETECT_SINGLETONS = FALSE) # Adapt if global DADA options are to be changed


# Step 2: Loading data ----
# Input data: paired-end fastq files that have been split (or “demultiplexed”) by sample and from which the primers, barcodes/adapters have already been removed

fldr <- "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/DADA2/SeqData/PrimerClipped_fastq" # Manual directory override
path <- fldr # Set the path to the data folder.
list.files(path)

system(paste("cd ", fldr, "&& bunzip2 *", sep = "")) # Unzip
system(paste("cd ", fldr, "&& ls *.fastq | xargs -I {} rename 's/-/_/g' {}", sep = "")) # Remove dashes and replace to underscores

# Read in files: Forward and reverse fastq filenames have format: SAMPLENAME_R1.fastq and SAMPLENAME_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnFs
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
fnRs

# Sample names require some messing around with

# Extract sample names, assuming filenames have format: SAMPLENAME_R1/2.fastq
sample.names <- basename(fnFs) # load forward sequence names
#sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 3) # only retains the 3th value separated by underscores
sample.names <- gsub("341F_785R_P.._..._CMET._", "", sample.names)  # Used to remove artifacts from more complicated strings
sample.names <- gsub("341F_785R_P.._..._CMETEX_", "", sample.names)
sample.names <- gsub("_R..fastq", "", sample.names) # Used to remove artifacts from more complicated strings

sample.names

# Create contigs
dacontigs.df <- do.call(rbind, Map(data.frame, fnFs = basename(fnFs), fnRs = basename(fnRs)))
rownames (dacontigs.df) <- sample.names
save(dacontigs.df, file = "DADA2/dacontigs.Rda")

kable(dacontigs.df)


# Step 3: Inspecting read quality ----
## Read quality
# Plot the quality profile of the reads. 
dada2::plotQualityProfile(fnFs[1:5]) + geom_vline(aes(xintercept = 240), linetype = "dotted")
dada2::plotQualityProfile(fnRs[1:5]) + geom_vline(aes(xintercept = 240), linetype = "dotted")

## Extract read quality
quality_scores_Fs <- list()
mean_quality_Fs <- list()

for (i in 1:length(fnFs)) {
  quality_scores_Fs[[i]] <- dada2::plotQualityProfile(fnFs[i])
  
  mean_quality_Fs[[i]] <- list(
    data = as.data.frame(quality_scores_Fs[[i]]$plot_env$means),
    name = gsub(".*_FM", "FM", quality_scores_Fs[[i]]$data$file[1])
  )
  
  colnames(mean_quality_Fs[[i]]$data)[colnames(mean_quality_Fs[[i]]$data) == "V1"] <- "Mean_quality_score"
  mean_quality_Fs[[i]]$data$Cycle <- rownames(mean_quality_Fs[[i]]$data)
  mean_quality_Fs[[i]]$name <- gsub(".fastq", "", mean_quality_Fs[[i]]$name)
}

quality_scores_Rs <- list()
mean_quality_Rs <- list()

for (i in 1:length(fnRs)) {
  quality_scores_Rs[[i]] <- dada2::plotQualityProfile(fnRs[i])
  
  mean_quality_Rs[[i]] <- list(
    data = as.data.frame(quality_scores_Rs[[i]]$plot_env$means),
    name = gsub(".*_FM", "FM", quality_scores_Rs[[i]]$data$file[1])
  )
  
  colnames(mean_quality_Rs[[i]]$data)[colnames(mean_quality_Rs[[i]]$data) == "V1"] <- "Mean_quality_score"
  mean_quality_Rs[[i]]$data$Cycle <- rownames(mean_quality_Rs[[i]]$data)
  mean_quality_Rs[[i]]$name <- gsub(".fastq", "", mean_quality_Rs[[i]]$name)
}

saveRDS(object = mean_quality_Fs, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/mean_quality_Fs.rds")
saveRDS(object = mean_quality_Rs, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/mean_quality_Rs.rds")

saveRDS(object = quality_scores_Fs, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/quality_scores_Fs.rds")
saveRDS(object = quality_scores_Rs, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/quality_scores_Rs.rds")


## Read cleanup
# Reads were trimmed based on quality plots and trunQ quality score cut-off. In addition reads with ambiguous bases (maxN = 0) and reads with more than 2 maxEE (expected errors) were filtered out. The result of this filter is plotted below:

# Step3: Remove primers, filter and trim reads
# Filtering 

## Place filtered files in filtered/ subdirectory and because the compress option will be used in the filtering step below, file names will have the fastq.gz extension
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

## Locate the primers in the files in notepad ++. Primers are not retrieved => they have been trimmed already 

#341F: CCTACGGGNGGCWGCAG =>  CCTACGGG[ACGT]GGC[AT]GCAG
#785R: GACTACHVGGGTATCTAAKCC =>  GACTAC[ACT][ACG]GGGTATCTAA[GT]CC => 785R rev comp: GG[CA]TTAGATACCC[TGC][TGA]GTAGTC
#FWD_PRIMER_LEN <- 19 
#REV_PRIMER_LEN <- 20

## Remove primers, trim read lenght based on quality plots and trunQ quality score cut-off and filter out reads with ambiguous bases (maxN=0), reads with higher than maxEE "expected errors" and reads that match against the phiX genome. The trunQ value comes from The Illumina manual (page 30): If a read ends with a segment of mostly low quality (Q15 or below), then all of the quality values in the segment are replaced with a value of 2 (encoded as the letter B in Illumina's text-based encoding of quality scores)... This Q2 indicator does not predict a specific error rate, but rather indicates that a specific final portion of the read should not be used in further analyses.
out <- filterAndTrim(fnFs,filtFs,fnRs,filtRs,trimLeft=c(0,0),truncLen=c(240,240),maxN=0,maxEE=c(2,2),truncQ=2,rm.phix=TRUE,compress=TRUE,multithread=TRUE) 
head(out)

exists <- file.exists(filtFs) & file.exists(filtRs)
filtFs <- filtFs[exists]
filtRs <- filtRs[exists]

save(filtFs,file="DADA2/filtFs.Rda")
save(filtRs,file="DADA2/filtRs.Rda")


# Step 3b ----
load(file="DADA2/filtFs.Rda")
load(file="DADA2/filtRs.Rda")
## Plot the quality profile of the reads. 
dada2::plotQualityProfile(filtFs[1:5])
dada2::plotQualityProfile(filtRs[1:5])

## Extract read quality
quality_scores_Fs_filt <- list()
mean_quality_Fs_filt <- list()

for (i in 1:length(filtFs)) {
  quality_scores_Fs_filt[[i]] <- dada2::plotQualityProfile(filtFs[i])
  
  mean_quality_Fs_filt[[i]] <- list(
    data = as.data.frame(quality_scores_Fs_filt[[i]]$plot_env$means),
    name = gsub(".*_FM", "FM", quality_scores_Fs_filt[[i]]$data$file[1])
  )
  
  colnames(mean_quality_Fs_filt[[i]]$data)[colnames(mean_quality_Fs_filt[[i]]$data) == "V1"] <- "Mean_quality_score"
  mean_quality_Fs_filt[[i]]$data$Cycle <- rownames(mean_quality_Fs_filt[[i]]$data)
  mean_quality_Fs_filt[[i]]$name <- gsub("_filt.fastq.gz", "", mean_quality_Fs_filt[[i]]$name)
}

quality_scores_Rs_filt <- list()
mean_quality_Rs_filt <- list()

for (i in 1:length(filtRs)) {
  quality_scores_Rs_filt[[i]] <- dada2::plotQualityProfile(filtRs[i])
  
  mean_quality_Rs_filt[[i]] <- list(
    data = as.data.frame(quality_scores_Rs_filt[[i]]$plot_env$means),
    name = gsub(".*_FM", "FM", quality_scores_Rs_filt[[i]]$data$file[1])
  )
  
  colnames(mean_quality_Rs_filt[[i]]$data)[colnames(mean_quality_Rs_filt[[i]]$data) == "V1"] <- "Mean_quality_score"
  mean_quality_Rs_filt[[i]]$data$Cycle <- rownames(mean_quality_Rs_filt[[i]]$data)
  mean_quality_Rs_filt[[i]]$name <- gsub("_filt.fastq.gz", "", mean_quality_Rs_filt[[i]]$name)
}

saveRDS(object = mean_quality_Fs_filt, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/mean_quality_Fs_filt.rds")
saveRDS(object = mean_quality_Rs_filt, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/mean_quality_Rs_filt.rds")

saveRDS(object = quality_scores_Fs_filt, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/quality_scores_Fs_filt.rds")
saveRDS(object = quality_scores_Rs_filt, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/quality_scores_Rs_filt.rds")

# Calculate area under mean quality
library(pracma)

area_Fs_filt <- data.frame()

for (i in 1:length(mean_quality_Fs_filt)) {
  area_Fs_filt[i, 1] <- mean_quality_Fs_filt[[i]]$name
  
  quality_sample <- mean_quality_Fs_filt[[i]]$data
  quality_sample$Cycle <- as.numeric(quality_sample$Cycle)
  area_Fs_filt[i, 2] <- pracma::trapz(quality_sample$Cycle, quality_sample$Mean_quality_score)
  
}

area_Rs_filt <- data.frame()

for (i in 1:length(mean_quality_Rs_filt)) {
  area_Rs_filt[i, 1] <- mean_quality_Rs_filt[[i]]$name
  
  quality_sample <- mean_quality_Rs_filt[[i]]$data
  quality_sample$Cycle <- as.numeric(quality_sample$Cycle)
  area_Rs_filt[i, 2] <- pracma::trapz(quality_sample$Cycle, quality_sample$Mean_quality_score)
  
}

area_filt <- data.frame(
  Sample = gsub("_F", "", area_Fs_filt$V1),
  Area_F = area_Fs_filt$V2,
  Area_R = area_Rs_filt$V2
)

saveRDS(object = area_filt, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/RDS_objects/area_quality_filt.rds")


# Step 4: Dereplicate ---- 
# Dereplicate
## Here we take the unique reads and generate a count file (sort of)    
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names[exists]
names(derepRs) <- sample.names[exists]
# Verify dereplication results   
derepFs  
derepRs


# Step 5: Estimate errors with DADA algorithm ----
# Error estimation
# The error rates for each possible transition (eg. A->C, A->G, …) are shown.

# Learn the error rates. This takes a while. 
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
# Look at the error rates used as input in every step of the selfconsist algorithm
dadaFs.lrn[[1]]$err_in
dadaRs.lrn[[1]]$err_in     
# Obtain error tables     
errF <- dadaFs.lrn[[1]]$err_out 
head(errF)
errR <- dadaRs.lrn[[1]]$err_out 
head(errR) 

save(dadaFs.lrn,file="DADA2/dadaFs.lrn.Rda")
save(dadaRs.lrn,file="DADA2/dadaRs.lrn.Rda")

# Step 5b
load(file="DADA2/dadaFs.lrn.Rda")
load(file="DADA2/dadaRs.lrn.Rda")
# Plot the error rates  
plotErrors(dadaFs.lrn, nominalQ=TRUE,err_in=TRUE)
plotErrors(dadaRs.lrn, nominalQ=TRUE,err_in=TRUE)

# Check DADA options
getDadaOpt()


# Step 6: Building an error model ----
# Sample inference
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaFs
dadaFs[[1]]$err_out
dadaFs[[1]]$err_in
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
dadaRs
dadaRs[[1]]$err_out
dadaRs[[1]]$err_in

# Step 6b: Pooled building of error model ----
dadaFS_pooled <- dada(derepFs, err = errF, multithread = TRUE, pool = TRUE)
dadaFS_pooled
dadaFS_pooled[[1]]$err_out
dadaFS_pooled[[1]]$err_in
dadaRS_pooled <- dada(derepRs, err = errF, multithread = TRUE, pool = TRUE)
dadaRS_pooled
dadaRS_pooled[[1]]$err_out
dadaRS_pooled[[1]]$err_in


# Step 7: Merging paired end reads ----
# Merge paired ends. Merge forward and reverse complement of the reverse read (denoised sequences). Merged sequences are only output if the forward and reverse reads overlap by at least 12 bases, and are identical to each other in the overlap region   
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
head(mergers[[1]])  


# Step 8: Creating ASV table ----
# Construct sequence table with columns = ASVs and rows = samples. The lenght of the sequences varies between 401-428, which is as expected given the used 341F/785R primer pair 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab))) # Distribution of length of sequences
head(seqtab)

write.csv2(seqtab, file = "/Projects1/Fabian/Oral_microbiome/StrainRecognitionFCM/DADA2/DADA2/seqtab.csv")


# Step 9: Removing chimeras ----
# Chimera removal
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

# Calculate percentage of chimeric ASVs and reads
(dim(seqtab)[2]-dim(seqtab.nochim)[2])/dim(seqtab)[2] # 65% of ASVs were chimeric
1-sum(seqtab.nochim)/sum(seqtab) # Only 3.6 % of the reads were chimeric 


# Step 10: Check reads through pipeline ----
# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
out.df <- as.data.frame(out)
rownames(out.df) <- sample.names
processed <- as.data.frame(cbind(sapply(derepFs, getN), sapply(derepRs, getN), sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.nochim !=0)))
track <- merge(out.df,processed, by=0, all=TRUE)
rownames(track) <- track$Row.names
track <-  subset(track, select = -c(Row.names) )
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "derepFs", "derepRs","denoisedF", "denoisedR", "merged", "nonchim", "ASVs")
track[is.na(track)] <- 0
track <- rbind(track, colSums(track))
rownames(track) <- c(sample.names, "Sum")
track


xlsxloc <- "DADA2/DADA2_processor.xlsx"
wb <- openxlsx::write.xlsx(x = track, file = "DADA2/DADA2_processor.xlsx", sheetName = "track", rowNames = TRUE)


# Step 11: Format ASV table ----
# Format Output: amplicon sequence variant (ASV) table, a higher-resolution analogue of the traditional OTU table, which records the number of times each exact amplicon sequence variant was observed in each sample. 
ASVtable <- seqtab.nochim
ASVfasta <- colnames(ASVtable)
ASV <- paste("ASV",seq(1,length(ASVfasta),1),sep="")
ASVfasta <- cbind(ASV, ASVfasta)
openxlsx::addWorksheet(wb,sheetName = "ASVfasta")
openxlsx::writeData(wb,sheet="ASVfasta",x=ASVfasta)
openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
colnames(ASVtable) <- ASV



# Step 12: Assign taxonomy ----
# Assign taxonomy: More info: https://benjjneb.github.io/dada2/training.html. Available DADA2 formatted reference taxonomies: Silva, RDP and greengenes are available. The Naive Bayesian Classifier is used for taxonomy assignments. The reference taxonomy files are copied to the DADA2 folder. A 50% cut-off is used as the bootstrap confidence for displaying taxonomy. Additionally, species level assignments can be done based on exact matching between ASVs and sequenced reference strains. Recent analysis suggests that exact matching (or 100% identity) is the only appropriate way to assign species to 16S gene fragments. Currently, species-assignment training fastas are available for the Silva and RDP 16S databases. 
if(database=="silva") {
  taxa_silvav138 <- assignTaxonomy(seqtab.nochim,"/Taxonomies/DADA2tax/silva_nr99_v138.1_train_set.fa", multithread=TRUE, minBoot = 50,tryRC=TRUE,outputBootstraps=TRUE)
  taxa_silvav138.df <- as.data.frame(taxa_silvav138)
  taxa_silvav138.df$seq <- rownames(taxa_silvav138.df)
  rownames(taxa_silvav138.df) <- ASV
  names(taxa_silvav138.df) = gsub(pattern = "tax.", replacement = "", x = names(taxa_silvav138.df))
  taxa_silvav138.genus <- select(taxa_silvav138.df, -contains("boot"))
  taxa_silvav138.genus <- select(taxa_silvav138.genus, -contains("seq"))
  
  if(species=="yes") {
    taxa_silvav138.species <- data.frame(addSpecies(taxa_silvav138[[1]], "/Taxonomies/DADA2tax/silva_species_assignment_v138.1.fa.gz",tryRC=TRUE,allowMultiple=3))
    rownames(taxa_silvav138.species) <- ASV
  }
}  
if(database=="unite") {
  taxa_UNITE <- assignTaxonomy(seqtab.nochim,"/Taxonomies/DADA2tax/UNITE/sh_general_release_s_all_10.05.2021/sh_general_release_dynamic_s_all_10.05.2021.fasta", multithread=TRUE, minBoot = 50,tryRC=TRUE,outputBootstraps=TRUE)
  taxa_UNITE.df <- as.data.frame(taxa_UNITE)
  taxa_UNITE.df$seq <- rownames(taxa_UNITE.df)
  rownames(taxa_UNITE.df) <- ASV
  names(taxa_UNITE.df) = gsub(pattern = "tax.", replacement = "", x = names(taxa_UNITE.df))
  taxa_UNITE.genus <- select(taxa_UNITE.df, -contains("boot"))
  taxa_UNITE.genus <- select(taxa_UNITE.genus, -contains("seq"))
}
if(database=="midas") {
  taxa_midas <- assignTaxonomy(seqtab.nochim,"/Taxonomies/DADA2tax/DADA2_file_MiDAS_4.8.1.fa", multithread=TRUE, minBoot = 50,tryRC=TRUE,outputBootstraps=TRUE)
  taxa_midas.df <- as.data.frame(taxa_midas)
  taxa_midas.df$seq <- rownames(taxa_midas.df)
  rownames(taxa_midas.df) <- ASV
  names(taxa_midas.df) = gsub(pattern = "tax.", replacement = "", x = names(taxa_midas.df))
  taxa_midas.genus <- select(taxa_midas.df, -contains("boot"))
  taxa_midas.genus <- select(taxa_midas.genus, -contains("seq"))
}

# Step 13: Combine taxonomy and ASV table ----
# Combine taxonomy with abundance ASV table
if(database=="silva") {
  if(species=="yes") {
    ASVtax138 <- dplyr::select(taxa_silvav138.species,-Species)
    ASVtax138.species <- cbind(taxa_silvav138.species[colnames(ASVtable),],t(ASVtable))
  } else {
    ASVtax138 <- taxa_silvav138.genus
  }
  ASVtax138 <- cbind(ASVtax138[colnames(ASVtable),],t(ASVtable))
  openxlsx::addWorksheet(wb,sheetName = "ASVsilva138")
  openxlsx::writeData(wb,sheet="ASVsilva138",x=ASVtax138, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
  if(species=="yes") {
    openxlsx::addWorksheet(wb,sheetName = "ASVsilva138.species")
    openxlsx::writeData(wb,sheet="ASVsilva138.species",x=ASVtax138.species, rowNames = TRUE)
    openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
  }
  openxlsx::addWorksheet(wb,sheetName = "Bootstrap_values")
  openxlsx::writeData(wb,sheet="Bootstrap_values",x=taxa_silvav138.df, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
}
if(database=="unite") {
  ASVtaxUNITE <- taxa_UNITE.genus
  ASVtaxUNITE <- cbind(ASVtaxUNITE[colnames(ASVtable),],t(ASVtable))
  openxlsx::addWorksheet(wb,sheetName = "ASVsilvaUNITE")
  openxlsx::writeData(wb,sheet="ASVsilvaUNITE",x=ASVtaxUNITE, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
  openxlsx::addWorksheet(wb,sheetName = "Bootstrap_values_UNITE")
  openxlsx::writeData(wb,sheet="Bootstrap_values_UNITE",x=taxa_UNITE.df, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
}
if(database=="midas") {
  ASVtaxmidas <- taxa_midas.genus
  ASVtaxmidas <- cbind(ASVtaxmidas[colnames(ASVtable),],t(ASVtable))
  openxlsx::addWorksheet(wb,sheetName = "ASVtaxmidas")
  openxlsx::writeData(wb,sheet="ASVtaxmidas",x=ASVtaxmidas, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
  openxlsx::addWorksheet(wb,sheetName = "Bootstrap_values_midas")
  openxlsx::writeData(wb,sheet="Bootstrap_values_midas",x=taxa_midas.df, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
}


# Step 14: Combine taxonomy and proportional ASV table ----
# Combine taxonomy with abundance ASV table
if(database=="silva") {
  if(species=="yes") {
    ASVtax <- dplyr::select(taxa_silvav138.species,-Species)
  } else {
    ASVtax <- taxa_silvav138.genus
  }
  ASVtax_prop <- sweep(ASVtable,1,rowSums(ASVtable),'/')*100
  ASVtax_prop_138 <- cbind(ASVtax[colnames(ASVtax_prop),],t(ASVtax_prop))
  openxlsx::addWorksheet(wb,sheetName = "ASVsilva138_prop")
  openxlsx::writeData(wb,sheet="ASVsilva138_prop",x=ASVtax_prop_138, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
}
if(database=="unite") {
  ASVtax <- taxa_UNITE.genus
  ASVtax_prop <- sweep(ASVtable,1,rowSums(ASVtable),'/')*100
  ASVtax_prop_UNITE <- cbind(ASVtax[colnames(ASVtax_prop),],t(ASVtax_prop))
  openxlsx::addWorksheet(wb,sheetName = "ASVsilvaUNITE_prop")
  openxlsx::writeData(wb,sheet="ASVsilvaUNITE_prop",x=ASVtax_prop_UNITE, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
}
if(database=="midas") {
  ASVtax <- taxa_midas.genus
  ASVtax_prop <- sweep(ASVtable,1,rowSums(ASVtable),'/')*100
  ASVtax_prop_midas <- cbind(ASVtax[colnames(ASVtax_prop),],t(ASVtax_prop))
  openxlsx::addWorksheet(wb,sheetName = "ASVmidas_prop")
  openxlsx::writeData(wb,sheet="ASVmidas_prop",x=ASVtax_prop_midas, rowNames = TRUE)
  openxlsx::saveWorkbook(wb,file="DADA2/DADA2_processor.xlsx",overwrite=TRUE)
}


# Overview

# Step 15: Handoff to Phyloseq ----
# Combine taxonomy with abundance ASV table
if(database=="silva") {
  if(species=="yes") {
    taxexport <- dplyr::select(taxa_silvav138.species,-Species)
    write.csv(taxa_silvav138.species,"DADA2/ASVtax138.species.csv")
  } else {
    taxexport <- taxa_silvav138.genus
  }
  write.csv(taxexport,"DADA2/ASVtax138.csv")
  ASVtable.t <- t(ASVtable)
  write.csv(ASVtable.t,"DADA2/ASVtable.csv")
  write.csv(ASVfasta,"DADA2/ASVfasta.csv", row.names=FALSE)
  write.csv(track,"DADA2/track.csv", row.names=TRUE)
  save(track,file="DADA2/track.Rda")
  kable(track)
}
if(database=="unite") {
  taxexport <- taxa_UNITE.genus
  write.csv(taxexport,"DADA2/ASVtaxUNITE.csv")
  ASVtable.t <- t(ASVtable)
  write.csv(ASVtable.t,"DADA2/ASVtable.csv")
  write.csv(ASVfasta,"DADA2/ASVfasta.csv", row.names=FALSE)
  write.csv(track,"DADA2/track.csv", row.names=TRUE)
  save(track,file="DADA2/track.Rda")
  kable(track)
}

if(database=="midas") {
  taxexport <- dplyr::select(taxa_midas.genus,-Species)
  write.csv(taxa_midas.genus,"DADA2/ASVtaxmidas.species.csv")
  write.csv(taxexport,"DADA2/ASVtaxmidas.csv")
  ASVtable.t <- t(ASVtable)
  write.csv(ASVtable.t,"DADA2/ASVtable.csv")
  write.csv(ASVfasta,"DADA2/ASVfasta.csv", row.names=FALSE)
  write.csv(track,"DADA2/track.csv", row.names=TRUE)
  save(track,file="DADA2/track.Rda")
  kable(track)
}
