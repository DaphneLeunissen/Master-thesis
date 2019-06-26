### Master thesis 
### R script for result comparison
### Galaxy workflow vs JSI data analysis
### Author: Daphne Leunissen
### Date: 12-04-2019

#-----------------------------------------------------------------------------#
# Install and load required packages
#-----------------------------------------------------------------------------#
# first install required packages if needed

library(tidyverse)  
library(xlsx)       

options(stringsAsFactors = FALSE)

#-----------------------------------------------------------------------------#
# Set working directory and import data files 
#-----------------------------------------------------------------------------#

setwd("C:\\Users\\daphn\\OneDrive\\Documenten\\Master Systems Biology\\Master thesis\\R\\Results")

# select the data subject for comparison 
file = "1284"

# import JSI results
data_JSI <- read.delim(paste0(file, "_JSI.txt"), stringsAsFactor = F, header=TRUE , sep="\t")

# import annovar results 
data_annovar <- read.delim(paste0(file,"_annovar.txt"), stringsAsFactor = F, header=TRUE , sep="\t")

#-----------------------------------------------------------------------------#
# Data preparation 
#-----------------------------------------------------------------------------#
# Annovar data preparation
#-----------------------------------------------------------------------------#
# remove columns in annovar data containing only NA values
data_annovar1 <- data_annovar[,colSums(is.na(data_annovar))<nrow(data_annovar)]

# remove columns in annovar data which are not needed/won't be used:
# reference, observed, CHROM, ID, FILTER, FORMAT, unknown, UCSCKnownGene_Func, UCSCKnownGene_Gene, UCSCKnownGene_ExonicFunc, 
# UCSCKnownGene_AAChange, EnsemblGene_Func, EnsemblGene_ExonicFunc, EnsemblGene_AAChange
data_annovar2 <- data_annovar1[, -c(4, 5, 6, 8, 12, 14, 15, 20, 21, 22, 23, 24, 26, 27)]

# only remain the variants with an exonic effect
data_annovar3 <- subset(data_annovar2, RefSeq_Func == "exonic")

## separate the information given in the info field of annovar data
variant_info <- as.data.frame(strsplit(data_annovar3$INFO, ";"))
variant_info <- as.data.frame(t(variant_info))
row.names(variant_info) <- row.names(data_annovar3)

# set column names 
colnames(variant_info) <- c("AB", "ABP", "AC", "AF", "AN", "AO", "CIGAR", "DP", "DPB", "DPRA", "EPP", "EPPR", "GTI", "LEN", "MEANALT", "MQM",
                            "MQMR", "NS", "NUMALT", "ODDS", "PAIRED", "PAIREDR", "PAO", "PQA", "PQR", "PRO", "QA", "QR", "RO", "RPL", "RPP",
                            "RPPR", "RPR", "RUN", "SAF", "SAP", "SAR", "SRF", "SRP", "SRR", "TYPE")
# create an empty matrix with same size as variant_info
variant_info2 <- data.frame(matrix(nrow = nrow(variant_info), ncol = ncol(variant_info)))

# only remain with the values in the columns, remove the parameter name and the equal sign
for (cname in colnames(variant_info)) {
  variant_info2[, cname] <- gsub("^[^=]+=", "", variant_info[, cname])
}

# remove columns with NA values 
variant_info2 <- variant_info2[,colSums(is.na(variant_info2))<nrow(variant_info2)]

#-----------------------------------------------------------------------------#
# JSI data preparation 
#-----------------------------------------------------------------------------#
# remove columns in JSI data containing only NA values
data_JSI <- data_JSI[,colSums(is.na(data_JSI))<nrow(data_JSI)]

# remove all FP, 3' UTR and 5'UTR.
data_JSI1 <- subset(data_JSI, mut.Effect != "5' UTR")
data_JSI1 <- subset(data_JSI1, mut.Effect != "3' UTR")
data_JSI1 <- subset(data_JSI1, mut.Effect != "FP")

# remove rows containing these gene: BIRC3, SLC7A8, ZNF2.
# these are msi markers and are not taken into account during data analysis
data_JSI1 <- subset(data_JSI1, Gene != "BIRC3")
data_JSI1 <- subset(data_JSI1, Gene != "SLC7A8")
data_JSI1 <- subset(data_JSI1, Gene != "ZNF2")
data_JSI1 <- subset(data_JSI1, Gene != "MSH2")

# remove Index, Type and, Nuc Change, Weighting column 
rm_col <- c("Index", "Type", "Nuc.Change", "Weighting")
data_JSI2 <- data_JSI1[, ! names(data_JSI1) %in% rm_col]

# split Pos. (position) column to obtain chromosome number and genomic location 
JSI_loc <- do.call(rbind, strsplit(data_JSI2$Pos., " |:"))
JSI_loc <- as.data.frame(JSI_loc)

# remove [ in V5 column
JSI_chr <- gsub("\\[", "", JSI_loc$V5)
JSI_chr <- as.matrix(JSI_chr)

# add chromosome to data_JSI2
data_JSI2 <- add_column(data_JSI2, JSI_chr, .before="Gene")
colnames(data_JSI2)[1] <- "chromosome"

# obtain genomic position, first remove g. 
JSI_pos <- gsub("\\g.", "", JSI_loc$V6)
JSI_pos <- as.matrix(JSI_pos)

# add to data_JSI2 as second column and change column name
data_JSI2 <- add_column(data_JSI2, JSI_pos, .before="Gene")
colnames(data_JSI2)[2] <- "POS"

# remove Pos. column
data_JSI2 <- data_JSI2[, ! names(data_JSI2) %in% "Pos."]

#-----------------------------------------------------------------------------#
# Calculate mutation load/frequency
#-----------------------------------------------------------------------------#
# obtain the AO and RO from variant_info data frame
info_AO <- as.data.frame(str_remove(variant_info$AO, "AO="))
colnames(info_AO) <- "AO"
info_RO <- as.data.frame(str_remove(variant_info$RO, "RO="))
colnames(info_RO) <- "RO"

# combine the two columns together 
info_mut_load <- cbind(info_AO, info_RO)

# calculate the mutation load: AO / AO + RO
mut_load <- (as.numeric(info_mut_load$AO) / (as.numeric(info_mut_load$RO) + as.numeric(info_mut_load$AO))) * 100

# insert the mutation load into annovar data file
data_annovar3 <- add_column(data_annovar3, mut_load, .before ="RefSeq_AAChange")
colnames(data_annovar3)[12] <- "Mutation_load"

#-----------------------------------------------------------------------------#
# Data comparison
#-----------------------------------------------------------------------------#

# merge data_annovar and data_JSI based on chromosome and POS (genomic position) columns. 
final_results <- merge(data_annovar3, data_JSI2, by=c("chromosome", "POS"))

#-----------------------------------------------------------------------------#
# Save results
#-----------------------------------------------------------------------------#
# first select the columns to keep for the report 
keeps <- c("chromosome", "POS", "RefSeq_Gene", "RefSeq_ExonicFunc", "Mutation_load", "Coverage", "c..HGVS")
report_results <- final_results[keeps]
colnames(report_results)[6] <- "Coverage JSI"
# save in excel file in the first sheet
write.xlsx(report_results, paste0(file,"_results.xlsx"), sheetName="report")     

# mutations missed in annovar that were called with JSI:
JSI_var <- data_JSI2$POS[!(data_JSI2$POS %in% final_results$POS)]
JSI_calls <- data_JSI2[is.element(data_JSI2$POS, JSI_var),]
# save in the same excel file in second excel sheet
write.xlsx(JSI_calls, paste0(file,"_results.xlsx"), sheetName="JSI", append=TRUE)

# mutations missed in JSI that were called with annovar:
annovar_var <- data_annovar3$POS[!(data_annovar3$POS %in% final_results$POS)]
annovar_calls <- data_annovar3[is.element(data_annovar3$POS, annovar_var),]
# save in the same excel file in third excel sheet 
write.xlsx(annovar_calls, paste0(file,"_results.xlsx"), sheetName="ANNOVAR", append=TRUE)

