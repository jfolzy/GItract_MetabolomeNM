#######################################
# source base configuration
#######################################

# load libraries
library(ggplot2)
library(tidyverse)
library(here)
library(corrplot)
library(cowplot)
library(readxl)
library(Hmisc)
library(ggpubr)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(dynamicTreeCut)
library(ggrepel)
library(pvclust)
library(parallel)
library(ggbeeswarm)
library(FactoMineR)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ropls")
library(ropls) 
library(ggdendro)
library(dendextend)
library(scales)
library(lme4)
library(lmerTest)
library(glmnet)

#--------------------------------------------
# define directories
#--------------------------------------------
CleanDataPath <- "1_data/clean_data/"

#--------------------------------------------
# Set colors for capsule types & stool
#--------------------------------------------

col<-c('#2364AA','#3DA5D9','#134611','#90D733','#743C2F')
CapTypeAndStoolColors <- col
CapType <- c(col[1:4])
CapAndStoolColors <- c(col[1],col[5])
CapQCsStoolColors <- c(col[1:4],"red",col[5])

#--------------------------------------------
#Load Data Locations
#--------------------------------------------
SampleMetaData_path <- paste0(CleanDataPath, "SampleMetaData_March2022.RDS")
FoodMetaData_path <- paste0(CleanDataPath, "FoodMetaData_March2022.RDS")
NonNorm_MetDF_path <- paste0(CleanDataPath, "NoNormalization_MetaboliteDataFrame.RDS")
iSTDNorm_MetDF_path <- paste0(CleanDataPath, "iSTDNormalization_MetaboliteDataFrame.RDS")
mTICNorm_MetDF_path <- paste0(CleanDataPath, "mTICNormalization_MetaboliteDataFrame.RDS")
SERRFNorm_MetDF_path <- paste0(CleanDataPath, "Serrf_MetaboliteDataFrame.RDS")
BileAcidsDF_path <- paste0(CleanDataPath, "BileAcidsOnly_DataFrame.RDS")

#ColorSchemeTables
SampleTypeColorScheme_path <- paste0(CleanDataPath, "SampleTypeColorScheme.RDS")
SubjectColorScheme_path <- paste0(CleanDataPath, "SubjectColorScheme.RDS")

#Normalized Data and Subsets
MS_MetaData_path <- paste0(CleanDataPath, "MS_MetaDataFrame.RDS")
MetOntology_path <- paste0(CleanDataPath, "MetaboliteOntology_March2022.RDS")
NormMetDF_path <- paste0(CleanDataPath, "NormalizedMetaboliteDataFrame.RDS")
NormMetDF_Unknowns_path <- paste0(CleanDataPath, "NormalizedMetaboliteDataFrame_Unknowns.RDS")
tNormMetDF_CSQ_path <- paste0(CleanDataPath, "tNormalizedLogMetaboliteDataFrame_CapStoolQCs.RDS")
tNormMetDF_CS_path <- paste0(CleanDataPath, "tNormalizedLogMetaboliteDataFrame_CapsStool.RDS")
tNormMetDF_C_path <- paste0(CleanDataPath, "tNormalizedLogMetaboliteDataFrame_Capsules.RDS")
tNormMetDF_CS_Uk_path <- paste0(CleanDataPath, "tNormalizedLogMetaboliteDataFrame_Uknowns_CS.RDS")
tNormMetDF_C_Uk_path <- paste0(CleanDataPath, "tNormalizedLogMetaboliteDataFrame_Unknowns_C.RDS")
logNormMetDF_CSQ_path <- paste0(CleanDataPath, "LogNormalizedMetaboliteDataFrame_CapStoolQCs.RDS")
logNormMetDF_CS_path <- paste0(CleanDataPath, "LogNormalizedMetaboliteDataFrame_CapsStool.RDS")
logNormMetDF_C_path <- paste0(CleanDataPath, "LogNormalizedMetaboliteDataFrame_Capsules.RDS")
logNormMetDF_Uk_CS_path <- paste0(CleanDataPath, "LogNormalizedMetaboliteDataFrame_Uk_CS.RDS")
logNormMetDF_Uk_C_path <- paste0(CleanDataPath, "LogNormalizedMetaboliteDataFrame_Uk_C.RDS")
FiltLogNormMetDF_CSQ_path <- paste0(CleanDataPath, "FilteredLogNormalizedMetaboliteDataFrame_CapsulesStoolQCs.RDS")
FiltLogNormMetDF_CS_path <- paste0(CleanDataPath, "FilteredLogNormalizedMetaboliteDataFrame_CapsulesStool.RDS")
FiltLogNormMetDF_C_path <- paste0(CleanDataPath, "FilteredLogNormalizedMetaboliteDataFrame_Capsules.RDS")
FiltLogNormMetDF_Uk_CS_path <- paste0(CleanDataPath, "FilteredLogNormalizedMetaboliteDataFrame_Uk_CS.RDS")
FiltLogNormMetDF_Uk_C_path <- paste0(CleanDataPath, "FilteredLogNormalizedMetaboliteDataFrame_Uk_C.RDS")

#Caffeine Corr Tables
CorrToCaffeine_path <- paste0(CleanDataPath,"TopCaffeineCors.RDS")