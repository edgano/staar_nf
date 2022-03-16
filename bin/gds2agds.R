rm(list=ls())
gc()

##########################################################################
#           Input
##########################################################################

### gds file
dir_geno <- "/lustre/scratch123/hgi/teams/hgi/mo11/associations/GDS/"
gds_file_name_1 <- "Interval_WGS_chr20_TF_binding_site_test.gds"

### annotation file
dir_anno <- "/lustre/scratch123/hgi/teams/hgi/mo11/associations/"
anno_file_name_1 <- "Anno_chr"
anno_file_name_2 <- "_STAARpipeline.csv"

chr <- as.numeric(commandArgs(TRUE)[1])
chr <- as.numeric('20')

###########################################################################
#           Main Function 
###########################################################################

### load required package
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(readr)

### read annotation data
FunctionalAnnotation <- read_csv(paste0(dir_anno,"chr",chr,"/",anno_file_name_1,chr,anno_file_name_2),
col_types=list(col_character(),col_double(),col_double(),col_double(),col_double(),
col_double(),col_double(),col_double(),col_double(),col_double(),
col_character(),col_character(),col_character(),col_double(),col_character(),
col_character(),col_character(),col_character(),col_character(),col_double(),
col_double(),col_character()))

dim(FunctionalAnnotation)

## rename colnames
colnames(FunctionalAnnotation)[2] <- "apc_conservation"
colnames(FunctionalAnnotation)[7] <- "apc_local_nucleotide_diversity"

## open GDS
gds.path <- paste0(dir_geno,gds_file_name_1)
genofile <- seqOpen(gds.path, readonly = FALSE)

Anno.folder <- addfolder.gdsn(index.gdsn(genofile, "annotation/info"), "FunctionalAnnotation")
add.gdsn(Anno.folder, "FunctionalAnnotation", val=FunctionalAnnotation, compress="LZMA_ra", closezip=TRUE)

seqClose(genofile)