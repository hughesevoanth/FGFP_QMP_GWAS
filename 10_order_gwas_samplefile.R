##############################################################################
## Script to: 
##  1) reorder the GWAS sample.file we made with all of the data
##   	to match the order in the newly generated bgen files
## By: David A Hughes; Date: Dec 5th 2020
##
##############################################################################

#######################################
##
## (I) record todays date
##
#######################################
today = Sys.Date()
today = gsub("-","_",today)

#######################################
##
## (II) Read in Parameter file
##
#######################################

## read in parameter file (specified on command line)
p = commandArgs(trailingOnly=T)[1]
pfile = read.table(p, sep = "=", as.is = TRUE)

# read in and set parameters
project <- as.character( pfile[1,2] )
print(paste0("Project: ",project))
working_dir <- as.character( pfile[2,2] )
print(paste0("Working directory: ",working_dir))
data_input_dir <- as.character( pfile[3,2] )
print(paste0("Input data directory: ",data_input_dir))
data_processed_dir <- as.character( pfile[4,2] )
print(paste0("Processed data directory: ",data_processed_dir))

### Path and names of the PRIMARY data IN and the primary data OUT
data_file_in = paste0(data_processed_dir,"samplefile/fgfp_qmp.sample")
model_samplefile = paste0(data_processed_dir,"genotypes/bgen_gwas/data_chr01.sample")

## OUT FILE
data_file_out = paste0(data_processed_dir,"samplefile/fgfp_qmp_sample_ordered.sample")


######################
## Read in the data
######################
gwasdata = read.table(data_file_in, header = TRUE, as.is = TRUE)
model = read.table(model_samplefile, header = TRUE, as.is = TRUE)

## match sample order
m = match(model$ID_1, gwasdata$ID_1)
gwasdata = gwasdata[m, ]

write.table(gwasdata, file = data_file_out, col.names = TRUE, row.names = FALSE, sep = " ", quote = FALSE)

