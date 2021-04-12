##############################################################################
## Script to: 
##  1) read in QMP MTs Data
##  2) perform linear models on continuous traits
##  3) extract residuals
##  4) write a sample file
##
## By: David A Hughes; Date: Dec 5th 2020
##
## Built off of work by: Laura Corbin, Caroline Bull, David Hughes
##                        Date begun: Feb 29 2020
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
data_file_in = paste0(data_processed_dir,"qmp/fgfp_qmp_gwas_dataset_step2_transformed.txt")
data_file_out = paste0(data_processed_dir,"qmp/fgfp_qmp_gwas_dataset_step3_residuals.txt")

## SAMPLE files 
original_sam_file = paste0(data_input_dir,"genotypes/bgen/data.imputed.sample")
n = paste0("mkdir -p ", data_processed_dir,"samplefile")
system(n)
sam_file_out = paste0(data_processed_dir,"samplefile/fgfp_qmp.sample")
samples_in_gwas = paste0(data_processed_dir,"samplefile/samples_in_gwas.txt")

#######################################
##
## (III) source our FUNCTIONS file
##  
#######################################
source( paste0(working_dir, "/scripts/fgfp_qmp_gwas/FUNCTIONS.R"))

#######################################
##
## (IV) Start a log file
##
#######################################
cat( paste0("Starting a log file - it will be in your processed data metabolites directory.\n") )


logfilename = paste0(data_processed_dir, "qmp/", project, "_08_data_lm_residuals_", today,  "_logfile.txt")
sink(file = logfilename , split = TRUE  )

## send information to log file
cat( paste0("Fitting a linear model to continuous traits and extracting residuals.\n") )
cat( paste0("\t-- Todays date is: ", today, ".\n\n") )

#######################################
##
## (V) Read in the data
##
#######################################


########################
## (Va) Read in metabolite data 
########################
cat( paste0("I) Reading in QMP & Clinical/Covariate Data.\n") )
#### read in the data
mdata = read.table( data_file_in , header = TRUE, as.is = TRUE)

## Sample and metabolite count feedback
suffix = c("_AB","_truncAB","_PA")
mt_locations = sapply(suffix, function(x){
        w = grep(x, colnames(mdata))
        return(w)
        })
mt_locations = unlist(mt_locations)
#####
f_count = length( mt_locations )
cat(paste0("\ta) There are ", f_count, " microbiome triats (MTs) in the QMP & clinical data file \n"))
cat(paste0("\tb) There are ", nrow(mdata), " samples in the QMP & clinical data file.\n"))


########################
## (Vb) Factor conversion
########################
cat( paste0("   -- convert / insure binary covariates (sex) and PA traits are factors.\n") )
mdata$sex = as.factor(mdata$sex)
####
pa = grep("_PA", colnames(mdata))
for(x in pa){
  mdata[, x] = as.factor(mdata[, x])
}

#######################################
##
## (VI - A) Linear Model Fitting
##  for cell_counts_per_g
#######################################
fit = lm(cell_counts_per_g ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
    data = mdata,
    na.action = na.exclude)
res = residuals(fit)
mdata$cell_counts_per_g = res

#######################################
##
## (VI) Linear Model Fitting
##
#######################################
cat( paste0("II) Running the linear models on continuous traits.\n") )
cat( paste0("\ti)covariates include:\n") )
cat( paste0("\t\t- age\n") )
cat( paste0("\t\t- sex\n") )
cat( paste0("\t\t- PC1-PC10\n") )

## identify traits to run
w = c( grep("_AB", colnames(mdata)), grep("_truncAB", colnames(mdata)))
MTcols = colnames(mdata)[w]

###
lm_models = list()
resdata = c()
## run a loop
for(mt in MTcols){
fit = lm(mdata[,mt] ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, 
    data = mdata,
    na.action = na.exclude)
  ### extract beta, se, and pvals
  coef = summary(fit)$coef
  ## extract adjusted R squared
  adjR2 = summary(fit)$adj.r.squared
  ## place coef and R2 in a list
  lm_models[[mt]] = list(coef = coef, adjR2 = adjR2)
  ### extract residuals
  res = residuals(fit)
  ## place residuals in a matrix
  resdata = cbind(resdata, res)

}
## redefine column | metabolite names
colnames(resdata) = names(lm_models)

#### write linear model coefs and R2 estimates to an Robj file
cat( paste0("\tii) linear models complete\n") )
cat( paste0("\t\t(a) model coefficents and adjusted R squared values, for ", length(lm_models)," MTs, saved as a list (lm_models), are written to an Robj file.\n") )
n = paste0(data_processed_dir,"qmp/lm_models.Robj")
save(lm_models, file = n)

## change owner of the file
system( paste0( "chown :nic_timpson_grp ", n) )

## change  permission of the file
system( paste0( "chmod 770 ", n) )


#######################################
##
## (VIII) WRITE DATA TO FILE
##
#######################################
cat( paste0("IV) Writing new metabolon data set to a .txt file. \n") )

cat( paste0("\ti) placing residual data into our dataset.\n") )

m = match(colnames(resdata), colnames(mdata) )
mdata[,m] = resdata

## write generated data set to flat text file
write.table(mdata, file = data_file_out , row.names = TRUE, col.names = TRUE , quote=F)

## change owner of the file
system( paste0( "chown :nic_timpson_grp ", data_file_out) )

## change  permission of the file
system( paste0( "chmod 770 ", data_file_out) )


#######################################
##
## (IX) WRITE A SAMPLE FIlE to file
##
#######################################
cat( paste0("V) Writing a SAMPLEFILE of the dataset. \n") )

## !!!!!!!  VERY IMPORTANT STEP  !!!!!!!
sam_original = read.table(original_sam_file, h = TRUE)
### EXAMPLE HEAD
# ID_1 ID_2 missing father mother sex plink_pheno
# 0 0 0 D D D C
# VDP.01054_C10 VDP.01054_C10 0 0 0 2 -9
# VDP.03187_A01 VDP.03187_A01 0 0 0 1 -9

## SAMPLE FILE format definition for row two
# 0 (for the first identifier column)
# D (for a column containing discrete values, e.g. a set of strings)
# P or C - for columns containing continuous value - each value must be numerical, or a missing value.
# B - for a column containing a binary trait. The values in this column must be '0', '1', 'control', or 'case'.

## List of samples
ids_in_gwas = data.frame(ids = mdata$genetic_id)
write.table(ids_in_gwas, file = samples_in_gwas , row.names = FALSE, col.names = FALSE , quote=F, sep = " ")

## change owner of the file
system( paste0( "chown :nic_timpson_grp ", samples_in_gwas) )

## change  permission of the file
system( paste0( "chmod 770 ", samples_in_gwas) )

### NEW SAMPLE FILE
new_sam_file = data.frame(
  ID_1 = mdata$genetic_id, 
  ID_2 = mdata$id, 
  missing = rep(0, nrow(mdata)), 
  father = rep(0, nrow(mdata)), 
  mother = rep(0, nrow(mdata)), 
  sex = mdata$sex,
  age = mdata$age,
  bmi = mdata$BMI,
  height = mdata$height,
  weight = mdata$weight,
  bacterial_cell_count_per_g = mdata$cell_counts_per_g
  )

## add PCs
w = grep("PC", colnames(mdata))
new_sam_file = cbind(new_sam_file, mdata[,w])

## add MTs
suffix = c("_AB","_truncAB","_PA")
mt_locations = sapply(suffix, function(x){
        w = grep(x, colnames(mdata))
        return(w)
        })
mt_locations = sort( unlist(mt_locations) )

new_sam_file = cbind(new_sam_file, mdata[,mt_locations])

# add variable specification line
cat( paste0("\ti) adding variable specification line. \n") )
samfile = new_sam_file

vs = rep("C", ncol(samfile))
names(vs) = colnames(samfile)

## redefine IDs
w = grep("ID", names(vs))
vs[w] = "0"

## redefine missing
w = grep("missing", names(vs))
vs[w] = "0"

## redefine father as discrete
w = grep("father", names(vs))
vs[w] = "D"

## redefine mother as discrete
w = grep("mother", names(vs))
vs[w] = "D"

## redefine sex as discrete
w = grep("sex", names(vs))
vs[w] = "D"

## redefine all continuous MT as Phenotypes
w = grep("_truncAB", names(vs))
vs[w] = "P"

## redefine all continuous MT as Phenotypes
w = grep("_AB", names(vs))
vs[w] = "P"


## redefine all P|A metabolites as Binary
w = grep("_PA", names(vs))
vs[w] = "B"

## redefine bacterial_cell_count_per_g as Phenotype
w = grep("bacterial_cell_count_per_g", names(vs))
vs[w] = "P"


temp = apply(samfile, 2, function(x){as.character(x)})
# merge df and variable specification line
samfile <- rbind(vs, temp)
rm(temp)

cat( paste0("\tii) writing samplefile to file. \n") )
### write samfple file to file.
write.table(samfile, file = sam_file_out , row.names = FALSE, col.names = TRUE , quote=F, sep = " ")

## change owner of the file
system( paste0( "chown :nic_timpson_grp ", sam_file_out) )

## change  permission of the file
system( paste0( "chmod 770 ", sam_file_out) )

#######################################
##
## (X) Build the sumission array pairings
##
#######################################
cat( paste0("VI) Build submision array pairing for qsub job array. \n") )

## extract MT names
suffix = c("_AB","_truncAB","_PA")
mt_locations = sapply(suffix, function(x){
        w = grep(x, colnames(samfile))
        return(w)
        })
mt_locations = unlist(mt_locations)
MTs = c("bacterial_cell_count_per_g", colnames(samfile)[mt_locations])

## which ones are P|A
pa = grep("_PA", MTs )

## set up a vector of P|A
PA = MTs[pa]
## set up a vector of Continuous traits
CONT = MTs[-pa]

###################
## Continuous trait
##  pairs
###################
cat( paste0("\ti) building continuous trait array pairing file. \n") )
out = c()
for(i in CONT){
  for(j in 1:22){
      if( j  < 10 ){
      j = paste0("0", j)
    }
    o = c(j, i)
    out = rbind(out, o)
  }
}
#####
cat( paste0("\tii) writing continuous_traits_jobarray.txt to file. \n") )
n = paste0(working_dir,"scripts/fgfp_qmp_gwas/continuous_traits_jobarray.txt")
write.table(out, file = n, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

## change owner of the log file
system( paste0( "chown :nic_timpson_grp ", n) )

## change  permission of the log file
system( paste0( "chmod 770 ", n) )



###################
## P|A trait
##  pairs
###################
cat( paste0("\tiii) building presence/abcense trait array pairing file. \n") )
out = c()
for(i in PA){
  for(j in 1:22){
      if( j  < 10 ){
      j = paste0("0", j)
    }
    o = c(j, i)
    out = rbind(out, o)
  }
}
#####
cat( paste0("\tiv) writing preabsence_traits_jobarray.txt to file. \n") )
n = paste0(working_dir,"scripts/fgfp_qmp_gwas/preabsence_traits_jobarray.txt")
write.table(out, file = n, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

## change owner of the log file
system( paste0( "chown :nic_timpson_grp ", n) )

## change  permission of the log file
system( paste0( "chmod 770 ", n) )


#########################
###### END SESSION ######
#########################
cat(paste0( "\n\n **** FINISHED **** \n\n\nR Session Info:\n") )

# print out details
print(sessionInfo())

## stop writing to log file
sink()


## change owner of the log file
system( paste0( "chown :nic_timpson_grp ", logfilename) )

## change  permission of the log file
system( paste0( "chmod 770 ", logfilename) )
