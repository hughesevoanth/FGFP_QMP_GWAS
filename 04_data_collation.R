##############################################################################
## Rscript to: 
##  1) read in QMP microbiome Data, QC that data
##  2) then merge that data with clinical covariates
##  3) then produce a *.sample file 
##
## By: David Hughes
## Date begun: Dec 4th 2020
##
##############################################################################
library(readxl)
########################
## record todays date
########################
today = Sys.Date()
today = gsub("-","_",today)

#########################
##
## 1. Read in Parameter file
##
#########################

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

### Path and names of the PRIMARY data OUT
data_file_out = paste0(data_processed_dir,"qmp/fgfp_qmp_gwas_dataset_step1_raw.txt")

########################
## 2. source our 
## 	  FUNCITONS	file
##
########################
source( paste0(working_dir, "/scripts/fgfp_qmp_gwas/FUNCTIONS.R"))


########################
##
## 3. Start a log file
##
########################
logfilename = paste0(data_processed_dir, "qmp/", project, "_04_data_collation_", today,  "_logfile.txt")
sink(file = logfilename , split = TRUE  )

## send information to log file
cat( paste0("Processing QMP and clinical phenotype data to generate a sample file.\n") )
cat( paste0("  -- Todays date is: ", today, ".\n\n") )


############################
##
## 4. READ IN:
##
############################

############################
## 4a. Sample linker file
############################
cat( paste0("I) Reading in linker file.\n") )

###  Read in the data
n = paste0(data_input_dir, "linker/fgfp_sampleid_linker_20200211.txt")
linker = read.table(n, header = TRUE, sep = "\t", as.is = TRUE)
linker$id = gsub("VDP.0","VDP.", linker$fgfp_id)

############################
## 4b. QMP data
############################
cat( paste0("II) Reading in qmp data.\n") )

# read in QMP data
#qmp_data_raw <- read.table(file=paste0(data_input_dir,"qmp/QMP.MR_new.alltax.abundances.tsv"), h=T, sep="\t", stringsAsFactors=F, as.is = TRUE)
qmp_data_raw <- read.table(file=paste0(data_input_dir,"qmp/release3/QMP_GWAS.allranks.final.tsv"), h=T, sep="\t", stringsAsFactors=F, as.is = TRUE)
############################
## 4c. QMP batch data
############################
cat( paste0("II) Reading in qmp batch data.\n") )

# read in QMP data
qmp_batch_data <- read.table(file=paste0(data_input_dir,"qmp/VDP3000_cells.csv"), h=T, sep=",", stringsAsFactors=F, as.is = TRUE)


############################
## 4c. Clinic data
############################
cat( paste0("III) Read in the clinical data.\n") )

clinic = read_excel( paste0(data_input_dir,"clinical/bristol_smoking_data_vdp3000.xlsx"), sheet = 1)
k = c("vdp_ids","gender","age", "BMI","length","weight")
clinic <- clinic[, k]
colnames(clinic) = c("fgfp_id", "sex","age", "BMI","height","weight")
o = order(clinic$fgfp_id)
clinic = clinic[o,]

# read in clinic data
#clinic <- read.table(file= paste0(data_input_dir,"clinical/fgfp_2482_full_metadata_updated.tsv") , h=T, sep="\t", stringsAsFactors=F, as.is = TRUE)
#clinic$fgfp_id = gsub("\\.","\\.0", rownames(clinic))

# extract desired phenotypes
#cat( paste0("     a) Parse the clinical data down to the 8 variables to keep.\n") )
#clinic <- clinic[,c("fgfp_id","Gender","age", "BMI","Height","Weight","Gluc_nuchter_mg.dL","Insuline_nuchter_mUL","Ureum_mgdL", "stoole_score")]
#colnames(clinic) = c("fgfp_id", "sex","age", "BMI","height","weight","fasting_glucose_mg.per.dl","fasting_insulin_mU_per_L", "Ureum_mg_per_dl", "stool_score")

############################
## 4d. Clinic data
############################
cat( paste0("IV) Read in the snpsex data from plink file.\n") )
snpsex = read.table( paste0(data_input_dir,"genotypes/plink/plink.sexcheck"), header = TRUE)

############################
##
## 5.  Merge DATA 
##
############################
cat( paste0("VI) Merge clincical covariate data, qmp and qmp batch data.\n") )

## Who has genotype data?
w = which(linker$genotypedata_good_for_gwas == 1)
ind2keep = linker[w, ]
## remove the duplicates
cat( paste0("VI -) Screen for duplicate samples.\n") )

duplicated_samples = names( which(table( ind2keep$id )>1) )
cat( paste0("     a) There are ",length( duplicated_samples )," duplicated samples.\n" ))
# drop a duplicate
for( i in 1:length(duplicated_samples) ){
	w = which(ind2keep$id == duplicated_samples[i] )
	###
	if(length(w) > 1){
		r = w[-1]
		ind2keep = ind2keep[ -r ,]
	}
}
	
####
cat( paste0("     a) There are ",nrow( ind2keep )," samples with genotype data available.\n" ))
o = order(ind2keep$id)
ind2keep = ind2keep[o, c("id","fgfp_id","genetic_id")]

## Who also has QMP data
w = which( as.character(ind2keep$id) %in% as.character(rownames(qmp_data_raw)) )
cat( paste0("     a) There are ",length( w )," samples who also have QMP data.\n" ))
# x = names( which( table( c(as.character(ind2keep$id), as.character(rownames(qmp_data_raw)) ) ) == 2 ) )
qmp_gwas_samples = ind2keep[w, ]
qmp_gwas_samples = qmp_gwas_samples[, c("id","fgfp_id","genetic_id")]

## Who also has clinical data
w = which( as.character(qmp_gwas_samples$id) %in% as.character(clinic$fgfp_id) )
cat( paste0("     a) There are ",length( w )," samples who also have clinic covariate data.\n" ))
# x = names( which( table( c(as.character(ind2keep$id), as.character(rownames(qmp_data_raw)) ) ) == 2 ) )
qmp_gwas_samples = qmp_gwas_samples[w, ]

#########################
## add in clinical data
#########################
m = match(qmp_gwas_samples$id, clinic$fgfp_id )
qmp_gwas_samples = cbind(qmp_gwas_samples, clinic[m, -1 ])

## add in snpsex data
m = match(qmp_gwas_samples$genetic_id, snpsex$FID )
qmp_gwas_samples = cbind(qmp_gwas_samples, snpsex$SNPSEX[m])
colnames(qmp_gwas_samples)[ncol(qmp_gwas_samples)] = "snpsex"

## add in qmp bath data
m = match(qmp_gwas_samples$id, qmp_batch_data$VDP_ID )
qmp_gwas_samples = cbind(qmp_gwas_samples, qmp_batch_data[m, -1])

## add in qmp data
m = match(qmp_gwas_samples$id, rownames(qmp_data_raw) )
qmp_gwas_samples = cbind(qmp_gwas_samples, qmp_data_raw[m, ])

############################
##
## 6. discard individual
## with mismatched sex
##
############################
cat( paste0("VII) QC: discard individuals with mismatched sex.\n") )
qmp_gwas_samples$sex = as.factor(qmp_gwas_samples$sex)
####
qmp_gwas_samples$snpsex[qmp_gwas_samples$snpsex == 2] = "f"
qmp_gwas_samples$snpsex[qmp_gwas_samples$snpsex == 1] = "m"
qmp_gwas_samples$snpsex = as.factor(qmp_gwas_samples$snpsex)
##
w = which(qmp_gwas_samples$sex != qmp_gwas_samples$snpsex)
if(length(w)>0){
	cat(paste0("     a) A total of  ", length(w), " individuals are being removed for mismatched sex.\n"))
	qmp_gwas_samples = qmp_gwas_samples[-w, ]	
}


############################
##
## 7. check and report dimensions
##
############################
cat( paste0("VIII) Final Report on data set sample size.\n") )

cat(paste0("     a) There are ", ncol(qmp_data_raw), " features in the qmp dataset.\n"))
cat(paste0("     b) There are ", nrow(qmp_gwas_samples), " samples in the QMP dataset after merging datasets removing duplicates.\n"))


############################
##
## 10. WRITE Raw Data set to file
##
############################
cat( paste0("IX) Writing metabolon data set a .txt file. \n") )

## write generated data set to flat text file
write.table(qmp_gwas_samples, file = data_file_out, row.names=FALSE, col.names = TRUE , quote=F)

## change owner of the file
system( paste0( "chown :nic_timpson_grp ", data_file_out) )

## change  permission of the file
system( paste0( "chmod 770 ", data_file_out) )

#########################
##
## 11.  END SESSION 
##
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

