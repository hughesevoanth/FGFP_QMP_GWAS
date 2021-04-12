##############################################################################
## Script to: 
##  1) read in QMP Data
##  2) identify taxa that need to go through a hurdle model
##  3) perform the data transformations
## By: David A Hughes Date: Dec 5th
## 		Modified script from By: Laura Corbin, Caroline Bull, David Hughes
## 								Date begun: Feb 10 2020
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
data_file_in = paste0(data_processed_dir,"qmp/fgfp_qmp_gwas_dataset_step1_raw.txt")
data_file_out = paste0(data_processed_dir,"qmp/fgfp_qmp_gwas_dataset_step2_transformed.txt")

#######################################
##
## (III) source our FUNCTIONS file
## 	
#######################################
source( paste0(working_dir, "scripts/fgfp_qmp_gwas/FUNCTIONS.R"))


#######################################
##
## (IV) Start a log file
##
#######################################
cat( paste0("Starting a log file - it will be in your processed data metabolites directory.\n") )


logfilename = paste0(data_processed_dir, "qmp/", project, "_06_data_transformations_", today,  "_logfile.txt")
sink(file = logfilename , split = TRUE  )

## send information to log file
cat( paste0("Transforming QMP microbiome data.\n") )
cat( paste0("  -- Todays date is: ", today, ".\n\n") )

#######################################
##
## (V) Read and merge the data
##
#######################################

########################
## (Va) Read in QMP data 
########################
cat( paste0("I) Reading in QMP Data.\n") )
#### read in the data
mdata = read.table( data_file_in, header = TRUE, as.is = TRUE)
rownames(mdata) = mdata$id

## feedback to user
prefix = c("P_","C_","O_","F_", "G_")
mt_locations = sapply(prefix, function(x){
	ss = substr(colnames(mdata), 1, nchar(x))
	w = grep(x, ss)
	return(w)
	})
mt_locations = unlist(mt_locations)
####
f_count = length( mt_locations )
cat(paste0("     a) There are ", f_count, " QMP features in the QMP & clinical data file \n"))
cat(paste0("     b) There are ", nrow(mdata), " samples in the metabolomics & clinical data file.\n"))

########################
## (Vb) Read PCs
########################
cat( paste0("II) Reading in principal component (PC) data.\n") )
####
n = paste0(data_processed_dir,"genotypes/plink/pc_data/data.eigenvec")
pcdata = read.table(n, header = FALSE, as.is = TRUE)
colnames(pcdata) = c("GENOTYPE_IID","GENOTYPE_FID", paste0("PC", 1:( ncol(pcdata)-2 ) )  )

## feedback to user
cat(paste0("     a) There are ", ncol(pcdata)-2, " PCs in the PC data file.\n"))
cat(paste0("     b) There are ", nrow(pcdata), " samples in the PC data file.\n"))


########################
## (Vc) merge PC and MT data
########################
cat( paste0("III) Merge Metabolite & clinical data with the PC data.\n") )
####
m = match(mdata$genetic_id, pcdata$GENOTYPE_IID)
mdata = cbind(mdata[, -c(mt_locations) ], pcdata[m, -c(1:2) ], mdata[, mt_locations ])

#######################################
##
## (VI) Exclude individuals PC outliers
##    
#######################################
cat( paste0("IV) Exclude individuals based on PCs.\n") )
####
cuttoffs = apply(mdata[, c("PC1","PC2")], 2, function(x){
	m = mean(x, na.rm = TRUE)
	sd = sd(x, na.rm = TRUE) * 5
	cuttoffs = c( m - sd, m + sd)
	})

### Identify any outliers and remove
indtoremove = which( mdata$PC1 < cuttoffs[1,1] | mdata$PC1 > cuttoffs[2,1] |  mdata$PC2 < cuttoffs[1,2] | mdata$PC2 > cuttoffs[2,2] )
cat( paste0("     a) Screen for PC outliers: +/- 5SD of PC1  | PC2 mean.\n") )
cat( paste0("        i) There are ", length(indtoremove), " individuals that are PC1-2 outliers.\n") )
## if we need to remove some individuals
if(length(indtoremove) > 0){
	 mdata = mdata[ -indtoremove , ]
	 cat( paste0("           ** ", length(indtoremove), " individuals have been removed from the meatbolite & clinical data set.\n") )
}

#### plot PCA 
n = paste0(data_processed_dir, "genotypes/plink/pc_data/PC1-2.pdf")
cat( paste0("     b) PC1-2 have been plotted and printed to file\n        - ",n,".\n") )
####
pdf(file = n, width = 6, height =6)
plot(mdata$PC1, mdata$PC2, pch = 21, 
	col = "steelblue", cex = 1.25,
	main = "FGFP metabolite GWAS PCA", 
	xlab = "PC1", ylab = "PC2")
abline(v = cuttoffs[, 1], col = "red", lwd = 3 )
abline(h = cuttoffs[, 2], col = "red", lwd = 3 )
dev.off()

#######################################
##
## (VII) Identify taxa to retain for GWAS 
##
#######################################
cat( paste0("V) Identify taxa to retain for GWAS. Those with at least 100 observations.\n") )

## redefine MT locations; as PC addition changes indexes
prefix = c("P_","C_","O_","F_", "G_")
mt_locations = sapply(prefix, function(x){
	ss = substr(colnames(mdata), 1, nchar(x))
	w = grep(x, ss)
	return(w)
	})
mt_locations = unlist(mt_locations)

## Number of samples with non-zero data
obs_count = apply( mdata[, mt_locations] , 2, function(x){ sum(x != 0) })
## taxa counts
taxa2keep = names( which(obs_count >= 100) )
taxa2remove = names( which(obs_count < 100) )

cat( paste0("\ta) A total of  ", length(taxa2keep), " MTs have >= 100 non-zero observations and will be kept in the GWAS.\n") )
cat( paste0("\tb) A total of  ", length(taxa2remove), " MTs have < 100 non-zero observations and will be discarded.\n") )

# w = which(colnames(mdata) %in% taxa2remove)
# mdata = mdata[, -w]


#######################################
##
## (VIII) Build the GWAS traits
##   1. Identify those taxa that need a hurdle model
##
#######################################
cat( paste0("VI) Identify taxa that need a hurdle model and build GWAS traits.\n") )
cat( paste0("\t- Hurdle taxa are defined as those with >=100 zero counts and at least >=100 presence counts.\n") )
cat( paste0("\t- Continuous taxa are defined as those with < 50% zero counts (", nrow(mdata)/2 , ").\n") )
cat( paste0("\t- Continuous and truncated-continuous MTs will be ranknormal transformed.\n") )
## define a new data set
d = mdata[, taxa2keep]
cat( paste0("   a) Starting with ", nrow(d), " samples and ", ncol(d) ," QMP MTs\n") )
###
gwas_qmp_data = c()
for(i in 1:ncol(d)){
  taxa_name = colnames(d)[i]
  x = d[,i]
  #############
  ## PA (Hurdle) trait
  #############
  pa = x
  pa[pa>0] = 1
  pcount = sum(pa==1, na.rm = TRUE)
  acount = sum(pa==0, na.rm = TRUE)
  PAtrait = ifelse( min( c(pcount, acount) )>=100  , TRUE , FALSE)
  if(PAtrait == TRUE){
  	##########################
  	## presence|absence
  	##########################
  	percent_absent = round(acount/nrow(d), d=2)
  	###
    gwas_qmp_data = cbind(gwas_qmp_data, pa)
    n = paste0(taxa_name,"_",percent_absent,"_PA")
    colnames(gwas_qmp_data)[ncol(gwas_qmp_data)] = n
    ##########################
  	## truncated continuous
  	##########################
  	truncated_continuous = x
    truncated_continuous[truncated_continuous == 0] = NA
    truncated_continuous = rntransform(truncated_continuous)
    #
    gwas_qmp_data = cbind(gwas_qmp_data, truncated_continuous)
    n = paste0(taxa_name,"_",percent_absent,"_truncAB")
    colnames(gwas_qmp_data)[ncol(gwas_qmp_data)] = n
    }
  ####################
  ## Continuous
  ##  trait
  ###################
  ## continuous trait with no truncation
  if(acount < nrow(d)/2 ){
  	percent_absent = round(acount/nrow(d), d=2)
  	continuous = rntransform(x)
  	gwas_qmp_data = cbind(gwas_qmp_data, continuous )
  	n = paste0(taxa_name,"_",percent_absent,"_AB")
  	colnames(gwas_qmp_data)[ncol(gwas_qmp_data)] = n
	}
}
## Add rownames
rownames(gwas_qmp_data) = rownames(d)

########################################
## how may AB, truncAB and PA traits do
## we have?
########################################
ABcount = length( grep("_AB", colnames(gwas_qmp_data) ) )
truncABcount = length( grep("_truncAB", colnames(gwas_qmp_data) ) )
PAcount = length( grep("_PA", colnames(gwas_qmp_data) ) )

cat( paste0("   b) We built a total of ", ABcount, " abundance MTs, ", truncABcount ," truncated AB MTs, and ", PAcount ," presence|absence MTs\n") )


#######################################
##
## (VIII - B) Transform cell_counts_per_g 
##            data
#######################################
mdata$cell_counts_per_g = rntransform(mdata$cell_counts_per_g)

#######################################
##
## (IX) Merge in CLINIC DATA 
##
#######################################
cat( paste0("VII) Merge new transformations QMP MTs back with clincial/covariate data.\n") )

dataout = cbind( mdata[, -mt_locations],  gwas_qmp_data)


#######################################
##
## (X) WRITE DATA TO FILE
##
#######################################

cat( paste0("VIII) Writing metabolon data set a .txt file. \n") )

## write generated data set to flat text file
write.table(dataout, file = data_file_out , row.names = TRUE, col.names = TRUE , quote=F)


## change owner of the file
system( paste0( "chown :nic_timpson_grp ", data_file_out) )

## change  permission of the file
system( paste0( "chmod 770 ", data_file_out) )


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

