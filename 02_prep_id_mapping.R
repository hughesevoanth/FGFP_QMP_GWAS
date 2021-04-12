##############################################################################
## Script to: 
##  1) prepare an -map-id-data file for the addition of rsids
##     using QCtools.
##
## By: David Hughes, Caroline Bull, Laura Corbin
## Date begun: Friday March 28th 2020
##
##############################################################################
library(tidyverse)
library(data.table)

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

#######################################
##
## (III) Data Directories & Files
##
#######################################
## path to data directories
snpstats_dir = paste0(data_input_dir, "genotypes/bgen_snpstats/")
rsid_dir = paste0(data_input_dir, "genotypes/bgen_rsids/")

## make a new directory to place the data
map_dir = paste0(data_input_dir, "genotypes/bgen_rsid_map")
n = paste0("mkdir -p ", map_dir)
system(n)
map_dir = paste0(data_input_dir, "genotypes/bgen_rsid_map/")

## list of files to process
snpstats_files = list.files(snpstats_dir)
k = grep("_snp-stats.txt", snpstats_files)
snpstats_files = snpstats_files[k]

rsid_files = list.files(rsid_dir)
k = grep("hg19_snp138_dropped", rsid_files)
rsid_files = rsid_files[k]

## a check to be sure all of the data you need is present
if( length(rsid_files) != length(snpstats_files) ){
	stop("The number of snp-stats files does not match the number of rsid files.\nCheck the data in your snpstats and rsid directories to be sure it is complete.")
} else {
	cat( paste0( "Looking good: you have ", length(rsid_files), " rsid files and ", length(snpstats_files), " snp-stats files.\n\tWe will proceed to process them\n") )
}

#######################################
##
## (IV) Read in, match and produce map files
##
#######################################
for(i in 1:length(snpstats_files)){

	cat(paste0("Proccessing chromosome ", i, ".\n"))

	## read in snpstats file
	n = paste0(snpstats_dir, snpstats_files[i] )
	ss = fread(n, header = TRUE, skip = 10)

	## read in rsid file
	n = paste0(rsid_dir, rsid_files[i] )
	rs = fread(n, header = FALSE, sep = "\t")
	colnames(rs) = c("build","rsid","anno")

	## make a snpid in the rsid file for matching
	snpid = sapply(rs$anno, function(x){
		a = strsplit( as.character(x), split = " ")[[1]][c(1,2,4,5)]
		o1 = paste(a[1:2], collapse = ":")
		out = paste( c(o1, a[3:4]), collapse = "_")
		})

	## add to rs file
	rs$snpid = snpid

	## now match SNPID between ss and rs data tables
	##   NOTE: the rsid column in the SNP-STATS tables are SNPIDs
	m = match(ss$rsid, rs$snpid )

	## now make a mapid file
	## 6 original columns [snpid, rsid, chr, position, a1, a2] + 6 NEW columns [snpid, rsid, chr, position, a1, a2] 
	map = cbind( ss[, c(1,2,3,4,5,6)], ss[,2], rs[m,2], ss[, c(3,4,5,6) ] )

	n = c("snpid", "rsid", "chr", "position", "a1", "a2")
	colnames(map) = c( paste0("old_", n),  paste0("new_", n) )

	## replace any rsid NAs with the old snpid
	w = which( is.na( map$new_rsid ) )
	map$new_rsid[w] = map$new_snpid[w]

	## write data to file
	cat(paste0("\tWriting map file for chromosome ", i, ".\n"))

	newfilename =  paste0( gsub("_snp-stats.txt","", snpstats_files[i]), "_rsid_map.txt" )
	outname = paste0(map_dir, newfilename)
	write.table(map, file = outname, row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
}


cat(paste0("\n\n\t** DONE!! **\n"))








