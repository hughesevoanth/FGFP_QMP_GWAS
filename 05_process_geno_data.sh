#######################################################
## script to generate PCs
## PCs then added to .sample file
## by: Laura Corbin, Caroline Bull, and David Hughes
## date: Feb 12th 2020
#######################################################
## this script needs to be run from scripts folder

####################
##  SET UP 
####################
#  module add apps/plink/2.0.0

# call file with file paths defined
# source parameter_files/pfile.txt
source $1

# move to working directory
cd $workingdirectory


###################################
##  GENERATE PCS IN PLINK 
###################################

genfile=${data_input_directory}genotypes/plink/data

## make the necessary directory folders
## and out shortcut names 
mkdir -p ${data_processed_directory}genotypes
mkdir -p ${data_processed_directory}genotypes/plink
mkdir -p ${data_processed_directory}genotypes/plink/pruned_data
outfile1=${data_processed_directory}genotypes/plink/pruned_data/data

mkdir -p ${data_processed_directory}genotypes/plink/pc_data
outfile2=${data_processed_directory}genotypes/plink/pc_data/data

## insure that ownership and permissions are OK
chown :nic_timpson_grp -R ${data_processed_directory}/genotypes
chmod 770 -R ${data_processed_directory}/genotypes

###############################################
## Generate a list qc'd geno IDs 
###############################################
##  - use Sample_Exclusion_Criteria_out.txt from FGFP 2017-03-01 data release. Each "1" in the matrix identifies an individual who has been flagged for a given QC step (N = 2165 excluding all flagged)
# awk '{ if(($2 == 0) && ($3 == 0) && ($4 == 0) && ($5 == 0) && ($6 == 0) && ($7 == 0)) { print $1, $1 } }' $data_input_dir/genetics/Sample_Exclusion_Criteria_out.txt > $data_input_dir/genetics/qcgenokeep.txt
# sed -i 1,2d $data_input_dir/genetics/qcgenokeep.txt
# cut -d " " -f 1-10 $data_input_dir/genetics/qcgenokeep.txt | head

## Using the sample file, which already contains the genotype IDs for individuals with Metabolon and Genotype data of quality
awk '{ print $3, $3}' ${data_processed_directory}qmp/fgfp_qmp_gwas_dataset_step1_raw.txt | tail -n +2 > ${data_processed_directory}qmp/genokeep.txt

## insure that ownership and permissions are OK
chown :nic_timpson_grp ${data_processed_directory}qmp/genokeep.txt
chmod 770 ${data_processed_directory}qmp/genokeep.txt

###############################################
## Identify and extract a genotype data set of 
## LD pruned SNPS (pairwise LD pruning)
###############################################

## Exclude high LD regions ** FROM GIB gist.github.com/explodecomputer/e4438771d04534e058ad#file-highldregionsb37-awk **
## Add statement to restrict PC generation to those with metabolomic data - defined in 01a_merge_pheno.R (exclusion_list_doubleid.txt)
# awk -f $data_input_dir/genetics/highldregionsb37.awk $data_input_dir/genetics/plink/data.bim > $data_input_dir/genetics/highldregionsb37.txt

## generate LD pruned data set
plink2 \
--bfile $genfile \
--exclude 'range' ${data_input_directory}genotypes/LongRange_LD_b37_20200212.txt \
--keep ${data_processed_directory}qmp/genokeep.txt \
--indep-pairwise 10000 5 0.1 \
--autosome \
--maf 0.05 \
--geno 0.01 \
--hwe 1e-05 \
--make-bed \
--out $outfile1


## Generate PCs
plink2 \
--bfile $outfile1 \
--pca \
--autosome \
--out $outfile2 \

## change ownership and permissions
chown :nic_timpson_grp $data_processed_directory/genotypes/plink/pruned_data/*
chmod -R 770 $data_processed_directory/genotypes/plink/pc_data


chown :nic_timpson_grp $data_processed_directory/genotypes/plink/pruned_data/*
chmod -R 770 $data_processed_directory/genotypes/plink/pruned_data
