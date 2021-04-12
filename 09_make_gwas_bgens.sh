###################################
## Produce Study BGENS
## 1) extract SNPs from the SNPSTATS file to exclude
## 2) produce a bgen with only study samples and SNPs with MAC > 10
##
## By: David Hughes; Date: dec 5th 2020
###################################
## Read in paramater file
source $1

## array id
## Define the chromosome number
CHR=$2
if ((${#CHR} == 1))
then
        CHR=0$2
fi

## define paths to input and output
IN=${data_processed_directory}genotypes/bgen_snpstats/data_chr${CHR}.snpstats
OUT=${data_processed_directory}genotypes/bgen_snpstats/snps_2qc_remove/data_chr${CHR}.snps2remove

mkdir -p ${data_processed_directory}genotypes/bgen_snpstats/snps_2qc_remove

## the MINOR ALL FREQ
# MAC=10
# N=1666
# MAF = MAC/N
# MAF=0.00600240096

## COLUMNS
## HWE = 9
## maf = 14
## info =  17
tail -n +13 ${IN} | awk -F "\t" '{ if(($9 <= 0.000000000001) || ($14 < 0.00600240096) || ($17 < 0.3)) { print } }' |  cut -f 2 > ${OUT}


###################################
##  Now generate the new
##  Bgen files
##
###################################

## add QCtools to environment
module add apps/qctool/2.0.7

## define paths to input and output
BGEN=${data_input_directory}genotypes/bgen_w_rsids/data_chr${CHR}.bgen
SAMPLE=${data_input_directory}genotypes/bgen/data.imputed.sample
##
mkdir -p ${data_processed_directory}genotypes/bgen_gwas
OUTBGEN=${data_processed_directory}genotypes/bgen_gwas/data_chr${CHR}.bgen
OUTSAMPLE=${data_processed_directory}genotypes/bgen_gwas/data_chr${CHR}.sample
##
SNP2EXCLUDE=${OUT}
SAM2INCLUDE=${data_processed_directory}/samplefile/samples_in_gwas.txt

time qctool -g ${BGEN} -s ${SAMPLE} -og ${OUTBGEN} -os ${OUTSAMPLE} -excl-rsids ${SNP2EXCLUDE} -incl-samples ${SAM2INCLUDE}


