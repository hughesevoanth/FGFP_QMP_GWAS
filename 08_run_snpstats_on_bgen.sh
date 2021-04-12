##############################################################################
## Script to: 
##  1) Run SumStats on New Bgen with rsids
##  2) use only those samples that will be in the GWAS
##	By David Hughes, Dec 5th 2020
## 			Built off script By: Laura Corbin, Caroline Bull, David Hughes
## 					Date begun: Friday March 13th 2020
##
##############################################################################
# source parameter_files/pfile.txt
source $1

## Define the chromosome to work with
CHR=$2
if ((${#CHR} == 1))
then
	CHR=0$2
fi

## add needed module 
module add apps/qctool/2.0.7

## DIRECTORY WITH THE BGEN with rsids 
datadir=${data_input_directory}genotypes/bgen_w_rsids/

## DIRECTORY WHERE THE SNPSTATS FIELS SHOULD GO
##  just make sure there is a bgen_snpstats directory
mkdir -p ${data_processed_directory}genotypes/bgen_snpstats
outdir=${data_processed_directory}genotypes/bgen_snpstats/

## SAMPLE FILE FOR THE SOURCE BGEN FILES
SAMIN=${data_input_directory}genotypes/bgen/data.imputed.sample

## THE GWAS SAMPLE FILE (REDUCED NUMBER OF INDIVIDUALS FROM SOURCE BGEN)
## NOT USEING THIS HERE
gwasfile=${data_processed_directory}samplefile/fgfp_qmp.sample

## SAMPLES IN THE GWAS, I.E. THE SAMPLES TO USE TO ESTIMATE SNPSTATS
SAM2KEEP=${data_processed_directory}samplefile/samples_in_gwas.txt


## insure that ownership and permissions are OK
chown :nic_timpson_grp -R ${outdir}
chmod 770 -R ${outdir}

#########################
## RUN qctools !
#########################
### LOOP
# ls -1 ${datadir} | grep .bgen | while read file
# do
# 	chr=${file:8:2}
# 	qctool -g ${datadir}${file} -s ${SAM} -assume-chromosome ${chr#0} -snp-stats -incl-samples ${samingwas} -osnp ${outdir}${file%.bgen}_snp-stats.txt 
# done

qctool -g ${datadir}data_chr${CHR}.bgen -s ${SAMIN} -assume-chromosome ${CHR#0} -snp-stats -incl-samples ${SAM2KEEP} -osnp ${outdir}data_chr${CHR}.snpstats


## insure that ownership and permissions are OK
chown :nic_timpson_grp ${outdir}data_chr${CHR}.snpstats
chmod 770 ${outdir}data_chr${CHR}.snpstats

