##########################################
##  what: A script to extract rsids
##        using ANNOVAR
##  by: David A Hughes, Caroline Bull 
##      and Laura Corbin
### date: Mar 27th 2020
## 
##########################################
# source parameter_files/pfile.txt
source $1

## Define the chromosome to work with
CHR=$2
if ((${#CHR} == 1))
then
	CHR=0$2
fi

## load the perl module into the environment
module add apps/qctool/2.0.7
module add lang/perl/5.30.0-bioperl-gcc

###########################
## (i) run SNPstats on bgen
###########################
## make a new directory if not there
mkdir -p ${data_input_directory}genotypes/bgen_snpstats

datadir=${data_input_directory}genotypes/bgen/
outdir=${data_input_directory}genotypes/bgen_snpstats/

#########################
## RUN qctools !
#########################
FILE=data_chr${CHR}.bgen
qctool -g ${datadir}${FILE} -s ${datadir}data.imputed.sample -assume-chromosome ${CHR#0} -snp-stats -osnp ${outdir}${FILE%.bgen}_snp-stats.txt 


###########################
## (ii) prepare data for ANNOVAR
###########################

## make a new directory to place the annovar input files into
mkdir -p ${data_input_directory}genotypes/bgen_snpstats/annovar_input
mkdir -p ${data_input_directory}genotypes/bgen_rsids

## extract SNP information from the SNPSTATS files
## we need: chromosome start stop allele1 allele2

in=${data_input_directory}genotypes/bgen_snpstats/${FILE%.bgen}_snp-stats.txt
out=${data_input_directory}genotypes/bgen_snpstats/annovar_input/chr${CHR}_annovarinput.txt
tail -n +12 ${in} | cat | awk '{ OFS=" "; print $3, $4, $4, $5, $6 }' > ${out}

## insure that ownership and permissions are OK
chown -R :nic_timpson_grp ${data_input_directory}genotypes/bgen_snpstats/annovar_input
chmod -R 770 ${data_input_directory}genotypes/bgen_snpstats/annovar_input


###########################
## (iii) now run ANNOVAR on each 
##      chromosome of data
###########################

## move to proper directory --- the ANNOVAR directory
cd ${workingdirectory}../../annovar/
##  iterate over the chr split annovar input files
in=${data_input_directory}genotypes/bgen_snpstats/annovar_input/chr${CHR}_annovarinput.txt
out=${data_input_directory}genotypes/bgen_rsids/chr${CHR}
perl annotate_variation.pl -filter \
	-out ${out} \
	-build hg19 \
	-dbtype snp138 \
	${in} \
	humandb


## insure that ownership and permissions are OK
chown -R :nic_timpson_grp ${data_input_directory}genotypes/bgen_rsids
chmod -R 770 ${data_input_directory}genotypes/bgen_rsids





