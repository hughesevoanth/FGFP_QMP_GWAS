#!/bin/bash

#PBS -N qmpGWAS
#PBS -l walltime=1:00:00 
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-7
#PBS -e runlog/error/
#PBS -o runlog/out/

## There are 5478 jobs (4224 in Apr 2021 run)
## source a paramater file holding argument paths
source ${pfile}

## array index from PBS -J above
i=${PBS_ARRAY_INDEX}

## extract the CHR to be run on this job
CHR=$(sed "${i}q;d" ${CONT_JOBS} | awk -F " " '{ print $1 }')

## extract the trait to be run on this job
TRAIT=$(sed "${i}q;d" ${CONT_JOBS} | awk -F " " '{ print $2 }') 

## give some feedback
echo "Processing chromsomes " ${CHR} " and trait " ${TRAIT}

## full path to the bgen file to use
GENDIR=${data_processed_directory}genotypes/bgen_gwas/
GEN=${GENDIR}data_chr${CHR}.bgen

## full path to the sample file
SAM=${data_processed_directory}samplefile/fgfp_qmp_sample_ordered.sample

## define the OUTDIRECTORY
OUTDIR=${results_directory}${TRAIT}

##  Make a directory for the trait being run
mkdir -p ${OUTDIR} # -p means dont remake it if it exists
mkdir -p ${OUTDIR}/log

## change ownership and permissions
chown :nic_timpson_grp ${OUTDIR}
chmod 770 ${OUTDIR}

chown :nic_timpson_grp ${OUTDIR}/log
chmod 770 ${OUTDIR}/log

FILEOUT=${OUTDIR}/chr${CHR}.txt
LOGFILEOUT=${OUTDIR}/log/chr${CHR}.log

### load my snptest module
module add apps/snptest/2.5.4

## run SNPTEST
time snptest_v2.5.4-beta3 \
	-data ${GEN} ${SAM} \
	-exclude_samples ${EX} \
	-o ${FILEOUT} \
	-log ${LOGFILEOUT} \
	-missing_code NA \
	-pheno ${TRAIT} \
	-use_raw_phenotypes \
	-frequentist 1 \
	-method expected \
	-hwe \
        -lower_sample_limit 100 \
        -analysis_name "FGFP QMP GWAS by David Hughes and Malte Rühlemann"


## change owner and permissions
chown :nic_timpson_grp ${FILEOUT}
chown :nic_timpson_grp ${LOGFILEOUT}
chmod 770 ${FILEOUT}
chmod 770 ${LOGFILEOUT}
