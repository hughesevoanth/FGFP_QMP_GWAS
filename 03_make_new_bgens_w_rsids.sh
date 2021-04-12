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

## PATHS for data in
datadir=${data_input_directory}genotypes/bgen/
mapdir=${data_input_directory}genotypes/bgen_rsid_map/

## PATH for data out
	## just make sure there is a place to put the new bgen files
mkdir -p ${data_input_directory}genotypes/bgen_w_rsids
outdir=${data_input_directory}genotypes/bgen_w_rsids/

## sample file
SAM=${datadir}data.imputed.sample

#########################
## RUN qctools !
#########################
# ls -1 ${datadir} | grep .bgen | while read file 
# do 
# 	chr=${file:8:2}
# 	in=${datadir}/${file}
# 	map=${mapdir}/${file%.bgen}_rsid_map.txt
# 	out=${outdir}/${file}
# 	qctool -g ${in} \
# 	-s ${SAM} \
# 	-assume-chromosome ${chr#0} \
# 	-og ${out} \
# 	-map-id-data ${map} 
# done

in=${datadir}/data_chr${CHR}.bgen
map=${mapdir}/data_chr${CHR}_rsid_map.txt
out=${outdir}/data_chr${CHR}.bgen
qctool -g ${in} \
-s ${SAM} \
-assume-chromosome ${CHR#0} \
-og ${out} \
-map-id-data ${map}


## insure that ownership and permissions are OK
chown -R :nic_timpson_grp ${outdir}
chmod -R 770 ${outdir}