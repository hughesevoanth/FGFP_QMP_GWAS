#######################################
##  A script to concatenate and format
##  SNPTEST OUTPUT
##  by: David Hughes
##  date: June 4th 2020
#######################################
echo -e "You are running a script to:\n\t1)concatenating SNPTEST results into one file\n\t2)reorder the columns to conform to a bed file format\n\t3)indexing the resultant file.\n"

## The file needs thre variables passed to it
## 1. the trait name
## 2. the SNPtest parent directory: there should be a folder in this directory holding the snptest results
## 3. the directory where the concatenated and indexed GWAS file should go
## 4. the directory where the top hits file should go

###
module add lib/htslib/1.10.2-gcc
###
echo -e "Adding the htslib library to your environment."
#######################################
## TRAIT name provided from command 
## line argument
#######################################
if [[ -z "$1" ]] ; then
    echo 'ERROR: you did not provide a trait name as argument 1'
    exit 1
fi

TRAIT=$1
###
echo -e "The trait you are processing is ${TRAIT}."

#######################################
## Directory containing SNPTEST results
## to process
#######################################
if [[ -z "$2" ]] ; then
    echo 'ERROR: you did not provide a SNPTEST results directory path as argument 2.'
    exit 1
fi

SNPTESTDIR=$2
###
echo -e "The directory of data you are processing is ${SNPTESTDIR}."

## **** MOVE INTO WORKING DATA DIRECTORY ****
cd ${SNPTESTDIR}${TRAIT}

#######################################
## Concatenated GWAS directory to place 
## processed/concatenated data
#######################################
if [[ -z "$3" ]] ; then
    echo 'ERROR: you did not provide a GWASDIR directory to place the results as argument 3.'
    exit 1
fi

GWASDIR=$3

GWASOUTFILE=${GWASDIR}${TRAIT}_allchr.txt

#######################################
## The top associated hits file
#######################################
if [[ -z "$4" ]] ; then
    echo 'ERROR: you did not provide a TOPHITS directory to place the associated SNPs as argument 4.'
    exit 1
fi

TOPHITSDIR=$4

TOPHITSFILE=${TOPHITSDIR}${TRAIT}_p5e8_hits.txt

### Give a little feedback to user
echo -e "Your new concatenated GWAS file will be ${GWASOUTFILE}."
echo -e "Your new associations at p<= 5e-8 will be ${TOPHITSFILE}."

#######################################
### Confirm that there are 22 chr files
#######################################
echo -e "Confirming that the data directory has 22 chromosome files."
##
chr_count=$(ls -1 | grep chr -c)

## WRAP THE FUNCITON IN AN IF STATEMENT
if [ $chr_count -ge 22 ]; then 
	###
	echo -e "\tCONFIRMED"
	#######################################
	### Identify the number of header lines
	#######################################
	## when do the column names start
	colnames_line=$(grep -n alternate_ids chr01.txt | cut -f1 -d:)
	## subtract one
	let data_starts_here=${colnames_line}+1

	#######################################
	### Concatenate together
	#######################################
	echo "Now concatenating all chromosome files together."
	#head -n ${colnames_line} chr01.txt > ${OUTFILE}
	head -n -1 chr01.txt | tail -n +${colnames_line} > temp.txt

	for i in {2..9} 
	do
		echo ${i} 
		head -n -1 chr0${i}.txt | tail -n +${data_starts_here} >> temp.txt
	done

	for i in {10..22} 
	do
		echo ${i} 
		head -n -1 chr${i}.txt | tail -n +${data_starts_here} >> temp.txt
	done

	#######################################
	## reorder columns
	#######################################	
	echo "Now reordering the columns."
	##
	#if[ "${TRAIT: -2}" == "AB" ]; then
		## CONT (AB) files
	#	awk 'BEGIN {FS=" "; OFS="\t"} {print $3, $4, $4, $2, $1, $5, $6, $14, $15, $16, $18, $19, $21, $9, $24, $25, $22}' temp.txt > ${GWASOUTFILE}
		### find ignificant sites SNPs
	#	awk '$17 <= 0.00000005' ${GWASOUTFILE} > ${TOPHITSFILE}
	#else	
		## PA files
		awk 'BEGIN {FS=" "; OFS="\t"} {print $3, $4, $4, $2, $1, $5, $6, $14, $15, $16, $18, $19, $20, $21, $22, $23, $24, $25, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $9, $47, $48, $45}' temp.txt > ${GWASOUTFILE}
		### find significant sites SNPs
		awk '$40 <= 0.00000005' ${GWASOUTFILE} > ${TOPHITSFILE}
	#fi
	###
	# remove trash
	rm temp.txt
	### BG zip the new data file
	echo "Now bgzipping the file."
	bgzip ${GWASOUTFILE}
	### index (tabix) the new data file
	echo "Now indexing the file."
	# tabix -p bed -S 1 ${OUTFILE}.gz
	tabix -c "#" -p bed -S 1 ${GWASOUTFILE}.gz
	### chown and chmod the new files
	echo "Now changing the group ownership and permissions."
	chown :nic_timpson_grp ${GWASOUTFILE}*
	chmod 770 ${GWASOUTFILE}*
	chown :nic_timpson_grp ${TOPHITSFILE}
	chmod 770 ${TOPHITSFILE}
	
	### move the data
	# echo "Now moving it to the declared results|output directory."
	# mv ${TRAIT}_allchr.* ${OUTPUTDIR}
	# mv ${TRAIT}_topsites.txt ${OUTPUTDIR}
	
	## TABIX extraction examples
	# zcat ${TRAIT}_allchr.txt.gz | head -n 1 > trpv6.txt
	# tabix ${TRAIT}_allchr.txt.gz 07:142568960-142583490 > trpv6.txt

	# zcat all_chrs.txt.gz | head -n 1 > trpv6.txt
	# tabix all_chrs.txt.gz 07:142568960-142583490 >> trpv6.txt

else "There are not 22 chromosomes in this directory"

fi
echo "FINISHED"

###################################
## first instance of chromosome 2
###################################
# grep -n -m 1 "02" all_chrs.txt
# cut -f 1 all_chrs.txt | grep -n -m 1 "02"

###################################
## print a range of lines to screen
###################################
# sed -n '593959,593979 p' all_chrs.txt | column -t |  less -S

###################################
## SNPTEST COLUMNS
###################################
# 1 alternate_ids 
# 2 rsid 
# 3 chromosome 
# 4 position 
# 5 alleleA 
# 6 alleleB 
# 7 index 
# 8 average_maximum_posterior_call 
# 9 info 
# 10 cohort_1_AA 
# 11 cohort_1_AB 
# 12 cohort_1_BB 
# 13 cohort_1_NULL 
# 14 all_AA 
# 15 all_AB 
# 16 all_BB 
# 17 all_NULL 
# 18 all_total 
# 19 all_maf 
# 20 missing_data_proportion 
# 21 cohort_1_hwe 
# 22 frequentist_add_pvalue 
# 23 frequentist_add_info 
# 24 frequentist_add_beta_1 
# 25 frequentist_add_se_1 
# 26 comment

