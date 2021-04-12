#!/bin/sh

#PBS -N metabGWAS_dataprep_pipeline
#PBS -l select=1:ncpus=4:mem=5gb
#PBS -l walltime=36:00:00
#PBS -o ../runlog/run_pipeline.log
#PBS -e ../runlog/run_pipeline.log
#PBS -j oe

### loaded needed modules
module add lang/r/3.6.0-gcc
module add apps/plink/2.0.0

## define my working directory
wdir=$WORK/fgfp/scripts/dh/fgfp_metabolomics_gwas

## step1 collate all of my data metabolites, clinical, and
##       GWAS availability
time Rscript $wdir/01_data_collation.R $wdir/parameter_files/pfile.txt
echo -e "\n\t** End of script 01_data_collation.R **\n"

## step2 process the QC'd genotyped data into a data set to
##       estimate principal components 
time sh $wdir/02_process_geno_data.sh $wdir/parameter_files/pfile.txt
echo -e "\n\t** End of script 02_process_geno_data.sh **\n"

## step3 transform metabolite data according to (1) Claudia Langenberg's protocols
##       and (2) the Team Timson best practices
time Rscript $wdir/03_data_transform.R $wdir/parameter_files/pfile.txt
echo -e "\n\t** End of script 03_data_transformation.R **\n"

## step4 sex and age imputation with KNN
##       and write to file
time Rscript $wdir/04_covariate_knn_impute.R $wdir/parameter_files/pfile.txt
echo -e "\n\t** End of script 04_covariate_knn_impute.R **\n"

## step5 fit a linear model to all continuous traits
##       and write to file
time Rscript $wdir/05_data_lm_residuals.R $wdir/parameter_files/pfile.txt
echo -e "\n\t** End of script 05_data_lm_residuals.R **\n"

## step6 summary statistics and plotting data distributions
##       and write to file
time Rscript $wdir/06_sumstats_and_plots.R $wdir/parameter_files/pfile.txt
echo -e "\n\t** End of script 06_sumstats_and_plots.R **\n"

## step7 estimate -snpstats on the raw|source bgen files with QCtools
# sh $wdir/07_run_snpstats_on_bgen.sh $wdir/parameter_files/pfile.txt
# echo -e "\n\t** End of script 07_run_snpstats_on_bgen.sh **\n"

## step8 extract the SNPIDs for SNPs to filter or remove from the bgen files
##       this will reduce the number of uneccesarry analysis performed in the GWAS
# Rscript $wdir/08_extract_SNPids_2_filter.R $wdir/parameter_files/pfile.txt
# echo -e "\n\t** End of script 08_extract_SNPids_2_filter.R **\n"

## step09 get rsids with for our SNPs using ANNOVAR
##        this uses the source SNP-STATS files to extract mapping info
# sh $wdir/09_annotate_bgen_snpIDs.sh $wdir/parameter_files/pfile.txt
# echo -e "\n\t** End of script 09_annotate_bgen_snpIDs.sh **\n"

## step10 Make rsid map files
##       to add rsid to the source bgen files
#time Rscript $wdir/10_prep_id_mapping.R $wdir/parameter_files/pfile.txt
#echo -e "\n\t** End of script 10_prep_id_mapping.R **\n"

## step11 make new SOURCE bgen files with rsids
# sh $wdir/11_make_new_bgens_w_rsids.sh $wdir/parameter_files/pfile.txt
# echo -e "\n\t** End of script 11_make_new_bgens_w_rsids.sh **\n"

## step12 make STUDY bgen files
# sh $wdir/12_produce_study_bgen.sh $wdir/parameter_files/pfile.txt
# echo -e "\n\t** End of script 12_produce_study_bgen.sh **\n"



