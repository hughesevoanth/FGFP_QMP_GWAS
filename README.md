# FGFP QMP GWAS
	by: David Hughes
	date: Dec 2020
	project portal: wt1/wp2/004

## Description

This git repository holds the scripts used to set up and run the FGFP quantitative microbiome profiling (QMP) GWAS. 

## Workflow

	1.  sh:	run ANNOVAR to annotate SNPs, with just SNPIDs, to also have rsids
	2.  R:		produce a mapping file to add newly annotated rsids to bgens
	3.  sh:	produce new bgens with rsids
	4.  R:		collate all of the data needed for the GWAS and make a list of individuals that will be included
	5.  sh: 	used the hard call array plink genotype data and produce new PCs
	6.  R:		identify abundance phenotypes to run in GWAS and transform that data into microbial traits (MTs)
	7.  R:		run linear models against covariates and extract residuals & produce the study sample-file
	8.  sh:	run snp-stats on the GWAS samples
	9.  sh:	make GWAS/study bgens excluding SNPs and individuals not included
	10. R:		order the study sample file to match the order of samples in the bgen
	11. sh:	run the continuous (AB and truncAB) trait GWAS
	12. sh:	run the presence|absence (PA) GWAS
	13. sh:	format the SNPTEST results: concatenate chromosomes, format into a BED file, index for tabix, extract SNPs with P>5e-8.
	 

## Directory Tree Structure

	fgfp
	|--QMP_GWAS
		|--scripts
		|	|--fgfp_qmp_gwas			** this git repo **
		|	   |--parameter_files
		|		`--runlog
		|--results
		|  	|--snptest_results		** gwas output goes here **
		|	|--concat_gwas		** concatenated version of GWAS here **
		|	`--tophits
		`--data
		   |--source
		   |	|--clinical		** study covariates **
		   |	|--genotypes
		   |	|   |--plink		
		   |	|   |--bgen
		   |	|   |--bgen_rsid_map		** map files to make new bgens **
		   |	|   |--bgen_rsids		** ANNOVAR output ** 
		   |	|   |--bgen_w_rsids		** new bgen files with rsids-for ALL samples **
		   |	|   |--LongRange_LD_b37_20200212.txt
		   |	|   |--README
		   |    |   `--text_bgens
		   |	|--linker
		   |    '--qmp		** source QMP files **
		   `--processed
				|--genotypes
				|   |--bgen_gwas		** bgen files used in QMP GWAS is here **
				|   |--bgen_snpstats
				|   `--plink
				|		|--pc_data		** estimated PCs for study is here **
				|		`--pruned_data		** data used to generate PCs **
				|--qmp		** each step of QMP data (MT) processing here **
				`--samplefile		** the QMP GWAS study sample file **
				

## How to run the pipeline prep on BP
	1. run the pipeline
		> qsub run_pipeline.sh
	
	5. Alternatively, you can edit run_pipeline.sh to remove the "wdir" paramater line and run the script as
		> qsub -v wdir=fullpath/to/workingdir/ run_pipeline.sh


## How to submit GWAS job array on BP

### submit continuous (rnt) trait jobs 
	qsub -v pfile=fullpath/scripts/fgfp_qmp_gwas/parameters/pfile.txt 11_run_continuous_GWAS.sh

### submit presence|absence (PA) trait jobs
	qsub -v pfile=fullpath/scripts/fgfp_qmp_gwas/parameters/pfile.txt 12_run_PA_GWAS.sh


## How to use tabix on GWAS results file

	The results for each GWAS will be found in a single BED formatted bgzip file.
	All chromosomes have been concatenated together in this one file and indexed with TABIX.
	You can extract SNPs from this file using a range extraction. As such, you will have to know the mapping data (chr and bp) of your SNP to extract it using TABIX. 
	This can be done as such:

	zcat ${TRAIT}_allchr.txt.gz | head -n 1 > trpv6.txt
	tabix ${TRAIT}_allchr.txt.gz 07:142568960-142583490 >> trpv6.txt

## Pipeline Descriptions

I) ---- 01_annotate_bgen_snpIDs.sh ----

	1) Run snp-stats, using qctools, on each raw saource bgen file. 
	2) use the snp-stats file to build the input file for ANNOVAR
		a) this includes chromosome, start, end, A1, A2
	3) then run ANNOVA
	
	- This scripts runs once for each chromosome
	- a job array was set up to run theses jobs for each chromosome. 

II) ---- 02_prep_id_mapping.R ----

	1) Use ANNOAVAR (hg19_snp138_dropped) output
		- produce a snp annotation map file for qctools

III) ---- 03_make_new_bgens_w_rsids.sh ----

	1) Run qctools with the map files made in step two to make new bgens that include rsids.

IV) ---- 04_data_collation.R ----

	1) Pull together phenotype data for GWAS
		- FGFP linker file
			- includes a binary for genotype availability
		- QMP source file
		- QMP batch data
		- clinic/covariate data
	2) Merge data into a single file
	3) Filter down to those with QMP and genotype data
		- remove those with sex and snpsex mismatches

V) ---- 05_process_geno_data.sh ----

	1) Use hard call array data used for imputation to:
		- produce LD pruned plink file
		- estimate PCs for study

VI) ---- 06_data_transform.R ----

	1) Read in the data from step 4
	2) Read in the PC data from step 5
	3) merge phenotype and PC data
	4) Sample exclusion based on PC1 and PC2
		- +/- 5SD from mean on PC1 and PC2
			- NO exclusions in QMP GWAS (this study) made here
	5) Identify taxa to retain for GWAS
		- If >= 100 non-zero observations made for a taxa it was retained.
	6) Identify Hurdle (truncAB + PA) taxa.
		- those where: PRESENCE count >= 100 & ABSENCE count >= 100
	7) Identify abundance (AB) taxa.
		- those where: ABSENCE count < half the sample size
			- NOTE: this means that up to 50% of the data may be zero counts.
				-  zero counts that WILL BE RANDOMLY RANKED.
	8) Also transformed "cell_counts_per_g" as a GWAS phenotype
	9) All AB and truncAB traits were rank normal transformed with tied values randomly split.

VII) ---- 07_data_lm_residuals.R ----

	1) Abundance traits fit to linear model
		- MT ~ age + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
	2) Residuals extracted to form the GWAS phenotype
	3) A new QCTOOLS/SNPTEST sample file made with study covariates and phenotypes. 

VIII) ---- 08_run_snpstats_on_bgen.sh ----

	1) Run snptest 
		- on bgen files with rsids
		- using only those samples to be included in the GWAS

IX) ---- 09_make_gwas_bgens.sh ----

	1) Identify SNPs to exclude from GWAS
		- HWE P<1e-12
		- INFO < 0.3
		- MAC < 10
	2) Produce new Study bgen files
		- bgens with rsids
		- only GWAS study samples
		- only SNPs that passed filters above
		
X) ---- 10_order_gwas_samplefile.R ----

	1) Order the sample file made in step 7 
		- to match the sample file made in step 9

XI) ---- 11_run_continuous_GWAS.sh ----

	1) Run all of the AB and trunAB MTs (~10 mins)

XII) ---- 12_run_PA_GWAS.sh ----

	1) Run all of the PA MTs
		- these run longer as they include all of our covariates (~5 hours)

XII) ---- 13_format_snptest_results.sh ----

	1) Concatenate all of the chromosome together
	2) Identify SNPs with P<5e-8 and write to a "tophits" file
	3) Convert concatenated file into BED format
	4) Index the BED file with tabix. 
	
		*Will look to use GWASvcf when I can get all of its dependencies installed* 


## Software

	1) R: r/3.6.1
	2) qctools: qctool/2.0.7
	3) snptest: snptest/2.5.4
	4) plink: plink/2.0.0
	5) GCTA: gcta/1.93.2-beta
	6) tabix: htslib/1.10.2-gcc