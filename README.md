# mbGWAS_regenie
Microbiome GWAS data preparation starting from a phyloseq object. Careful about the setup files, those need to be tweaked
Here is a quick tutorial on how to run it with Nextflow Regenie and some utility functions provided by Michele Filosi:

1. Create the working directory and `cd` there (i.e. "/scratch/work/mbGWAS")

```
mkdir /scratch/work/mbGWAS
cd /scratch/work/mbGWAS
``` 

2. clone Michele Filosi's repository

```
git clone https://gitlab.gm.eurac.edu/mfilosi/regenie_pipeline
```

3. Launch the phyloseq preparation file. This is rather case-specific, so make sure you modify it accordingly:

	- covariates to include
	- their transformation
	- other stuff...
	
```
Rscript 00_data_generation_script_CHRISMB.R
	-b "/scratch/work/mbGWAS"
	-i "path/to/phyloseq.Rds"
	-l "Genus"
	-g "TOPMed"
	-t "IRN"
	-d 10
	-p 0.2
``` 

4. Check sanity of phenotype, covariates and configuration files in the working directory

5. Now you are ready to go:
```
sbatch regenie_pipeline/bin/run_pipeline.sh
```
