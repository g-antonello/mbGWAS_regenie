sbatch -c8 --mem=32G -pslow /shared/bioinf/R/bin/Rscript-4.3-BioC3.17 render-gene-level-report.R -s "ADD-SKATO-ACAT" -t 16

# -p is the pvalue option, which by default makes a calculation as shown in line 96-109 of the `.Rmd` file
