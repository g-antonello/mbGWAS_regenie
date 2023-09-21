#!/usr/bin/env Rscript


#---------------------- SETUPS STEPS ------------------------------
library(optparse)

# to play with parameters, try from here: https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/

option_list = list(
  make_option(c("-s", "--statistical_test"), type="character", default="ADD-SKATO-ACAT", help = "type of test to extract. the default is 'ADD-SKATO-ACAT'"),
  make_option(c("-t", "--threads"), type="integer", default=8, help = "Number of threads given to the process, default is 8"),
  make_option(c("-p", "--p_value_threshold"), type="integer", default= "autoFind", help = "-log10(Pvalue) threshold. Values above it will be kept"),
  make_option(c("-l", "--log10"), type="logical",  default=TRUE, action = "store_true", help = "Are P-values stored as -log10 or as crude values?. The default is: yes, they are stored as -log10"),
  make_option(c("-v", "--verbose"), type="logical",  default=FALSE, action = "store_true", help = "Makes the script a little more chatty"),
  
  
)

opt_parser <- OptionParser(option_list=option_list)

opts_args <- parse_args(opt_parser)

library(rmarkdown)
library(data.table)
suppressWarnings(suppressMessages(library(tidyverse)))

t0 <- Sys.time()

# Set the working directory and extract GWAS runs info from mother GWAS 

base.path <- file.path(here::here(), "../")

cat("Base path", base.path, sep = ": ")
config0 <- readLines(file.path(base.path, "conf", "mygwas_gene_level.conf"))

project <- grep("project", config0, value = TRUE) %>%
  strsplit(" = ", fixed = T) %>% .[[1]] %>% .[2] %>%
  str_replace(pattern = "\\'", replacement = "") %>%
  str_replace(pattern = "\\'", replacement = "")


taxonomy_info <- taxaTab <- fread(file.path(base.path, ".." , "regenie_input", "taxaTable.txt")) %>% 
  as.data.frame()

taxLvl <- colnames(taxonomy_info)[ncol(taxonomy_info)]

outfile <- file.path(base.path, paste0("Gene-Level-Report_", project, "_", taxLvl, ".html"))

genetic_dataset <- list("hrc" = "HRC",
                        "wes" = "WES",
                        "topmed" = "TOPMed")[[sapply(c("hrc","wes", "topmed"), function(d) str_extract(tolower(base.path), d)) %>% 
                                                .[!is.na(.)]]]


basic_setup_gene_level <- grep("regenie_base_conf", config0, value = T) %>% 
  gsub('includeConfig \\"', "", .) %>% 
  gsub('\\"', '', .) %>% 
  readLines()
#-------------------------------------------------------------------------------


cat("\nGenetic dataset used:", genetic_dataset, sep = "\t")

print(list.files(here::here(), recursive=T, pattern="gene-level-report", full.names = TRUE))

render(file.path(here::here(), "gene-level-report.Rmd"),
       output_file = outfile)# p value threshold under which to filter out the results. this is useful to minimize the amount of dots.


t1 <- Sys.time()

cat("Time taken for completion:\n")

print(t1-t0)
