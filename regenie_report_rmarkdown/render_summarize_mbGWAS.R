#!/usr/bin/env Rscript

library(optparse)

# to play with parameters, try from here: https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/

option_list = list(
  make_option(c("-f", "--fast"), type="logical", default=FALSE, action = "store_true", help = "Whether it should make per-trait report. Careful: it takes a long time"),
  make_option(c("-y", "--plotly"), type="logical", default=FALSE, action = "store_true", help = "Whether it should make interactive plots with plotly"),
  make_option(c("-t", "--threads"), type="integer", default=8, help = "Number of threads given to the process"),
  make_option(c("-v", "--verbose"), type="logical",  default=FALSE, action = "store_true", help = "Makes the script a little more chatty")
  )

opt_parser <- OptionParser(option_list=option_list)

opts_args <- parse_args(opt_parser)

library(rmarkdown)
suppressWarnings(suppressMessages(library(tidyverse)))

t0 <- Sys.time()

# base.path <- file.path(dirname(this.path::this.path()), "../")
base.path <- file.path(here::here(), "../")


config0 <- readLines(file.path(base.path, "/conf/mygwas.conf"))

project <- grep("project", config0, value = TRUE) %>%
	strsplit(" = ", fixed = T) %>% .[[1]] %>% .[2] %>%
	str_replace(pattern = "\\'", replacement = "") %>%
	str_replace(pattern = "\\'", replacement = "")

outfile <- file.path(base.path, paste0("Full-Report_", project, ".html"))

genetic_dataset <- list("hrc" = "HRC",
     "wes" = "WES",
     "topmed" = "TOPMed")[[sapply(c("hrc","wes", "topmed"), function(d) str_extract(tolower(base.path), d)) %>% 
  .[!is.na(.)]]]

cat("\nGenetic dataset used:", genetic_dataset, sep = "\t")

render(list.files(base.path, recursive=T, pattern="summarise_mbGWAS", full.names = TRUE),
       output_file = outfile)# p value threshold under which to filter out the results. this is useful to minimize the amount of dots.


t1 <- Sys.time()

cat("Time taken for completion:\n")

print(t1-t0)
