t0 <- Sys.time()

library("optparse")

option_list = list(
  make_option(c("-b", "--base.dir"), type="character", default=NA, help = "Base directory"),
  make_option(c("-l", "--level"), type="character", default="ASV", help = "Taxonomic level to aggregate to"),
  make_option(c("-g", "--genotype_dataset"), type="character", default=NA, help = "Genetic Dataset to use for the GWAS, no defaults allowed"),
  make_option(c("-t", "--transform"), type="character", default="irn", help = "Transformation to apply to your dataset"),
  make_option(c("-d", "--detectionFilter"), type="numeric", default=10, help = "Transformation to apply to your dataset"),
  make_option(c("-p", "--prevalenceFilter"), type="numeric", default=0.01, help = "Transformation to apply to your dataset"),
  make_option(c("-H", "--HHID"), type="logical", default=FALSE, action = "store_true", help = "Pick only one person per household?")
)

opt_parser = OptionParser(option_list=option_list)

opts = parse_args(opt_parser)

suppressWarnings(suppressMessages(library(gautils)))

cat("Parameters Chosen:\n")
print(opts)

#+++++++++++++++++++++ IMPORTANT CONFIGURATION PARAMETERS +++++++++++++++++++++

genetic_data <- grep(opts$genotype_dataset, c("WES", "HRC", "TOPMed"), ignore.case = T, value = T) %>% paste0("13K")

proj_title <-   paste("mbGWAS", "on", genetic_data, sep = "__")

dir.create(file.path(opts$base.dir, "regenie_input"), showWarnings = FALSE, recursive = T)
system(paste0("git clone https://gitlab.gm.eurac.edu/mfilosi/regenie_pipeline ", opts$base.dir, "/regenie_pipeline"))
# save opts as object to replicate potential errors
save(opts, file = file.path(opts$base.dir, "regenie_input", "setupObjects.Rda"))
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#----------- LOAD MICROBIOME BASIC DATA AND CHRIS METADATA --------------------
chrismb_phy <- readRDS("/shared/statgen/microbiome/CHRISMB/processed_data/microbiome_tables/phyloseq_built.Rds") 
sample_names(chrismb_phy) <- chrismb_phy@sam_data$AID

chrismb_phy_core <- chrismb_phy %>% 
  core(
    detection = opts$detectionFilter, 
    prevalence = opts$prevalenceFilter
  )

#-------------------------------------------------------------------------------


#-------- Generate the number of taxa needed for 95% variance explained --------
tmp_for_PCoA <- chrismb_phy_core %>% 
  phy_tax_glom(opts$level) %>% 
  phy_transform("compositional") %>% # here the compositionality should be respected, so no IRN 
  ordinate(method = "PCoA", distance = "bray")

saveRDS(min(which(cumsum(tmp_for_PCoA$values$Relative_eig)> 0.95)), paste0(opts$base.dir,"/regenie_input/n_components_0.95.rds"))

rm(tmp_for_PCoA)
#-------------------------------------------------------------------------------
# ----- Create Metadata object to merge into the phyloseq object

library(chrisData)

possible_data_files <- chrisDataFile()

selected_data_files <- c(
  "general_information",
  "household",
  "clinical_traits",
  "interview",
  "labdata",
  "substudies",
  "drugs_summary",
  "selfadmin"
)

data_selected.list <- sapply(selected_data_files, function(x) chrisData(paste0(x, ".txt.gz")), USE.NAMES = T) %>% 
  suppressWarnings() %>% 
  suppressMessages()

data_selected.joined.df <- purrr::reduce(data_selected.list, full_join, by = "AID")
data_selected.joined.df <- data_selected.joined.df[match(chrismb_phy_core@sam_data$AID, data_selected.joined.df$AID),]

# Columns with all NAs must be necessarily removed 
all_NAs <- apply(data_selected.joined.df, 2, function(x) all(is.na(x)))
data_selected.joined.df <- data_selected.joined.df[, !all_NAs]

# Columns of class "chron" are problematic, they must be removed to be able to work with phyloseq and metadata data.frame back and forth
data_selected.joined_NODates.df <- select(data_selected.joined.df, -contains(c("x0bp08", "x0bc03", "x0bc06")))

rm(data_selected.list)

chrismb_phy_core <- phy_add_metadata_variables(chrismb_phy_core, df = data_selected.joined_NODates.df, by = "AID")

#------------ Get basic phyloseq (already filtered) ready for other fancy adventures -------------
phyloseq_ready_for_further_stuff <-  chrismb_phy_core %>% 
  phy_tax_glom(opts$level)
trait_names <- taxa_names(phyloseq_ready_for_further_stuff) %>% 
  make.names() %>% gsub("\\.$", "", .)
taxa_names(phyloseq_ready_for_further_stuff) <- trait_names
#-------------------------------------------------------------------------------


#------------- Generate Traits and Covariates files ---------------------------
genetic_pc_CHRIS13k <- read.table("/shared/statgen/CHRIS13K/CHRIS13K.GT.evecs", header = FALSE, colClasses = c("character", rep("numeric", 20))) %>%
  set_names(c("IID", paste0("PC", 1:20))) %>% 
  select(IID, all_of(paste0("PC", 1:10)))

traits_and_covs_df <- phyloseq_ready_for_further_stuff %>% 
  #phy_alpha_div_into_phy_metadata(5000, measure = c("Observed", "Shannon")) %>% 
  phy_transform(opts$transform) %>% 
  phy_OtuMetaTable() %>% 
  #start generating covariates
  mutate(
    FID = AID,
    IID = AID,
    SMOKING = case_when(grepl("Current", x0sm51) ~ 1, 
                        TRUE ~ 0),
    TEETH = as.integer(x0oh01),
    ABX = case_when(x0dd51 == "No" ~ 0,
                    x0dd51 == "Yes" ~ 1),
    x0_sex = ifelse(x0_sex == "Male", 1, 2)
  ) %>% 
  relocate(FID, IID) %>% 
  inner_join(genetic_pc_CHRIS13k, by = "IID") %>%
  filter(x0_sex == 1)




covariates_final <- select(traits_and_covs_df, 
                           FID, 
                           IID, 
                           #x0_sex,
                           x0_ager,
                           SMOKING, 
                           TEETH, 
                           ABX, 
                           contains("PC",ignore.case = F)
)

traits_final <- select(traits_and_covs_df, FID, IID, all_of(trait_names))

# write input files

write.table(
  covariates_final,
  file.path(opts$base.dir, "regenie_input", "covariates.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  traits_final,
  file.path(opts$base.dir, "regenie_input", "traits.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
write.table(
  phy_BasicStats(
    phyloseq_ready_for_further_stuff %>% filter_sample_data(AID %in% traits_final$IID),
    "compositional"
  ),
  file = file.path(opts$base.dir, "regenie_input", "taxaStats.txt"),
  quote = FALSE,
  row.names = FALSE
)

write.table(phyloseq_ready_for_further_stuff,
            file = file.path(opts$base.dir, "regenie_input", "taxaTable.txt"),
            sep = "\t", 
            quote = FALSE,
	    row.names = FALSE
)

#---------------------- write configuration file ------------------------------

source(file.path(opts$base.dir, "regenie_pipeline", "bin", "base_functions.R"))

dir.create(file.path(opts$base.dir, "conf"), showWarnings = FALSE)

# configuration file and parameters for traits
myconffile <- write.conf(
  x = colnames(traits_final)[3:ncol(traits_final)],
  filename = file.path(opts$base.dir, "conf", "mygwas.conf"),
  binary = FALSE, # Change to TRUE if phenotypes are binary (0/1 coded)
  overwrite = TRUE,
  phenofile = file.path(opts$base.dir,"regenie_input", "traits.txt"),
  project_name = proj_title, # this will be prepended to each of the outputfile
  outpath = file.path(opts$base.dir, "results") # specify where the results directory will go
)

# configuration parameter for covariates

write_covariates2conf(
  x = colnames(covariates_final)[3:ncol(covariates_final)],
  covariate_file = file.path(opts$base.dir,"regenie_input", "covariates.txt"),
  conf_file = myconffile
)

# add genotype data configuration
add_baseconfig2conf(
  dataset = genetic_data,
  conf_file = myconffile,
  conf.path = file.path(opts$base.dir, "regenie_pipeline", "conf")
)


cat("\nR setup part finished")
cat("\n")
print(Sys.time() - t0)
