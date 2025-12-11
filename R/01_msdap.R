# Mar-2025 ==== author: Brooke Genovese
##============================================================================##
# * MS-DAP WITH BUILT-IN LINEAR REGRESSION ON TOP OF LIMMA *
# * ERB specimens + batch covariate *
##============================================================================##

# --- load libraries --- #
library(dplyr)
library(msdap) 
library(writexl)
library(pcaMethods)

# load data files, output from upstream raw data processor and the exact same fasta file(s) used there
dataset =  msdap::import_dataset_spectronaut(filename = "") # 02_preprocess-ctrl.tsv
dataset = import_fasta(
  dataset,
  files = c("") # ERB fasta and universal contaminents files

# time to step away from R for a sec, and edit this template file in Excel;
# - describe the sample group of each sample in the "group" column
# - further documentation is available in the "instructions" tab within the Excel file

# load sample metadata 
dataset = import_sample_metadata(dataset, filename = "") # metadata excel file

# describe a statistical contrast
dataset = add_contrast(
  dataset,
  colname_condition_variable = "group",
  values_condition1 = c("female"),
  values_condition2 = c("male"),
  colname_additional_variables = c("batch")
)

# main function that runs the entire pipeline
dataset = analysis_quickstart(
  dataset,
  filter_min_detect = 2,            
  filter_min_quant = 0,            
  filter_fraction_detect = 0,    
  filter_fraction_quant = 0,     
  filter_min_peptide_per_prot = 1, 
  filter_by_contrast = TRUE,        
  norm_algorithm = c("vsn", "modebetween_protein"), 
  dea_algorithm = c("deqms"), 
  dea_qvalue_threshold = 0.05,                      
  dea_log2foldchange_threshold = NA,                
  diffdetect_min_peptides_observed = 1,
  diffdetect_min_samples_observed = 2,
  diffdetect_min_fraction_observed = 0,
  output_qc_report = FALSE ,                          
  output_abundance_tables = FALSE,                   
  output_dir = "msdap_results",                     
  output_within_timestamped_subdirectory = FALSE
)


# print a short summary of results at the end
print_dataset_summary(dataset)
dd <- dataset$dd_proteins
dea <- dataset$de_proteins


result_tally <- dea %>%
  summarise(
    qvalue_below_threshold = sum(qvalue < 0.05, na.rm = TRUE),
    abs_foldchange_above_threshold = sum(abs(foldchange.log2) > 1.5, na.rm = TRUE),
    both_criteria_met = sum(qvalue < 0.05 & abs(foldchange.log2) > 1.5, na.rm = TRUE)
  )
print(result_tally)

sig_result_dea <- dea %>%
  mutate(abs_foldchange = abs(foldchange.log2)
  ) %>% 
  filter(qvalue < 0.05 | abs_foldchange > 1.5)

#write_xlsx(list("dea_ctrl_sex" = dea, "dd_ctrl_sex" = dd, "sig_ctrl_log2FC" = sig_result_dea), path = output_file1)


##============================================================================##
# * EXTRACT PCA DATA *
##============================================================================##
# peptide level data matrix, using filtered+normalized peptides across all groups
tibw_noexclude = dataset$peptides %>% 
  select(peptide_id, sample_id, intensity_all_group) %>%
  filter(!is.na(intensity_all_group)) %>%
  pivot_wider(id_cols = peptide_id, names_from = sample_id, values_from = intensity_all_group)
matrix_sample_intensities = msdap:::as_matrix_except_first_column(tibw_noexclude)

# compute PCA
PPCA = pcaMethods::pca(base::scale(t(matrix_sample_intensities), center=TRUE, scale=FALSE), method="ppca", nPcs = 3, seed = 123, maxIterations = 2000)
mat_pca = PPCA@scores 
pca_var = array(PPCA@R2, dimnames = list(colnames(PPCA@scores))) # variation per dimension

# if you only want a table with the PCA coordinates, just print the mat_pca variable and you're done at this point



