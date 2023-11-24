# Script to preprocess CSF metabolomics data. 
# Generate files with qced and tranformed data 

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# libraries
library(openxlsx) # for excel reading and writing
library(maplet) # metabolomics toolbox
library(tidyverse) # %>%

# input files ----
csf_metabo <- 'input/BEAM_CSF_Nightingale_raw.xlsx'
anno_file <- 'input/BEAM_CSF_metadata.xlsx'

# output files ----
csf_qcdata <- 'results/tmp_BEAM_CSF_Nightingale_post_qc.xlsx'
csf_pdata <- 'results/tmp_BEAM_CSF_Nightingale_preprocessed.xlsx'

# preprocessing ----
D <- mt_load_nightingale (file=csf_metabo, format_type = "multiple_sheets_v2") %>%
  mt_anno_xls(file=anno_file, sheet=1, anno_type = 'samples', 
                     anno_id_col = "sample_ID", data_id_col = "Sample_id")  %>%
  # write out data
  mt_write_se_xls(file = csf_qcdata)

  # turn zeros to na 
D %<>% #mt_pre_zero_to_na() %>%
  # Filter out metabolites with more than 40% missingness
  mt_pre_filter_missingness(feat_max = 0.4) %>%
  # Filter out samples with more than 40% missingness
  mt_pre_filter_missingness(samp_max = 0.4) %>%
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
  # Log2 transformation
  mt_pre_trans_log() %>% 
  # impute missing values
  mt_pre_impute_knn() %>%
  # sample outlier removal
  mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>%
  # metabolic outlier detection followed by imputation
  mt_pre_outlier_to_na(use_quant = T, quant_thresh = 0.025) %>%
  # kNN imputation
  mt_pre_impute_knn() %>%
  # write out data
  mt_write_se_xls(file = csf_pdata)


## finished
print("Done! CSF metabolomics preprocessing.") 
print("Generated excel file with supplementary tables in results folder!") 