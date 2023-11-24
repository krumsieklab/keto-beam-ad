# Script to preprocess serum metabolomics data. 
# Generate files with qced and tranformed data 

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# libraries
library(openxlsx) # for excel reading and writing
library(maplet) # metabolomics toolbox
library(tidyverse) # %>%
library(compositions) # clr
# input files ----
serum_metabolomics <- 'input/BEAM_Serum_Nightingale_raw.xlsx'
anno_file <- 'input/BEAM_Serum_metadata.xlsx'
# output files ----
serum_qcdata <- 'results/tmp_BEAM_Serum_Nightingale_post_qc.xlsx'
serum_pdata <- 'results/tmp_BEAM_Serum_Nightingale_preprocessed.xlsx'

# load data
D <- mt_load_se_xls(file=serum_metabolomics)%>%
  # qc to remove low protein samples ----
  mt_modify_filter_samples(Low_protein=='0')  %>%
  # load sample annotations
  mt_anno_xls(file=anno_file, sheet=1, 
              anno_type="samples", anno_id="sample_id", data_id='sample_ids') %>%
  # write the qced data
  mt_write_se_xls(D, file=serum_qcdata) %>%
  {.}

# preprocessing ----
D1 <- D %>%
  # turn zeros to na and log it
  mt_pre_zero_to_na() %>%
  # Filter out metabolites with more than 40% missingness
  mt_pre_filter_missingness(feat_max = 0.40) %>%
  # Filter out samples with more than 40% missingness
  mt_pre_filter_missingness(samp_max = 0.4) %>% 
  # add missingness percentage as annotation to metabolites
  mt_anno_missingness(anno_type = "features", out_col = "missing")

# percentages - clr transform, imputation and outlier detection
Dpct <- D1[grep('_pct', rowData(D1)$CSV_column_name), ]
X <- Dpct %>% assay() %>% t() %>% as.matrix() %>% clr() %>% t()
assay(Dpct) <- X

Dpct %<>%
  # impute missing values
  mt_pre_impute_knn() %>% 
  # sample outlier removal
  mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>%
  # metabolic outlier detection followed by imputation
  mt_pre_outlier_to_na(use_quant = T, quant_thresh = 0.025) %>%
  # kNN imputation
  mt_pre_impute_knn()

# non pct biomarkers
DO <- D1[-grep('_pct', rowData(D1)$CSV_column_name), ]

# sizes - no transform, imputation, and outlier detection
Dsz <- DO[grep('_size', rowData(DO)$CSV_column_name), ]
Dsz %<>%
  # impute missing values
  mt_pre_impute_knn() %>% 
  # sample outlier removal
  mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>%
  # metabolic outlier detection followed by imputation
  mt_pre_outlier_to_na(use_quant = T, quant_thresh = 0.025) %>%
  # kNN imputation
  mt_pre_impute_knn()

# non size biomarkers
DO <- DO[-grep('_size', rowData(DO)$CSV_column_name), ]

DO %<>%  # Log2 transformation
  mt_pre_trans_log() %>%
  # impute missing values
  mt_pre_impute_knn() %>% 
  # sample outlier removal
  mt_pre_outlier_lof(seq_k = c(5, 10, 20, 30, 40, 50)) %>%
  # metabolic outlier detection followed by imputation
  mt_pre_outlier_to_na(use_quant = T, quant_thresh = 0.025) %>%
  # kNN imputation
  mt_pre_impute_knn()


# # combine the datasets back
samps <- intersect(intersect(DO$Sample_id, Dpct$Sample_id), Dsz$Sample_id)
DO %<>% mt_modify_filter_samples(filter = Sample_id %in%samps)
Dpct %<>% mt_modify_filter_samples(filter = Sample_id %in%samps)
Dsz %<>% mt_modify_filter_samples(filter = Sample_id %in%samps)
# merge two SEs
#all.equal(DO$Sample_id, Dpct$Sample_id); all.equal(Dpct$Sample_id, Dsz$Sample_id);all.equal(DO$Sample_id, Dsz$Sample_id) 
#[1] TRUE
X <- bind_rows(bind_rows(assay(DO) %>% data.frame(), 
               assay(Dpct) %>% data.frame()), 
               assay(Dsz) %>% data.frame())
Y <- bind_rows(bind_rows(rowData(DO)%>% data.frame(), 
               rowData(Dpct) %>% data.frame()),
               rowData(Dsz) %>% data.frame())
Z <- colData(DO) %>% data.frame()
rownames(Z) <- colnames(X)

D2 <- SummarizedExperiment(assay=X, 
                           rowData=Y,
                           colData=Z,
                           metadata=list())
# write the preprocessed data out
mt_write_se_xls(D2, file=serum_pdata)

## finished
print("Done! serum metabolomics preprocessing.") 
print("Generated excel file with supplementary tables in results folder!") 
