# Script to find serum metabolites that changed post MMKD and AHAD. 
# Generate files with association results

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))
# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# libraries
library(openxlsx) # for excel reading and writing
library(maplet) # MT 
library(lme4) # mixed effect models
library(tidyverse) # %>%
library(glue) #> formula
source("custom_functions.R") # customized functions

# define empty results list
res <- list()
# input files
data_file <- 'results/tmp_BEAM_Serum_Nightingale_preprocessed.xlsx'
# output files
serum_associations <- 'results/supplementary_table_X_BEAM_Serum_Nightingale_analysis.xlsx'
serum_aha <- 'results/tmp_BEAM_Serum_Nightingale_preprocessed_matched_AHAD.xlsx'
serum_keto <- 'results/tmp_BEAM_Serum_Nightingale_preprocessed_matched_MMKD.xlsx'

# data load
D <- mt_load_se_xls(file=data_file)  %>% 
  # remove samples without age and sex
  mt_modify_filter_samples(filter=!is.na(Age)) %>%
  mt_modify_filter_samples(filter=!is.na(Sex)) %>%
  {.}

# modify annotation columns
D1 <- D %>% 
  # modifying outcome columns
  mt_anno_mutate(anno_type = 'samples', col_name = 'diag', term=
                   case_when(MemDx=='NC' ~0,  MemDx=='NC-SC'~0, 
                             MemDx=='MCI' ~1,  MemDx=='MCI-A'~1,MemDx=='MCI-MDA'~1,  
                             TRUE~NA_real_))%>%
  mt_anno_mutate(anno_type = 'samples', col_name='diag', term=as.factor(as.matrix(diag))) %>%
  mt_anno_mutate(anno_type = 'samples', col_name='Age', term=as.numeric(as.matrix(Age))) %>%
  mt_anno_mutate(anno_type = 'samples', col_name='Sex', term=as.factor(as.matrix(Sex))) %>%
  # create group definitions
  mt_anno_mutate(anno_type = 'samples', col_name = 'keto_group', 
                 term=case_when((first_diet=="keto" & timepoint_NG_serum=='PREA-COG')~0, 
                                (first_diet=="keto" & timepoint_NG_serum=='POSTA-COG')~1,
                                (first_diet=="aha" & timepoint_NG_serum=='PREB-COG')~0, 
                                (first_diet=="aha" & timepoint_NG_serum=='POSTB-COG')~1,
                                TRUE~ NA_real_))%>% 
  mt_anno_mutate(anno_type = 'samples', col_name = 'aha_group', 
                 term=case_when((first_diet=="aha" & timepoint_NG_serum=='PREA-COG')~0, 
                                (first_diet=="aha" & timepoint_NG_serum=='POSTA-COG')~1,
                                (first_diet=="keto" & timepoint_NG_serum=='PREB-COG')~0, 
                                (first_diet=="keto" & timepoint_NG_serum=='POSTB-COG')~1,
                                TRUE~ NA_real_))

# select keto samples
D11 <- mt_modify_filter_samples(D1, filter=!is.na(keto_group)) 
D11 %<>% # select paired samples
  mt_modify_filter_samples(filter=SubjectID %in% (D11$SubjectID[which(duplicated(D11$SubjectID))]))%>% 
  # remove patient with no diagnosis  
  mt_modify_filter_samples(filter=!is.na(diag)) %>%
  # write out the samples analysed
  mt_write_se_xls(file=serum_keto) %>%
  {.}
tmp <- get_diet_fc(D=D11, id_col = "SubjectID", grp_col='keto_group')

# scaling
D11 %<>% mt_pre_trans_scale(center = T)
# keto analysis results
res[['keto_group']] <- association_analysis(D=D11, outcome = 'keto_group', outcome_type = 'twofactor',
                                            int_w_analyte = "diag",
                                            conf_formula = "Age + Sex  + (1|SubjectID)", 
                                            all_vals=T)

res$keto_group <- left_join(res$keto_group, tmp, by=c('name'='met'))

# select aha samples 
D11 <- mt_modify_filter_samples(D1, filter=!is.na(aha_group))
D11 %<>% # select paired samples
  mt_modify_filter_samples(filter=SubjectID %in% (D11$SubjectID[which(duplicated(D11$SubjectID))])) %>%
  # remove patient with no diagnosis
  mt_modify_filter_samples(filter=!is.na(diag)) %>%
  # write out the samples analysed
  mt_write_se_xls(file=serum_aha) %>%
  {.}

tmp <- get_diet_fc(D=D11, id_col = "SubjectID", grp_col='aha_group')

# scaling
D11 %<>% mt_pre_trans_scale(center = T)
# aha results
res[['aha_group']] <- association_analysis(D=D11, outcome = 'aha_group', outcome_type = 'twofactor',
                                           int_w_analyte = "diag",
                                           conf_formula = "Age + Sex  + (1|SubjectID)", 
                                           all_vals=T)
res$aha_group <- left_join(res$aha_group, tmp, by=c('name'='met'))
# write out the results
sheet_names <- list() 
sheet_names[['keto_group']] <- 'keto serum';  sheet_names[['aha_group']] <- 'aha serum';

wb <- openxlsx::createWorkbook()
# loop over outcomes
for(out in c('keto_group', 'aha_group')){
  # results of this outcome
  model <- 'metabolite ~ prepost + diagnosis + age + sex  + prepost:diagnosis + (1|subjectID)'
  this_res <- res [[out]] %>% .[order(.$adj_p), ] %>% mutate(model= model) %>%
    select(-analyte, -outcome, -covariates) %>% 
    select(Biomarker_name, Group, Subgroup, name, model, estimate, std_error, df, statistic, p_value, adj_p, outcome_lfc_ci_upper, outcome_lfc_ci_mean, outcome_lfc_ci_lower, diag1_estimate, diag1_std_error, diag1_df, diag1_statistic, diag1_p_value, diag1_adj_p,
           outcome.diag_estimate, outcome.diag_std_error, outcome.diag_df, outcome.diag_statistic, outcome.diag_p_value, outcome.diag_adj_p,
           Age_estimate, Age_std_error, Age_df, Age_statistic, Age_p_value, Age_adj_p, Sex2_estimate, Sex2_std_error, Sex2_df, Sex2_statistic, Sex2_p_value, Sex2_adj_p)
  names(this_res) <- c("Biomarker_name",	"Group",	"Subgroup",	"name",	"model",	"prepost_estimate", "prepost_std_error",	"prepost_df",	"prepost_statistic",	"prepost_p_value",	"prepost_adj_p",  "prepost_lfc_ci_upper", "prepost_lfc_ci_mean", "prepost_lfc_ci_lower", "diagnosis_estimate",	"diagnosis_std_error",	"diagnosis_df",	"diagnosis_statistic",	"diagnosis_p_value",	"diagnosis_adj_p",	"prepost:diagnosis_estimate",	"prepost:diagnosis_std_error",	"prepost:diagnosis_df",	"prepost:diagnosis_statistic",	"prepost:diagnosis_p_value",	"prepost:diagnosis_adj_p", "age_estimate",	"age_std_error",	"age_df", "age_statistic",	"age_p_value",	"age_adj_p",	"sex_estimate",	"sex_std_error",	"sex_df",	"sex_statistic",	"sex_p_value",	"sex_adj_p"	)
  out <- sheet_names[[out]]
  # create worksheet
  openxlsx::addWorksheet(wb,sprintf('%s', out))
  # write data
  openxlsx::writeData(wb, sprintf('%s', out), this_res,rowNames = F, colNames = T)
  # create and add a style to the column headers
  headerStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center', textDecoration = 'bold')
  # style for body
  bodyStyle <- createStyle(fontName = 'Arial', fontSize = 12, halign = 'center', valign = 'center')
  # apply style
  addStyle(wb, sheet = sprintf('%s', out), bodyStyle, rows = 1:(nrow(this_res)+1), cols = 1:ncol(this_res), gridExpand = TRUE)
  addStyle(wb, sheet = sprintf('%s', out), headerStyle, rows = 1, cols = 1:ncol(this_res), gridExpand = TRUE)
}

# write workbook
openxlsx::saveWorkbook (wb, file=serum_associations, overwrite=TRUE)

## finished
print("Done! serum metabolomics analysis with diet-related groups completed.") 
print("Generated excel file with supplementary tables in results folder!") 
