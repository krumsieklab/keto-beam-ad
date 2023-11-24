# Customized functions used in other scripts.
# statistical association analysis
association_analysis <- function (
    D,
    # SE
    outcome,
    # outcome to be used for response
    outcome_type,
    # outcome type - numeric/binary/ordinal
    conf_formula = NULL,
    # confounders to be corrected for
    int_w_analyte = NULL,
    # name of covariate that interacts with metabolite
    all_vals = F,
    # return all estimate and p-values from the model
    padj_method = 'BH'){# p-value adjustment method
  # merge data with sample info
  Ds <-
    D %>% maplet:::mti_format_se_samplewise () # NOTE: No explosion of dataset size, no gather() - 6/2/20, JK
  # metabolites in data
  mets <- D %>% assay () %>% rownames ()
  # loop over metabolites
  univ_stats <- lapply (mets, function(x) {
    # turn the outcome variable into factor
    Ds %<>% mutate(!!sym(outcome) := as.factor(as.matrix(!!sym(outcome))))
    # with covariates?
    if (is.null(conf_formula) == F) {
      # with interaction term too
      if (is.null (int_w_analyte) == F) {
        this_formula <- as.formula(glue('{x} ~ {outcome}*{int_w_analyte} + {conf_formula}'))
        # only covariates
      } else {
        this_formula <- as.formula(glue('{x} ~ {outcome} + {conf_formula}'))
      }
      # only metabolite? 
    } else{
      this_formula <- as.formula(glue('{x} ~ {outcome}'))
    }
    # logistic regression
    #this_fit <-glm(this_formula, data = Ds,family = binomial(link = 'logit'))
    # this_fit <- glmer(this_formula, data = Ds, family = binomial, 
    #                  control = glmerControl(optimizer = "bobyqa"), nAGQ = 10)
    this_fit <- lmerTest::lmer(this_formula, data = Ds)
    # results summary
    this_res <- this_fit %>% summary() %>% coefficients() %>% data.frame()
    # format results
    names(this_res) <- c('estimate', 'std_error', 'df', 'statistic', 'p_value')
    if (all_vals) {
      # formating output accordingly 
      tmp <- lapply(2:nrow(this_res), FUN = function(i) this_res[i, ])
      tmp_names <- lapply(rownames(this_res)[2:nrow(this_res)], FUN = function(x)
        paste0(x, sep = '_', c('estimate', 'std_error', 'df', 'statistic', 'p_value')))
      tmp_names[[1]] <- c('estimate', 'std_error', 'df', 'statistic', 'p_value')
      tmp_names[[length(tmp)]] <-paste0(paste0('outcome:', int_w_analyte),sep = '_',
                                        c('estimate', 'std_error', 'df', 'statistic', 'p_value'))
      this_res <-  do.call(cbind, tmp)
      names(this_res) <- unlist(tmp_names)
      # output only for the outcome and metabolite?
    }
    else {
      this_res <- this_res [2, ]
    }
    this_res %<>% mutate(analyte = rownames(this_res), outcome = outcome, covariates = conf_formula)
    # order output columns
    this_res %<>% select(analyte, outcome, covariates, everything())
    
  }) %>% # create data from of results
    do.call(rbind, .) %>% data.frame() %>%
    # bind rowdatas
    bind_cols(D %>% rowData() %>% data.frame() %>% select(name, Biomarker_name, Group, Subgroup))
  # adjust pvalues
  univ_stats %<>% mutate(adj_p = p.adjust(p_value, method = padj_method))
  # order by adjusted p-values
  univ_stats <- univ_stats[order(univ_stats$adj_p), ]
  # adjust the p-values 
  pvals <- names(univ_stats)[grep('_p_value', names(univ_stats))]
  adj_pvals <- sub('_p_value', '_adj_p', pvals)
  for(i in 1:length(pvals)){
    univ_stats %<>% dplyr::mutate(!!as.name(adj_pvals[i]) := p.adjust(!!as.name(pvals[i]), method = padj_method))
  }
  
  # return results
  return(univ_stats)
}

# Function to get fold change from pre- to post- diet
## Input is a summarized experiment
## Output is a dataframe with each metabolite's fold change for each subject
get_diet_fc <-function(D, id_col, grp_col, ci=T){
  
  # Combine subject and diet information with metabolite data
  orig_df <- dplyr::bind_cols((colData(D) %>% data.frame() %>% select(id_col, grp_col)),
                              t(assay(D)))
  
  # For each subject
  diff_df <- lapply(unique(orig_df[[id_col]]), function(i){
    sub_df <- orig_df[which(orig_df[[id_col]]==i), ] # Get only the subject's data
    post_row <- which(sub_df[[grp_col]] == 1) # Post-keto row
    pre_row <- which(sub_df[[grp_col]] == 0) # Pre-keto row
    
    # Subtract metabolte values pre-keto from post-keto
    diffs <- sub_df[post_row, rownames(D)] - sub_df[pre_row, rownames(D)] 
    data.frame(diffs)
    
  }) %>% bind_rows()
  
  if(ci){
    # Get 95% confidence intervaals for the fold changes
    ci_data <- data.frame(t(sapply(diff_df, Rmisc::CI)))
    names(ci_data) <- paste0('outcome_lfc_ci_', names(ci_data))
    ci_data$met <- D %>% rowData() %>% data.frame() %>% pull(name)
    return(ci_data)
  } else{
    return(diff_df)
  }
}
