# Script to generate figure panels from the paper. Panels of Figures 2, 3a, 4 (Figures 1, 3b, and 4 (partially) were manually generated) 
# Generate pdfs/html with plots in the results folder

# clear workspace
rm(list = setdiff(ls(),c("codes.makepath","data.makepath","results.makepath")))
# set working directory to location of source code
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# libraries
library(readxl)
library(maplet)
library(dplyr)
library(MetBrewer)
library(SummarizedExperiment)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(magrittr)
library(Hmisc)
library(Rmisc)
library(igraph)
library(visNetwork)
library(stringr)
source("custom_functions.R") # customized functions
# get colors
palette <- c(MetBrewer::met.brewer("Archambault"), "#0c7156")

# Read in more annotations for significant diet-associated serum metabolites
sig_met_annos <- read.csv("input/tmp_annotations_sig_serum_diet.csv")


##### Read in MMKD results ####
# Load in all the preprocessed (logged) values for serum
all_serum_D <- maplet::mt_load_se_xls(file="results/tmp_BEAM_Serum_Nightingale_preprocessed_matched_MMKD.xlsx")


# Load in all the preprocessed (logged) values for CSF
all_csf_D <- maplet::mt_load_se_xls(file="results/tmp_BEAM_CSF_Nightingale_preprocessed_matched_MMKD.xlsx")


# Read in serum data
interaction_res_serum <- readxl::read_xlsx("results/supplementary_table_X_BEAM_Serum_Nightingale_analysis.xlsx", sheet="keto serum")[,c(4, 10, 11, 6, 19, 15, 25, 21)] %>%
  data.frame()

# Rename columns
names(interaction_res_serum)<- c("Metabolite", "Diet_raw_p_val_serum","Diet_p_val_serum_adj","Diet_direction_serum", "Cog_p_val_serum_adj","Cog_direction_serum","Interaction_p_val_serum_adj","Interaction_direction_serum")


# Read in CSF data
interaction_res_csf <- readxl::read_xlsx("results/supplementary_table_X_BEAM_CSF_Nightingale_analysis.xlsx", sheet="keto csf")[,c(4, 10, 11, 6, 19, 15, 25, 21)] %>%
  data.frame()

# Rename columns
names(interaction_res_csf)<- c("Metabolite","Diet_raw_p_val_csf", "Diet_p_val_csf_adj", "Diet_direction_csf","Cog_p_val_csf","Cog_direction_csf","Interaction_p_val_csf","Interaction_direction_csf")


# Combine CSF and Serum results
interaction_res <- merge(interaction_res_serum, interaction_res_csf, by="Metabolite", all=T)
####

#### Load in original data (no logging or scaling!!) ####
orig_data <- maplet::mt_load_se_xls(file="results/tmp_BEAM_Serum_Nightingale_post_qc.xlsx")

# Get annotations for metabolites from original data
met_annos <- orig_data %>%
  rowData() %>%
  data.frame() %>%
  dplyr::select(Excel_column_name, Biomarker_name, Group, Subgroup)

# Combine statistical results with annotations
ps_w_annos <- merge(interaction_res, 
                    met_annos, 
                    by.x= "Metabolite", 
                    by.y = "Excel_column_name", 
                    all=T)


#### Preprocessed data loading ####

# Set ID to ID and keto 
colData(all_serum_D)$ID <- paste0(colData(all_serum_D)$SubjectID, colData(all_serum_D)$keto_group)
colnames(all_serum_D)<-colData(all_serum_D)$ID

# Set ID to ID and keto 
colData(all_csf_D)$ID <- paste0(colData(all_csf_D)$SubjectID, colData(all_csf_D)$keto_group)
colnames(all_csf_D)<-colData(all_csf_D)$ID

# Starting with diet, get significant metabolites (adj. p < 0.2)
sig_serum_diet <- ps_w_annos %>%
  dplyr::filter(Diet_p_val_serum_adj<=0.2) %>% 
  dplyr::mutate(anno = ifelse(grepl("olester", Group), "Cholesterol concentrations", Group)) %>%
  dplyr::mutate(anno = ifelse(grepl("ipoprot", anno), "Lipoprotein related metrics", anno)) %>%
  dplyr::mutate(anno = ifelse(anno == "Amino acids", sub(" ","<br>",anno), anno))


# For cognition, get significant metabolites (adj. p < 0.2)
# Add group annotations

sig_serum_cog <-  ps_w_annos %>%
  dplyr::filter(Cog_p_val_serum_adj<=0.2) %>% 
  dplyr::mutate(anno = ifelse(grepl("olester", Group), "Cholesterol concentrations", Group)) %>%
  dplyr::mutate(anno = ifelse(grepl("ipoprot", anno), "Lipoprotein related metrics", anno)) %>%
  dplyr::mutate(anno = ifelse(anno == "Amino acids", sub(" ","<br>",anno), anno))

# Subset serum summarized experiment for significant diet-associated metabolites 
sig_serum_D <- all_serum_D[rowData(all_serum_D)$Excel_column_name %in% sig_serum_diet$Metabolite,]
serum_sig_mets <- rownames(assay(sig_serum_D))

# Subset CSF summarized experiment for significant diet-associated metabolites 

## The significance cutoff for CSF is dependent on that of serum
# As there are fewer CSF metabolites, the higher p-val cutoff is more lenient
# and therefore we will make the significance cut at the raw p-value that
# aligns with a significant hit in serum

sig_csf_diet <- ps_w_annos %>%
  dplyr::filter(Diet_p_val_csf_adj<=0.05)

# Add annotations to metabolites measured in CSF
sig_csf_diet$annos <- c("BCAA Metabolism", "Lipid Oxidation", "BCAA Metabolism", "BCAA Metabolism", "BCAA Metabolism", 
                        "Energy Metabolism", "Energy Metabolism", "Energy Metabolism", "Amino Acids","Energy Metabolism",
                        "Amino Acids","Energy Metabolism", "Energy Metabolism", "Other", "Amino Acids","Amino Acids","Other", "Amino Acids")



sig_csf_D <- all_csf_D[rowData(all_csf_D)$Excel_column_name %in% sig_csf_diet$Metabolite,]
sig_mets_csf <- rownames(assay(sig_csf_D))

#### Figure 2 ####


##### Figure 2a: Forest Plots #####
serum_ci_data <-readxl::read_xlsx("results/supplementary_table_X_BEAM_Serum_Nightingale_analysis.xlsx", sheet="keto serum")[,c(4, 12, 13, 14)] %>%
  data.frame()
names(serum_ci_data) <- c('met', 'upper', 'mean', 'lower')
# Combine confidence intervals with annotations
ci_data <- serum_ci_data %>% 
  merge(data.frame(rowData(all_serum_D)), by.x = "met", by.y = "name") %>% 
  data.frame() %>% 
  merge(sig_met_annos, by = "Biomarker_name")

# Change annotation names
ci_data$Biomarker_name[which(ci_data$met =="Total_BCAA")]<- 
  "Total concentration of branched-chain amino acids"

ci_data$Biomarker_name[which(ci_data$met =="HDL_size")]<- 
  "Average diameter for HDL particles *"

# Figure 2a forest plot
fig_2a_plot <- ggplot(ci_data, aes(x=mean, y=Biomarker_name))+
  
  geom_vline(xintercept=0, color="black")+
  
  geom_linerange(aes(xmin=lower, xmax=upper),linewidth=2,color = palette[5])+
  
  geom_point(aes(x=mean), shape=21, size=5, fill = palette[5], color="black") +
  
  facet_grid(as.formula(sprintf("anno_pie~.")), scales = "free_y", space = "free_y") +
  
  xlab('Fold Change (95% Confidence Interval)')+
  
  ylab('Metabolite') +
  
  theme(strip.background =element_rect(fill=NA),
        strip.text = element_text(colour = 'black', face = "bold", size=12),
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text=element_text(size=10),
        panel.grid.major.y = element_line(color ="gray"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill=NA, color ="black"))+
  
  xlim(-0.6,4.3)

pdf("results/Figure_2a.pdf", width=10)
plot(fig_2a_plot)
dev.off()

##### Figure 2b: Serum significant hits #####

# Keto Diet top sig hits

##### BCAAs ######
bcaas <-  c("Val", "Leu", "Ile")

# Full names for plotting
full_names <- rowData(all_serum_D)$Biomarker_name[match(bcaas, rowData(all_serum_D)$Excel_column_name)]

names(full_names) <- bcaas

# Subset data to only bcaas
bcaas_df <- data.frame(t(assay(all_serum_D)[match(bcaas, rowData(all_serum_D)$Excel_column_name),]),
                       colData(all_serum_D)$keto_group,
                       colData(all_serum_D)$SubjectID)
names(bcaas_df) <- c(full_names, "Diet", "subj")

bcaas_df %<>% reshape2::melt(measure.vars = full_names)

# Significance annotations
annotation_df <- data.frame(
  variable = full_names,
  Diet = c(rep(c(0,1), length(bcaas))),
  start = c(rep(0.8,length(bcaas))),
  end = c(rep(2.3, length(bcaas))),
  y = .05 + unlist(lapply(as.vector(full_names), function(i) max(bcaas_df %>% filter(variable == i) %>% .$value))),
  pval = sig_serum_diet$Diet_p_val_serum_adj[match(bcaas, sig_serum_diet$Metabolite)] %>%  round(3)
) %>%
  dplyr::mutate(label = ifelse(pval < 0.05, "*", "n.s"))

# Figure 2b BCAA plot
fig_2b_plot<- ggplot(bcaas_df, aes(x=as.factor(Diet), y = value, fill = as.factor(Diet))) +
  
  geom_boxplot(alpha = 0.8) +
  
  geom_point(position = position_dodge(width = .75))+
  
  geom_path(aes(group = subj))+
  
  scale_fill_manual(name = "Diet", labels = c("Pre-Keto Diet","Post-Keto Diet"),values=c(palette[3], palette[1], palette[4])) +
  
  theme_bw() +
  
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        strip.text = element_text(size=24),
        legend.title = element_blank(),
        strip.background = element_rect(fill="seashell")) +
  
  ggtitle("Top metabolite changes in the Ketogenic Diet") +
  
  ylab("Normalized Concentration") +
  
  xlab("")+
  
  geom_signif(
    data = annotation_df,
    aes(xmin = start, xmax = end, annotations = label, y_position = y),
    textsize = 10, vjust = 0.5,
    manual = TRUE
  ) +
  
  facet_wrap(~factor(variable, levels = full_names), scales="free", nrow=1)+
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

pdf("results/Figure_2b.pdf")
plot(fig_2b_plot)
dev.off()

##### Figure 2d: HDL plot #####

# Cholesterol names
conc_names <- c("_FC","_TG","_PL","_CE", "_P")

# Function to find cholesterol names within metabolite names
fun <- function(x, y) {
  grepl(x, y)
}

# Get indices of lipoprotein lipid concentrations (don't want percentages)
conc_indices <- which(unlist(lapply(rownames(orig_data),
                                    function(i) sum(mapply(fun,conc_names, i))>0 &!(grepl("_pct",i)) &!(grepl("Total",i)))))

# Combine raw data of lipoprotein lipid concentrations (unlogged) with diet data
conc_df <- data.frame(t(assay(orig_data)[conc_indices,]),
                      SubjectID = colData(orig_data)$SubjectID,
                      sample_id = colData(orig_data)$sample_id) %>%
  merge(data.frame(colData(all_serum_D), by=c("SubjectID", "sample_id"), all.y=T)) %>% 
  select(all_of(c(rownames(orig_data[conc_indices,]), "SubjectID", "keto_group"))) %>% 
  distinct() %>% 
  dplyr::rename(Diet=keto_group) %>% 
  melt(id.var=c("SubjectID", "Diet")) 

# Get lipoprotein (HDL, LDL, VLDL) and fat (Cholesterol, FC, CE, etc) classes
classes <- bind_rows(lapply(conc_df$variable, function(i){
  split <- stringr::str_split(i,"_")[[1]]
  data.frame(lipo_class = paste(split[1:(length(split)-1)], collapse="_"),
             fat_class = split[length(split)])
}))


# Combine concentrations with their classes
conc_df <- cbind(conc_df, classes)

# Figure 2d - Cholesterol Esters in HDL particles
fig_2di_plot <- ggplot(conc_df %>% 
         dplyr::filter(fat_class == "CE") %>% 
         dplyr::filter(grepl("HDL", lipo_class)),
       aes(x = factor(lipo_class, levels = c("HDL", 'S_HDL','M_HDL','L_HDL','XL_HDL')), 
           y = value, 
           fill = as.factor(Diet)))+
  geom_boxplot() +
  theme_bw() +
  xlab("")+
  ylab("Cholesterol Ester Concentration (mmol/l)")+
  scale_fill_manual(name = "Diet", 
                    labels = c("Pre-Keto Diet","Post-Keto Diet"),
                    values=c(palette[3], palette[1], palette[4])) +
  scale_x_discrete(labels = c("All\nHDL", "Small\nHDL","Medium\nHDL",'Large\nHDL','X-Large\nHDL'))+
  theme(axis.text.x=element_text(size=25),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.background = element_rect(fill="seashell"))+
  geom_signif(annotations = c("*","*","*"),
              y_position = c(.5, 0.1, -0.4), 
              xmin=c(0.75,3.75, 4.75), 
              xmax=c(1.25,4.25,5.25), 
              textsize = 10,
              vjust=0.5)+
  scale_y_log10() 

# Figure 2d - Free Cholesterol in HDL particles
fig_2dii_plot <- ggplot(conc_df %>% 
         filter(fat_class == "FC") %>% 
         filter(grepl("HDL", lipo_class)),
       aes(x = factor(lipo_class, levels = c("HDL", 'S_HDL','M_HDL','L_HDL','XL_HDL')), 
           y = value, 
           fill = as.factor(Diet)))+
  geom_boxplot() +
  theme_bw() +
  xlab("")+
  ylab("Free Cholesterol Concentration (mmol/l)")+
  scale_fill_manual(name = "Diet", 
                    labels = c("Pre-Keto Diet","Post-Keto Diet"),
                    values=c(palette[3], palette[1], palette[4])) +
  scale_x_discrete(labels = c("All\nHDL", "Small\nHDL","Medium\nHDL",'Large\nHDL','X-Large\nHDL'))+
  theme(axis.text.x=element_text(size=25),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.background = element_rect(fill="seashell"))+
  
  
  scale_y_log10() +  geom_signif(annotations = c("*","*"),
                                 y_position = c(-0.04, -0.06), 
                                 xmin=c(0.75,3.75), 
                                 xmax=c(1.25,4.25), 
                                 textsize = 10,
                                 vjust=0.5)

# Figure 2d - Total concentration of HDL particles
fig_2diii_plot <- ggplot(conc_df %>% 
         filter(fat_class == "P") %>% 
         filter(grepl("HDL", lipo_class)),
       aes(x = factor(lipo_class, levels = c("HDL", 'S_HDL','M_HDL','L_HDL','XL_HDL')), 
           y = value, 
           fill = as.factor(Diet)))+
  geom_boxplot() +
  theme_bw() +
  xlab("")+
  ylab("Concentration of Particles (mmol/l)")+
  scale_fill_manual(name = "Diet", 
                    labels = c("Pre-Keto Diet","Post-Keto Diet"),
                    values=c(palette[3], palette[1], palette[4])) +
  scale_x_discrete(labels = c("All\nHDL", "Small\nHDL","Medium\nHDL",'Large\nHDL','X-Large\nHDL'))+
  theme(axis.text.x=element_text(size=25),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        legend.position = "bottom",
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        strip.background = element_rect(fill="seashell"))+
  geom_signif(annotations = c("*","*"),
              y_position = c(-2, -3), 
              xmin=c(3.75,4.75), 
              xmax=c(4.25,5.25), 
              textsize = 10,
              vjust=0.5)+
  scale_y_log10()

pdf("results/Figure_2di.pdf")
print(fig_2di_plot)
dev.off()

pdf("results/Figure_2dii.pdf")
print(fig_2dii_plot)
dev.off()

pdf("results/Figure_2diii.pdf")
plot(fig_2diii_plot)
dev.off()

#### Figure 3: CSF significant hits ####

###### CSF Forest Plots #####

# Get 95% confidence intervals for fold changes of csf metabolites
csf_ci_data <-readxl::read_xlsx("results/supplementary_table_X_BEAM_CSF_Nightingale_analysis.xlsx", sheet="keto csf")[,c(4, 12, 13, 14)] %>%
  data.frame()
names(csf_ci_data) <- c('met', 'upper', 'mean', 'lower')

# Combine confidence intervals with annotations
ci_data <- csf_ci_data%>% 
  merge(data.frame(rowData(all_csf_D)), by.x = "met", by.y = "CSV_column_name") %>% 
  data.frame() %>% 
  merge(ps_w_annos, by.x="met",by.y="Metabolite") %>% 
  dplyr::rename("Biomarker_name" = "Biomarker_name.x") %>% 
  dplyr::rename("Group" = "Group.x")

# Add metabolite annotations for CSF metabolites from HMDB and literature
# HMDB
ci_data$Group[ci_data$Biomarker_name %in% c("Creatinine", "Creatine")] <- 'ATP Cycle'
# literature
ci_data$Group[ci_data$Biomarker_name %in% c("2-Hydroxybutyrate", '2-Ketoisovalerate','2-Hydroxyisovalerate', '3-Hydroxyisobutyrate', '3-Hydroxyisovalerate')] <- "BCAA Metabolism"
ci_data$Group[ci_data$Group =="Miscellaneous"] <- 'Other'

# Merge serum CI data with annotations, map to CSF metabolite annotations
serum_ci_data_csf <- serum_ci_data %>% 
  merge(data.frame(rowData(all_serum_D)), by.x = "met", by.y = "name") %>% 
  data.frame() %>% 
  dplyr::filter(met %in% ci_data$met)

# Clean up names
serum_ci_data_csf$Group[serum_ci_data_csf$Biomarker_name %in% c("Creatinine", "Creatine")] <- 'ATP Cycle'
serum_ci_data_csf$Group[serum_ci_data_csf$Biomarker_name %in% c("2-Hydroxybutyrate", '2-Ketoisovalerate','2-Hydroxyisovalerate', '3-Hydroxyisobutyrate', '3-Hydroxyisovalerate')]<- "BCAA Metabolism"
serum_ci_data_csf$Group[serum_ci_data_csf$Group =="Miscellaneous"] <- 'Other'

# Plot figure 3a
fig_3a_plot <- ggplot(ci_data, aes(x=mean, y=Biomarker_name))+
  geom_vline(xintercept=0, color="black")+
  geom_linerange(aes(xmin=lower, xmax=upper),linewidth=2,color = palette[1])+
  geom_point(aes(x=mean), shape=21, size=5, fill = palette[1], color="black") +
  geom_point(data = serum_ci_data_csf, shape=21, size=5, fill = palette[5], color="black", alpha = 0.5) +
  geom_linerange(data = serum_ci_data_csf,aes(xmin=lower, xmax=upper), linewidth = 2, alpha = 0.7, color = palette[5])+
  
  facet_grid(as.formula(sprintf("Group~.")), scales = "free_y", space = "free_y") +
  
  xlab("Fold Change with 95% confidence interval") +
  
  ylab("Metabolite")+
  
  theme(strip.background =element_rect(fill=NA),
        strip.text = element_text(colour = 'black', face = "bold", size=15),
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.text=element_text(size=12),
        panel.grid.major.y = element_line(color ="gray"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill=NA, color ="black"))+
  xlim(-2, 4.5)

pdf("results/Figure_3a.pdf", width = 10)
print(fig_3a_plot)
dev.off()

##### Figure 4 #####

# Get overlapping samples between serum and CSF
sample_overlap <- intersect(colData(all_serum_D)$SubjectID, colData(all_csf_D)$SubjectID)

# Get serum deltas of significant metabolites, edit rownames to establish source
serum_diffs <- get_diet_fc(D=all_serum_D, id_col = "SubjectID", grp_col='keto_group', ci=F)[sample_overlap,rownames(sig_serum_D)]
colnames(serum_diffs) <- paste0(colnames(serum_diffs), "_serum")

# Get csf deltas of significant metabolites, edit rownames to establish source
csf_diffs <-  get_diet_fc(D=all_csf_D, id_col = "SubjectID", grp_col='keto_group', ci=F)[sample_overlap,rownames(sig_csf_D)]
colnames(csf_diffs) <- paste0(colnames(csf_diffs), "_csf")


# combine deltas
keto_deltas <- cbind(serum_diffs, csf_diffs)

#### Regular correlation plots of keto metabolites

# Calculate the correlation between the deltas
cor_info <- Hmisc::rcorr(matrix(unlist(keto_deltas), nrow = nrow(keto_deltas)))
cor_mat <- cor_info$r
rownames(cor_mat)<-colnames(cor_mat) <- colnames(keto_deltas)

# P-value adjustment of correlations
p_adj_mat <- matrix(1, nrow = nrow(cor_info$P), ncol = ncol(cor_info$P))

cor_ps <- cor_info$P[upper.tri(cor_info$P)] %>% 
  as.vector() %>% 
  p.adjust(method="fdr") 

p_adj_mat[upper.tri(p_adj_mat)] <- cor_ps
rownames(p_adj_mat)<-colnames(p_adj_mat) <- colnames(cor_mat)

# Make the CSF-Serum bipartite graph
p_adj_mat_bipartite <- matrix(1, nrow=nrow(p_adj_mat), ncol = ncol(p_adj_mat))
p_adj_mat_bipartite[1:nrow(sig_serum_D), (1+nrow(sig_serum_D)):ncol(p_adj_mat)] <- p_adj_mat[1:nrow(sig_serum_D), (1+nrow(sig_serum_D)):ncol(p_adj_mat)]
bipartite_adj_mat <- 1*(p_adj_mat_bipartite<0.05)
rownames(bipartite_adj_mat)<-colnames(bipartite_adj_mat) <- colnames(keto_deltas)

# Keep nodes with significant correlations within the bipartite graph
to_keep <- c(which(colSums(bipartite_adj_mat)>0),which(rowSums(bipartite_adj_mat)>0))

# Convert adjacency matrix to igraph object
keto_regcor_G<-igraph::graph_from_adjacency_matrix(bipartite_adj_mat[to_keep, to_keep],
                                                   add.rownames = T, mode="undirected")

# Set node features to source
igraph::V(keto_regcor_G)$NodeType <- c(rep("csf", sum(colSums(bipartite_adj_mat)>0)), 
                                       rep("serum",sum(rowSums(bipartite_adj_mat)>0)))

# Set x coordinate for serum vs. csf
igraph::V(keto_regcor_G)$x <- c(rep(2,  sum(colSums(bipartite_adj_mat)>0)),
                                rep(1, sum(rowSums(bipartite_adj_mat)>0)))

# Set y coordinate
igraph::V(keto_regcor_G)$y <- c( 2*1: sum(colSums(bipartite_adj_mat)>0),
                                 2*1:sum(rowSums(bipartite_adj_mat)>0))

# Set color for serum vs. csf
igraph::V(keto_regcor_G)$color <- c(rep(palette[3],  sum(colSums(bipartite_adj_mat)>0)), 
                                    rep(palette[5],sum(rowSums(bipartite_adj_mat)>0)))

igraph::V(keto_regcor_G)$size=rep(12, length(to_keep))

# Set labels of nodes manually
igraph::V(keto_regcor_G)$label <- c("Glutamine","Phenylalanine", "Tyrosine", "Isoleucine", "Fructose", "Glucose", 'Myo-Inositol', 
                            "3-Hydroxyisovalerate","Formate","Urea","2-Hydroxybutyrate", "Alanine", "Glutamine", "Glycine", 
                            "Total BCAAs","Valine", "3-Hydroxybutyrate","Acetoacetate","Acetone", 'Ratio of SFA to total FA',
                            "Ratio of Cholesterol to total lipids in Large VLDL", 'Ratio of CE to total lipids in Large VLDL',
                            "Ratio of FC to total lipids in Medium VLDL")


# Visualizing the bipartite graph
fig_4_plot <- visIgraph(keto_regcor_G,
          idToLabel =F,
          layout = "layout_nicely",
          physics = FALSE,
          smooth = FALSE,
          type = "square",
          randomSeed = NULL,
          layoutMatrix = NULL )

visSave(fig_4_plot, file = "results/Figure_4.html")


## finished
print("Done!") 
print("Generated pdfs and html with plots in results folder!") 

